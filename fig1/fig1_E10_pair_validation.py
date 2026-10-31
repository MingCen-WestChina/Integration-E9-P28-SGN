# Source: pair_validation.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
import csv
import gzip
import h5py
import json
import numpy as np
import re
import time

# Parameters
# Set each storage placeholder below to its own local file or directory.
ROOT = Path()
# Set the file or directory for raw.
RAW = ROOT
# Set the file or directory for temp.
TEMP = ROOT
# Set the manifest file; the released public manifests are in metadata/e10.
MANIFEST = ROOT
# Set the approved 73-cell E10.5 selection H5AD input.
SELECTED = ROOT
# Set the file or directory for results.
RESULTS = ROOT
# Set the file or directory for status.
STATUS = ROOT
WORKERS = 2
POLL_SECONDS = 15
QUALITY_BATCH_RECORDS = 10000
MAX_IDLE_SECONDS = 7200
BASES = re.compile(rb'[ACGTN]+')
QUALITY = re.compile(rb'[!-~]+')


def log(message):
    print(time.strftime('%H:%M:%S') + ' ' + message, flush=True)


def atomic_json(path, value):
    temporary = path.with_name(path.name + '.writing')
    temporary.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temporary.replace(path)


def column(node):
    if isinstance(node, h5py.Group):
        categories = node['categories'].asstr()[:]
        return [None if code < 0 else str(categories[code]) for code in node['codes'][:]]
    return node.asstr()[:].tolist()


# Approved membership and download state
def selected_cells():
    with h5py.File(SELECTED, 'r') as handle:
        obs = handle['obs']
        index = obs.attrs['_index']
        if isinstance(index, bytes):
            index = index.decode()
        ids = column(obs['original_cell_id'])
        globals_ = column(obs[index])
        plates = column(obs['plate_id'])
        stages = column(obs['stage_standardized'])
    if len(ids) != 73 or len(set(ids)) != 73 or set(stages) != {'E10.5'}:
        raise ValueError('The selected object does not contain the approved 73 unique E10.5 cells.')
    return {cell: (global_id, plate) for cell, global_id, plate in zip(ids, globals_, plates)}


def manifest_pairs(selected):
    with MANIFEST.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 146 or {row['original_cell_id'] for row in rows} != set(selected):
        raise ValueError('Download manifest membership differs from the selected cells.')
    if len({row['filename'] for row in rows}) != 146 or len({row['run_accession'] for row in rows}) != 73:
        raise ValueError('Download manifest contains duplicate filenames or runs.')
    pairs = {}
    for row in rows:
        cell = row['original_cell_id']
        if (row['cell_id_global'], row['plate_id']) != selected[cell]:
            raise ValueError('Cell identity mismatch: ' + cell)
        if row['mate'] not in {'1', '2'}:
            raise ValueError('Invalid mate: ' + cell)
        filename = row['run_accession'] + '_' + row['mate'] + '.fastq.gz'
        if row['filename'] != filename or (RAW / filename).resolve().parent != RAW.resolve():
            raise ValueError('Invalid destination filename: ' + cell)
        if not re.fullmatch(r'[0-9a-f]{32}', row['expected_md5']) or int(row['expected_bytes']) <= 0:
            raise ValueError('Invalid source checksum or size: ' + cell)
        if row['mate'] in pairs.setdefault(cell, {}):
            raise ValueError('Duplicate mate: ' + cell)
        pairs[cell][row['mate']] = row
    for cell, pair in pairs.items():
        if set(pair) != {'1', '2'}:
            raise ValueError('Incomplete pair: ' + cell)
        for field in ['GSM', 'experiment_accession', 'run_accession']:
            if pair['1'][field] != pair['2'][field]:
                raise ValueError('Conflicting paired identity: ' + cell)
    return pairs


def is_ready(pair):
    for row in pair.values():
        if row['status'] != 'verified':
            return False
        expected = int(row['expected_bytes'])
        if int(row['actual_bytes']) != expected or row['actual_md5'] != row['expected_md5']:
            raise ValueError('A verified manifest row has inconsistent verification evidence.')
        path = RAW / row['filename']
        if not path.is_file() or path.stat().st_size != expected:
            raise ValueError('A verified FASTQ is missing or has changed size: ' + str(path))
    return True


# Full paired FASTQ validation
def normalized_header(header, mate, run, record):
    if not header.startswith(b'@'):
        raise ValueError(f'Invalid FASTQ header at record {record}, mate {mate}.')
    tokens = header[1:].split()
    if not tokens:
        raise ValueError(f'Empty FASTQ identifier at record {record}.')
    normalized = []
    for token in tokens:
        if token.endswith((b'/1', b'/2')):
            if token[-1:] != str(mate).encode():
                raise ValueError(f'Incorrect mate suffix at record {record}, mate {mate}.')
            token = token[:-2]
        elif len(token) > 2 and token[:2] in (b'1:', b'2:'):
            if token[:1] != str(mate).encode():
                raise ValueError(f'Incorrect CASAVA mate at record {record}, mate {mate}.')
            token = token[2:]
        normalized.append(token)
    if normalized[0].startswith(b'SRR') and not normalized[0].startswith(run.encode() + b'.'):
        raise ValueError(f'Run accession mismatch at record {record}.')
    return tuple(normalized)


def validate_pair(cell, pair):
    rows = [pair['1'], pair['2']]
    run = rows[0]['run_accession']
    paths = [RAW / row['filename'] for row in rows]
    before = [(p.stat().st_size, p.stat().st_mtime_ns) for p in paths]
    lengths = [Counter(), Counter()]
    quality_hist = [np.zeros(94, dtype=np.int64), np.zeros(94, dtype=np.int64)]
    quality_buffer = [bytearray(), bytearray()]
    n_bases = [0, 0]
    records = 0
    first_id = last_id = None
    with gzip.open(paths[0], 'rb') as left, gzip.open(paths[1], 'rb') as right:
        streams = [left, right]
        while True:
            headers = [handle.readline() for handle in streams]
            if not any(headers):
                break
            if not all(headers):
                raise ValueError(f'Unequal mate record counts after {records} records.')
            records += 1
            identifiers = []
            for index, handle in enumerate(streams):
                header = headers[index].rstrip(b'\r\n')
                sequence_line, plus_line, quality_line = handle.readline(), handle.readline(), handle.readline()
                if not sequence_line or not plus_line or not quality_line:
                    raise ValueError(f'Truncated FASTQ at record {records}, mate {index + 1}.')
                sequence = sequence_line.rstrip(b'\r\n')
                plus = plus_line.rstrip(b'\r\n')
                quality = quality_line.rstrip(b'\r\n')
                if not plus.startswith(b'+'):
                    raise ValueError(f'Missing plus line at record {records}, mate {index + 1}.')
                if plus[1:] and plus[1:].split()[0] != header[1:].split()[0]:
                    raise ValueError(f'Conflicting plus-line identifier at record {records}.')
                if not BASES.fullmatch(sequence):
                    raise ValueError(f'Empty sequence or a non-ACGTN base at record {records}, mate {index + 1}.')
                if len(sequence) != len(quality) or not QUALITY.fullmatch(quality):
                    raise ValueError(f'Invalid Phred33 characters or sequence/quality length mismatch at record {records}, mate {index + 1}.')
                identifiers.append(normalized_header(header, index + 1, run, records))
                lengths[index][len(sequence)] += 1
                n_bases[index] += sequence.count(b'N')
                quality_buffer[index].extend(quality)
            if identifiers[0] != identifiers[1]:
                raise ValueError(f'Paired record identifiers differ at record {records}.')
            first_id = first_id or identifiers[0][0].decode('ascii')
            last_id = identifiers[0][0].decode('ascii')
            if records % QUALITY_BATCH_RECORDS == 0:
                for index in (0, 1):
                    values = np.frombuffer(quality_buffer[index], dtype=np.uint8)
                    quality_hist[index] += np.bincount(values - 33, minlength=94)
                    del values
                    quality_buffer[index].clear()
    if records == 0:
        raise ValueError('Empty paired FASTQ files: ' + cell)
    for index in (0, 1):
        if quality_buffer[index]:
            values = np.frombuffer(quality_buffer[index], dtype=np.uint8)
            quality_hist[index] += np.bincount(values - 33, minlength=94)
    after = [(p.stat().st_size, p.stat().st_mtime_ns) for p in paths]
    if before != after:
        raise ValueError('FASTQ input changed during validation: ' + cell)
    mates = []
    for index in (0, 1):
        bases = sum(length * count for length, count in lengths[index].items())
        hist = {str(q): int(count) for q, count in enumerate(quality_hist[index]) if count}
        if sum(hist.values()) != bases:
            raise ValueError('Quality histogram does not cover every base.')
        mates.append({'mate': index + 1, 'filename': paths[index].name, 'records': records,
                      'bases': bases, 'read_lengths': dict(sorted(lengths[index].items())),
                      'n_bases': n_bases[index], 'q30_bases': int(quality_hist[index][30:].sum()),
                      'q30_percent': round(100 * float(quality_hist[index][30:].sum()) / bases, 4),
                      'phred33_histogram': hist})
    return {'original_cell_id': cell, 'cell_id_global': rows[0]['cell_id_global'],
            'plate_id': rows[0]['plate_id'], 'run_accession': run, 'status': 'passed',
            'paired_records': records, 'gzip_crc': 'passed_both_full_streams',
            'first_record_id': first_id, 'last_record_id': last_id, 'mates': mates}


# Incremental execution and concise status
def summary(results, failures, running, ready):
    passed = [result for result in results.values() if result['status'] == 'passed']
    mates = [mate for result in passed for mate in result['mates']]
    bases = sum(mate['bases'] for mate in mates)
    q30 = sum(mate['q30_bases'] for mate in mates)
    return {'checked_pairs': len(passed), 'total_pairs': 73, 'ready_pairs': ready,
            'paired_records': sum(result['paired_records'] for result in passed),
            'total_reads': 2 * sum(result['paired_records'] for result in passed),
            'total_bases': bases, 'q30_bases': q30,
            'q30_percent': round(100 * q30 / bases, 4) if bases else None,
            'read_lengths': sorted({int(length) for mate in mates for length in mate['read_lengths']}),
            'failures': failures, 'running': running,
            'note': 'All-record structure and gzip CRC validation; download MD5 is not recomputed. Quality values describe the deposited FASTQ encoding.'}


def main():
    selected = selected_cells()
    TEMP.mkdir(parents=True, exist_ok=True)
    if RESULTS.exists() or STATUS.exists():
        raise FileExistsError('Validation output already exists; inspect it before restarting.')
    results, failures, scheduled, pending = {}, {}, set(), {}
    last_progress = time.monotonic()
    previous_ready = -1
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        while len(results) + len(failures) < 73:
            pairs = manifest_pairs(selected)
            ready = {cell for cell, pair in pairs.items() if is_ready(pair)}
            if len(ready) != previous_ready:
                last_progress, previous_ready = time.monotonic(), len(ready)
            for cell in sorted(ready - scheduled):
                pending[pool.submit(validate_pair, cell, pairs[cell])] = cell
                scheduled.add(cell)
            if pending:
                done, _ = wait(pending, timeout=POLL_SECONDS, return_when=FIRST_COMPLETED)
                for future in done:
                    cell = pending.pop(future)
                    try:
                        results[cell] = future.result()
                        log(f'PASSED {len(results)}/73 cell={cell} paired_records={results[cell]["paired_records"]}')
                    except Exception as error:
                        failures[cell] = f'{type(error).__name__}: {error}'
                        log(f'FAILED cell={cell} {failures[cell]}')
                    last_progress = time.monotonic()
                if done:
                    atomic_json(RESULTS, {'cells': results, 'failures': failures})
            else:
                time.sleep(POLL_SECONDS)
            state = summary(results, failures, True, len(ready))
            atomic_json(STATUS, state)
            if not pending and time.monotonic() - last_progress > MAX_IDLE_SECONDS:
                raise TimeoutError('No new verified pairs became available within the idle limit.')
    state = summary(results, failures, False, len(ready))
    atomic_json(STATUS, state)
    atomic_json(RESULTS, {'cells': results, 'failures': failures})
    if failures:
        raise RuntimeError(f'{len(failures)} pair(s) failed full validation.')
    log('COMPLETE ' + json.dumps(state, separators=(',', ':')))


if __name__ == '__main__':
    main()
