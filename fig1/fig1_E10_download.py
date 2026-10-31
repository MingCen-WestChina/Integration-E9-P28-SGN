# Source: download.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode, urlparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from collections import Counter
import csv
import gzip
import hashlib
import h5py
import io
import json
import re
import shutil
import threading
import time

# Set each storage placeholder below to its own local file or directory.
ROOT = Path()
# Set the file or directory for dest.
DEST = ROOT
# Set the file or directory for temp.
TEMP = ROOT
# Set the input or working directory used by this step.
SOURCE = ROOT
# Set the manifest file; the released public manifests are in metadata/e10.
MANIFEST = ROOT
# Set the file or directory for state.
STATE = ROOT
WORKERS = 4
ATTEMPTS = 8
TIMEOUT = 60
CHUNK = 1024 * 1024
GEO = 'GSE231513'
PROJECT = 'PRJNA966179'
lock = threading.Lock()
active = {}


def log(message):
    print(time.strftime('%H:%M:%S') + ' ' + message, flush=True)


def get_bytes(url):
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={'User-Agent': 'E10-FASTQ-download/1.0'}), timeout=TIMEOUT) as response:
                return response.read()
        except OSError:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


def h5_column(node):
    if isinstance(node, h5py.Group):
        category = node['categories']
        values = category.asstr()[:] if h5py.check_string_dtype(category.dtype) else category[:]
        return [None if code < 0 else str(values[code]) for code in node['codes'][:]]
    return node.asstr()[:].tolist()


def build_manifest():
    with h5py.File(SOURCE, 'r') as handle:
        obs = handle['obs']
        index_key = obs.attrs['_index']
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        ids = h5_column(obs['original_cell_id'])
        global_ids = h5_column(obs[index_key])
        plates = h5_column(obs['plate_id'])
        stages = h5_column(obs['stage_standardized'])
    if len(ids) != 73 or len(set(ids)) != 73 or set(stages) != {'E10.5'}:
        raise ValueError('Selected-cell membership differs from the approved 73 E10.5 cells.')
    selected = {cell: (global_id, plate) for cell, global_id, plate in zip(ids, global_ids, plates)}
    geo_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231513/soft/GSE231513_family.soft.gz'
    geo_text = gzip.decompress(get_bytes(geo_url)).decode('utf-8')
    geo_records = {}
    for block in re.split(r'(?m)^\^SAMPLE = ', geo_text)[1:]:
        lines = block.splitlines()
        gsm = lines[0].strip()
        fields = {}
        for line in lines[1:]:
            if line.startswith('^'):
                break
            if line.startswith('!Sample_') and ' = ' in line:
                key, value = line.split(' = ', 1)
                fields.setdefault(key, []).append(value.strip())
        cells = [x for x in fields.get('!Sample_description', []) if x in selected]
        if not cells:
            continue
        if len(cells) != 1 or cells[0] in geo_records or GEO not in fields.get('!Sample_series_id', []):
            raise ValueError('Ambiguous selected-cell GEO mapping: ' + gsm)
        relations = ' '.join(fields.get('!Sample_relation', []))
        srx = set(re.findall(r'\bSRX\d+\b', relations))
        biosample = set(re.findall(r'\bSAMN\d+\b', relations))
        if len(srx) != 1 or len(biosample) != 1:
            raise ValueError('Ambiguous experiment or BioSample: ' + gsm)
        cell = cells[0]
        if selected[cell][1] != GEO + '_plate_' + cell.rsplit('.', 1)[1]:
            raise ValueError('Plate identity conflict: ' + cell)
        geo_records[cell] = (gsm, next(iter(srx)), next(iter(biosample)))
    if set(geo_records) != set(selected):
        raise ValueError('GEO does not cover every selected cell.')
    fields = ['run_accession', 'experiment_accession', 'study_accession', 'sample_accession', 'sample_alias', 'scientific_name', 'library_strategy', 'library_layout', 'fastq_ftp', 'fastq_bytes', 'fastq_md5']
    ena_url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urlencode({'accession': PROJECT, 'result': 'read_run', 'fields': ','.join(fields), 'format': 'tsv', 'download': 'false'})
    ena_text = get_bytes(ena_url).decode('utf-8')
    experiments = {}
    for row in csv.DictReader(io.StringIO(ena_text), delimiter='\t'):
        experiments.setdefault(row['experiment_accession'], []).append(row)
    records = []
    for cell in ids:
        gsm, srx, biosample = geo_records[cell]
        runs = experiments.get(srx, [])
        if len(runs) != 1:
            raise ValueError('Run count differs from the verified single-run mapping: ' + cell)
        row = runs[0]
        if (row['study_accession'], row['sample_accession'], row['sample_alias'], row['scientific_name'], row['library_strategy'], row['library_layout']) != (PROJECT, biosample, gsm, 'Mus musculus', 'RNA-Seq', 'PAIRED'):
            raise ValueError('ENA identity or layout mismatch: ' + cell)
        urls, sizes, md5s = [row[key].split(';') for key in ['fastq_ftp', 'fastq_bytes', 'fastq_md5']]
        if not (len(urls) == len(sizes) == len(md5s) == 2):
            raise ValueError('Incomplete paired-end file metadata: ' + cell)
        mates = set()
        for raw_url, raw_size, md5 in zip(urls, sizes, md5s):
            url = raw_url if '://' in raw_url else 'https://' + raw_url
            parsed = urlparse(url)
            filename = Path(parsed.path).name
            match = re.fullmatch(re.escape(row['run_accession']) + r'_([12])\.fastq\.gz', filename)
            if parsed.scheme != 'https' or parsed.hostname != 'ftp.sra.ebi.ac.uk' or not match or not raw_size.isdigit() or not re.fullmatch(r'[0-9a-fA-F]{32}', md5):
                raise ValueError('Unexpected FASTQ URL or checksum: ' + url)
            mate = match[1]
            mates.add(mate)
            records.append({'cell_id_global': selected[cell][0], 'original_cell_id': cell, 'plate_id': selected[cell][1], 'GSM': gsm, 'experiment_accession': srx, 'run_accession': row['run_accession'], 'mate': mate, 'filename': filename, 'url': url, 'expected_bytes': int(raw_size), 'expected_md5': md5.lower(), 'actual_bytes': '', 'actual_md5': '', 'status': 'pending'})
        if mates != {'1', '2'}:
            raise ValueError('Invalid mate allocation: ' + cell)
    if len(records) != 146 or len({r['filename'] for r in records}) != 146 or len({r['run_accession'] for r in records}) != 73:
        raise ValueError('Final manifest does not match the approved 73 pairs.')
    return records


def digest(path):
    value = hashlib.md5()
    with path.open('rb') as handle:
        while data := handle.read(4 * CHUNK):
            value.update(data)
    return value.hexdigest()


def write_manifest(records):
    temporary = TEMP / 'manifest-writing.tsv'
    with temporary.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    temporary.replace(MANIFEST)


def transfer(row):
    target = DEST / row['filename']
    partial = TEMP / (row['filename'] + '.part')
    expected_size = row['expected_bytes']
    if target.exists():
        actual_md5 = digest(target)
        if target.stat().st_size != expected_size or actual_md5 != row['expected_md5']:
            raise ValueError('Existing destination does not match; left untouched: ' + str(target))
        return expected_size, actual_md5
    for attempt in range(1, ATTEMPTS + 1):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            if offset > expected_size:
                raise ValueError('Partial download is larger than the source metadata.')
            with lock:
                active[row['filename']] = offset
            if offset < expected_size:
                headers = {'User-Agent': 'E10-FASTQ-download/1.0', 'Accept-Encoding': 'identity'}
                if offset:
                    headers['Range'] = f'bytes={offset}-'
                with urlopen(Request(row['url'], headers=headers), timeout=TIMEOUT) as response:
                    if response.status == 206:
                        content_range = response.headers.get('Content-Range', '')
                        if not content_range.startswith(f'bytes {offset}-'):
                            raise ValueError('Unexpected resume offset: ' + content_range)
                        mode = 'ab' if offset else 'wb'
                    elif response.status == 200:
                        offset, mode = 0, 'wb'
                    else:
                        raise OSError('Unexpected HTTP status ' + str(response.status))
                    with partial.open(mode) as handle:
                        while data := response.read(CHUNK):
                            handle.write(data)
                            offset += len(data)
                            if offset > expected_size:
                                raise ValueError('Response exceeds the expected file size.')
                            with lock:
                                active[row['filename']] = offset
                if partial.stat().st_size != expected_size:
                    raise OSError('Transfer ended before reaching the expected size.')
            actual_md5 = digest(partial)
            if actual_md5 != row['expected_md5']:
                # Only this task's partial file may be reset after a checksum failure.
                partial.open('wb').close()
                raise OSError('MD5 mismatch; restarting this file from byte zero.')
            if target.exists():
                raise FileExistsError('Destination appeared during transfer; left untouched.')
            partial.rename(target)
            with lock:
                active.pop(row['filename'], None)
            return expected_size, actual_md5
        except (OSError, ValueError) as error:
            if attempt == ATTEMPTS:
                raise
            log(f'RETRY {row["filename"]} attempt={attempt} reason={type(error).__name__}: {error}')
            time.sleep(min(2 * attempt, 15))


def main():
    if not DEST.is_dir():
        raise FileNotFoundError('The verified destination directory is missing.')
    TEMP.mkdir(parents=True, exist_ok=True)
    if MANIFEST.exists():
        with MANIFEST.open(encoding='utf-8', newline='') as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
        for row in records:
            row['expected_bytes'] = int(row['expected_bytes'])
        if len(records) != 146 or len({r['original_cell_id'] for r in records}) != 73:
            raise ValueError('Existing manifest does not have the authorized membership.')
    else:
        records = build_manifest()
        write_manifest(records)
    total_bytes = sum(row['expected_bytes'] for row in records)
    if shutil.disk_usage(DEST).free < total_bytes + 1024**3:
        raise OSError('Insufficient free disk space.')
    log(f'START cells=73 files=146 bytes={total_bytes} workers={WORKERS} destination={DEST}')
    start = time.monotonic()
    completed = 0
    completed_bytes = 0
    failures = []
    next_status = start
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        pending = {pool.submit(transfer, row): row for row in records}
        while pending:
            done, _ = wait(pending, timeout=15, return_when=FIRST_COMPLETED)
            for future in done:
                row = pending.pop(future)
                try:
                    size, md5 = future.result()
                    row.update(actual_bytes=size, actual_md5=md5, status='verified')
                    completed += 1
                    completed_bytes += size
                    log(f'VERIFIED {completed}/146 {row["filename"]}')
                except Exception as error:
                    row['status'] = 'failed'
                    failures.append({'filename': row['filename'], 'error': str(error)})
                    with lock:
                        active.pop(row['filename'], None)
                    log(f'FAILED {row["filename"]}: {error}')
                write_manifest(records)
            now = time.monotonic()
            if now >= next_status or not pending:
                with lock:
                    partial_bytes = sum(active.values())
                state = {'verified_files': completed, 'total_files': len(records), 'verified_bytes': completed_bytes, 'active_partial_bytes': partial_bytes, 'total_bytes': total_bytes, 'elapsed_seconds': round(now-start, 1), 'failures': failures, 'running': bool(pending)}
                STATE.write_text(json.dumps(state, indent=2), encoding='utf-8')
                log(f'PROGRESS verified={completed}/146 received_GB={(completed_bytes+partial_bytes)/1e9:.3f}/{total_bytes/1e9:.3f} elapsed_s={now-start:.0f}')
                next_status = now + 30
    if failures:
        raise RuntimeError(f'{len(failures)} files remain incomplete; see state.json and manifest.')
    log('COMPLETE: all 146 files passed exact byte-size and MD5 checks. No read processing was performed.')


if __name__ == '__main__':
    main()
