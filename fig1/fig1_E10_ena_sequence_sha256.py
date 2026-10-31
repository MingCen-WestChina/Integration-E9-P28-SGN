# Source: ena_sequence_sha256.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import csv
import gzip
import hashlib
import json
import time

# Parameters
# Set each storage placeholder below to its own local file or directory.
ROOT = Path()
# Set the file or directory for raw.
RAW = ROOT
# Set the file or directory for temp.
TEMP = ROOT
# Set the file or directory for output.
OUTPUT = ROOT
WORKERS = 4
CANONICAL = 'ordered_uppercase_sequence_plus_LF_v1'
FIELDS = ['run_accession', 'mate', 'original_cell_id', 'cell_id_global', 'filename',
          'records', 'bases', 'sequence_sha256', 'normal_sra_total_spots',
          'matches_normal_spots', 'canonicalization']


def table(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


# Canonical sequence digest
def sequence_digest(row, expected_records, expected_bases, normal_spots):
    path = RAW / row['filename']
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    if before[0] != int(row['expected_bytes']):
        raise ValueError('Input size changed: ' + path.name)
    digest = hashlib.sha256()
    records = bases = 0
    with gzip.open(path, 'rb') as handle:
        while header := handle.readline():
            sequence_line, plus, quality_line = handle.readline(), handle.readline(), handle.readline()
            if not header.startswith(b'@') or not plus.startswith(b'+') or not sequence_line or not quality_line:
                raise ValueError('FASTQ framing changed: ' + path.name)
            sequence = sequence_line.rstrip(b'\r\n').upper()
            if len(sequence) != len(quality_line.rstrip(b'\r\n')):
                raise ValueError('FASTQ length mismatch: ' + path.name)
            digest.update(sequence + b'\n')
            records += 1
            bases += len(sequence)
    after = (path.stat().st_size, path.stat().st_mtime_ns)
    if before != after:
        raise ValueError('Input changed during hashing: ' + path.name)
    if records != expected_records or bases != expected_bases:
        raise ValueError('Counts differ from full paired validation: ' + path.name)
    if normal_spots is not None and records != normal_spots:
        raise ValueError('Read count differs from official normal SRA spots: ' + path.name)
    return {'run_accession': row['run_accession'], 'mate': row['mate'],
            'original_cell_id': row['original_cell_id'], 'cell_id_global': row['cell_id_global'],
            'filename': row['filename'], 'records': records, 'bases': bases,
            'sequence_sha256': digest.hexdigest(),
            'normal_sra_total_spots': normal_spots if normal_spots is not None else '',
            'matches_normal_spots': True if normal_spots is not None else '',
            'canonicalization': CANONICAL}


# Verified input contracts and execution
def main():
    if OUTPUT.exists():
        raise FileExistsError('Sequence digest output already exists; inspect before replacing.')
    # Uses download_manifest.tsv from fig1_E10_download.py.
    rows = table(ROOT)
    # Set the intermediate JSON file produced by the preceding E10 validation/recovery step.
    validation = json.loads((ROOT).read_text(encoding='utf-8'))
    cells = validation['cells']
    if validation['failures'] or len(cells) != 73 or len(rows) != 146:
        raise ValueError('The complete 73-pair validation is required.')
    if {row['original_cell_id'] for row in rows} != set(cells):
        raise ValueError('Cell membership differs from paired validation.')
    if len({row['filename'] for row in rows}) != 146 or len({(row['run_accession'], row['mate']) for row in rows}) != 146:
        raise ValueError('Duplicate FASTQ or run/mate keys.')
    # Set the file or directory for normal path.
    normal_path = ROOT
    normal = {}
    if normal_path.is_file():
        normal_rows = table(normal_path)
        normal = {row['run_accession']: row for row in normal_rows}
        if len(normal_rows) != 73 or len(normal) != 73 or set(normal) != {row['run_accession'] for row in rows}:
            raise ValueError('Normal SRA manifest run membership differs.')
    jobs = []
    for row in rows:
        cell = cells[row['original_cell_id']]
        mate = int(row['mate'])
        expected = next(item for item in cell['mates'] if item['mate'] == mate)
        if (row['status'] != 'verified' or row['actual_md5'] != row['expected_md5']
                or int(row['actual_bytes']) != int(row['expected_bytes'])):
            raise ValueError('Unverified source FASTQ: ' + row['filename'])
        if (row['run_accession'] != cell['run_accession'] or row['cell_id_global'] != cell['cell_id_global']
                or row['filename'] != expected['filename'] or row['filename'] != row['run_accession'] + '_' + row['mate'] + '.fastq.gz'):
            raise ValueError('Run, cell or mate identity differs from paired validation.')
        normal_spots = None
        if normal:
            source = normal[row['run_accession']]
            if (source['original_cell_id'], source['cell_id_global']) != (row['original_cell_id'], row['cell_id_global']):
                raise ValueError('Normal SRA cell identity differs.')
            normal_spots = int(source['total_spots'])
            if normal_spots != expected['records'] or int(source['total_bases']) != sum(item['bases'] for item in cell['mates']):
                raise ValueError('Normal SRA metadata read counts or base totals differ.')
        jobs.append((row, expected['records'], expected['bases'], normal_spots))
    print(f'START files={len(jobs)} workers={WORKERS} canonicalization={CANONICAL} normal_spots_checked={bool(normal)}', flush=True)
    results = {}
    start = time.monotonic()
    with OUTPUT.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        handle.flush()
        with ProcessPoolExecutor(max_workers=WORKERS) as pool:
            pending = {pool.submit(sequence_digest, *job): job[0]['filename'] for job in jobs}
            for future in as_completed(pending):
                result = future.result()
                results[result['filename']] = result
                writer.writerow(result)
                handle.flush()
                if len(results) % 10 == 0 or len(results) == len(jobs):
                    print(f'HASHED {len(results)}/{len(jobs)} elapsed_s={time.monotonic()-start:.0f}', flush=True)
    if len(results) != 146:
        raise ValueError('Incomplete sequence digest output.')
    written = table(OUTPUT)
    if len(written) != 146 or {row['filename'] for row in written} != set(results):
        raise ValueError('Saved sequence digest membership verification failed.')
    print('COMPLETE ' + json.dumps({'files': len(written), 'pairs': 73,
          'total_records_both_mates': sum(int(row['records']) for row in written),
          'total_bases': sum(int(row['bases']) for row in written),
          'normal_spots_checked': bool(normal), 'output': str(OUTPUT)}), flush=True)


if __name__ == '__main__':
    main()
