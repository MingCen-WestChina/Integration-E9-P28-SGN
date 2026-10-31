#!/usr/bin/env python3
# Source: recover_quality.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
"""Recover the 73 E10.5 full-quality FASTQ pairs; stop after validation and FastQC.

Use --download-only for network transfer; use --process-only in the compute job.
The default performs both phases. No trimming or alignment is performed.
The script directory is WORK unless --work is supplied. All temporary files remain
under WORK. Only Normalized SRA files matching the manifest size and MD5 are used.
Sequence SHA256 uses ordered_uppercase_sequence_plus_LF_v1: each uppercase
sequence plus LF in spot order; headers and qualities are excluded. Every mate
must exactly match the independently validated ENA sequence baseline.
"""

# Imports and fixed limits
import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter

# Set each storage placeholder below to its own local file or directory.
ROOT = Path()

SCHEMA = "e10-full-quality-recovery-v1"
CHUNK = 4 * 1024 * 1024
NETWORK_WORKERS = 16
ACTIVE_RUNS = 4
PER_RUN_CHUNKS = 4
TIMEOUT = 60
ATTEMPTS = 6
STOP = threading.Event()
IDENTITY = ("cell_id_global", "original_cell_id", "plate_id", "GSM",
            "experiment_accession", "run_accession")


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def check_stop():
    if STOP.is_set():
        raise InterruptedError("Recovery interrupted; resumable files were retained.")


def digest(path, algorithm="md5"):
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            check_stop()
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".writing")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


class Store:
    def __init__(self, work, manifest_sha, baseline_sha):
        self.work = work
        self.path = work / "quality_recovery_state.json"
        # Set full_quality_results.json, written here and read by fig1_E10_trim_pilot.py and fig1_E10_full_counts.py.
        self.results_path = ROOT
        self.lock = threading.RLock()
        self.token = {"schema": SCHEMA, "manifest_sha256": manifest_sha,
                      "baseline_sha256": baseline_sha}
        self.data = dict(self.token, created_utc=now(), runs={})
        self.results = dict(self.token, expected_runs=73, runs={})
        for path, attribute in ((self.path, "data"), (self.results_path, "results")):
            if path.exists():
                value = json.loads(path.read_text())
                require(all(value.get(k) == v for k, v in self.token.items()),
                        f"Unrecognized state file: {path}")
                setattr(self, attribute, value)
        self.save()

    def save(self):
        with self.lock:
            self.data["updated_utc"] = now()
            atomic_json(self.path, self.data)

    def update(self, run, **values):
        with self.lock:
            self.data["runs"].setdefault(run, {}).update(values, updated_utc=now())
            self.save()

    def get(self, run):
        with self.lock:
            return dict(self.data["runs"].get(run, {}))

    def result(self, run):
        with self.lock:
            return self.results["runs"].get(run)

    def put_result(self, run, value):
        with self.lock:
            self.results["runs"][run] = value
            self.results["updated_utc"] = now()
            self.results["completed_runs"] = sum(
                x.get("fastqc_complete", False) for x in self.results["runs"].values())
            atomic_json(self.results_path, self.results)

    def owned_directory(self, path):
        require(not path.is_symlink(), f"Refusing symlink directory: {path}")
        require(path.resolve().is_relative_to(self.work) and path.resolve() != self.work,
                f"Managed directory is outside WORK: {path}")
        path.mkdir(parents=True, exist_ok=True)
        marker = path / "recovery_owner.json"
        if marker.exists():
            require(json.loads(marker.read_text()) == self.token,
                    f"Directory ownership mismatch: {path}")
        else:
            require(not any(path.iterdir()), f"Refusing unmanaged nonempty directory: {path}")
            with marker.open("x", encoding="utf-8") as handle:
                json.dump(self.token, handle)
        return path

    def clear_owned_directory(self, path):
        self.owned_directory(path)
        # Only this run's explicitly marked scratch directory is removed.
        shutil.rmtree(path)
        return self.owned_directory(path)


# Manifest and tool checks
def load_manifest(path):
    raw = path.read_bytes()
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t"))
    require(len(rows) == 73, "The manifest must contain exactly 73 runs.")
    for key in ("run_accession", "cell_id_global", "GSM", "experiment_accession"):
        require(len({r[key] for r in rows}) == 73, f"Nonunique {key} in manifest.")
    for row in rows:
        run = row["run_accession"]
        require(re.fullmatch(r"SRR[0-9]+", run), f"Invalid run: {run}")
        require(row["normal_semantic_name"] == "SRA Normalized", f"Not Normalized: {run}")
        parsed = urllib.parse.urlparse(row["normal_sra_url"])
        require(parsed.scheme == "https" and parsed.netloc == "sra-pub-run-odp.s3.amazonaws.com"
                and parsed.path == f"/sra/{run}/{run}" and not parsed.query,
                f"Unexpected Normalized source URL: {run}")
        require(re.fullmatch(r"[0-9a-f]{32}", row["normal_sra_md5"]), f"Invalid MD5: {run}")
        require(row["normal_sra_md5"] != row["lite_sra_md5"]
                and row["normal_sra_url"] != row["lite_sra_url"], f"Lite source selected: {run}")
        require(row["normal_access_type"] == "anonymous" and row["normal_free_egress"] == "worldwide",
                f"Nonanonymous or charged source: {run}")
        for key in ("normal_sra_bytes", "total_spots", "total_bases", "read1_length", "read2_length"):
            row[key] = int(row[key])
        require(row["normal_sra_bytes"] > 0 and row["total_spots"] > 0,
                f"Invalid source size/count: {run}")
        require(row["read1_length"] == row["read2_length"] == 75
                and row["total_bases"] == 150 * row["total_spots"], f"Unexpected read layout: {run}")
    return rows, hashlib.sha256(raw).hexdigest()


def load_baseline(path, rows):
    raw = path.read_bytes()
    entries = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t"))
    require(len(entries) == 146, "The ENA baseline must contain exactly 146 mates.")
    expected = {(r["run_accession"], mate): r for r in rows for mate in (1, 2)}
    result = {}
    for entry in entries:
        key = (entry["run_accession"], int(entry["mate"]))
        require(key in expected and key not in result, f"Unexpected or duplicate baseline mate: {key}")
        entry["records"], entry["bases"] = int(entry["records"]), int(entry["bases"])
        require(entry["records"] == expected[key]["total_spots"]
                and entry["bases"] == 75 * entry["records"], f"Baseline count mismatch: {key}")
        require(re.fullmatch(r"[0-9a-f]{64}", entry["sequence_sha256"]),
                f"Invalid baseline SHA256: {key}")
        require(entry.get("canonicalization") == "ordered_uppercase_sequence_plus_LF_v1",
                f"Unexpected baseline hash algorithm: {key}")
        for identity in ("cell_id_global", "original_cell_id"):
            if identity in entry:
                require(entry[identity] == expected[key][identity],
                        f"Baseline cell mapping mismatch: {key}/{identity}")
        result[key] = entry
    require(set(result) == set(expected), "The baseline is missing selected mates.")
    return result, hashlib.sha256(raw).hexdigest()


def tools_config(work, env_bin, fastqc_dir):
    environment = os.environ.copy()
    environment["PATH"] = str(env_bin) + os.pathsep + environment.get("PATH", "")
    environment.update(OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for key in ("TMP", "TEMP", "TMPDIR"):
        environment[key] = str(work / "export_tmp")
    executables = {}
    for name in ("fasterq-dump", "pigz", "java", "perl"):
        value = shutil.which(name, path=environment["PATH"])
        require(value is not None, f"Required executable not found: {name}")
        executables[name] = value
    help_text = subprocess.check_output([executables["fasterq-dump"], "--help"],
                                        env=environment, stderr=subprocess.STDOUT, text=True)
    for option in ("--seq-defline", "--qual-defline", "--split-files", "--outfile"):
        require(option in help_text, f"fasterq-dump lacks required option: {option}")
    require((fastqc_dir / "fastqc").is_file(), "The official FastQC launcher is missing.")
    qc = [executables["perl"], str(fastqc_dir / "fastqc"), "--java", executables["java"],
          "--dir", str(work / "export_tmp"), "--nogroup", "--quiet"]
    qc_version = subprocess.check_output(qc + ["--version"], env=environment,
                                         stderr=subprocess.STDOUT, text=True).strip()
    require("0.12.1" in qc_version, f"Expected FastQC 0.12.1; found {qc_version}")
    versions = {"FastQC": qc_version}
    for name in ("fasterq-dump", "pigz"):
        versions[name] = subprocess.check_output([executables[name], "--version"],
                                                env=environment, stderr=subprocess.STDOUT,
                                                text=True).strip()
    return executables, qc, environment, versions


# Bounded, resumable range downloads
def fetch_chunk(row, start, end, path):
    expected = end - start + 1
    require(not path.is_symlink(), f"Refusing symlink partial: {path}")
    require(not path.exists() or path.stat().st_size <= expected, f"Oversized partial: {path}")
    last_error = None
    for attempt in range(ATTEMPTS):
        check_stop()
        have = path.stat().st_size if path.exists() else 0
        if have == expected:
            return
        offset = start + have
        request = urllib.request.Request(row["normal_sra_url"], headers={
            "Range": f"bytes={offset}-{end}", "Accept-Encoding": "identity",
            "User-Agent": "E10FullQualityRecovery/1"})
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                total = row["normal_sra_bytes"]
                wanted_range = f"bytes {offset}-{end}/{total}"
                if response.status == 206:
                    require(response.headers.get("Content-Range") == wanted_range,
                            f"Unexpected Content-Range for {row['run_accession']}: {response.headers}")
                else:
                    # A 200 response is allowed only when the request is the entire object.
                    require(response.status == 200 and offset == 0 and end == total - 1,
                            f"Server did not honor the bounded Range request: HTTP {response.status}")
                require(int(response.headers.get("Content-Length", "-1")) == end - offset + 1,
                        "Unexpected HTTP Content-Length.")
                require(response.headers.get("Content-Encoding", "identity") == "identity",
                        "Unexpected HTTP content encoding.")
                remaining = end - offset + 1
                with path.open("ab") as handle:
                    while remaining:
                        check_stop()
                        block = response.read(min(64 * 1024, remaining))
                        if not block:
                            raise OSError("Incomplete HTTP body; retrying the remaining range.")
                        handle.write(block)
                        remaining -= len(block)
                    handle.flush()
                    os.fsync(handle.fileno())
            require(path.stat().st_size == expected, "Chunk length mismatch.")
            return
        except (urllib.error.URLError, OSError, RuntimeError) as error:
            last_error = error
            if STOP.is_set():
                raise
            if isinstance(error, urllib.error.HTTPError) and error.code in (400, 401, 403, 404):
                raise
            if attempt + 1 < ATTEMPTS:
                STOP.wait(min(30, 2 ** attempt))
    raise RuntimeError(f"Chunk failed after {ATTEMPTS} attempts: {path}: {last_error}")


def download_run(row, work, store, chunk_pool, allow_network=True):
    check_stop()
    run = row["run_accession"]
    final = work / "normal_sra" / f"{run}.sra"
    require(not final.is_symlink(), f"Refusing symlink final SRA: {final}")
    if final.exists():
        require(final.stat().st_size == row["normal_sra_bytes"]
                and digest(final) == row["normal_sra_md5"],
                f"Existing SRA is not the expected verified Normalized file; not overwriting: {final}")
        store.update(run, download="verified", sra_path=str(final))
        return final
    require(allow_network, f"Verified SRA is required for --process-only: {final}")
    pieces = store.owned_directory(work / "normal_sra_parts" / run)
    store.update(run, download="in_progress")
    bounds = [(start, min(start + CHUNK, row["normal_sra_bytes"]) - 1)
              for start in range(0, row["normal_sra_bytes"], CHUNK)]
    iterator = iter(bounds)
    pending = set()

    def schedule():
        try:
            start, end = next(iterator)
        except StopIteration:
            return False
        path = pieces / f"{start:012d}-{end:012d}.part"
        pending.add(chunk_pool.submit(fetch_chunk, row, start, end, path))
        return True

    for _ in range(PER_RUN_CHUNKS):
        schedule()
    try:
        while pending:
            check_stop()
            finished, _ = cf.wait(pending, timeout=1, return_when=cf.FIRST_COMPLETED)
            for future in finished:
                pending.remove(future)
                future.result()
                schedule()
    except BaseException:
        for future in pending:
            future.cancel()
        raise
    assembled = pieces / f"{run}.assembling"
    h = hashlib.md5()
    size = 0
    with assembled.open("wb") as output:
        for start, end in bounds:
            path = pieces / f"{start:012d}-{end:012d}.part"
            require(path.stat().st_size == end - start + 1, f"Incomplete chunk: {path}")
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    check_stop()
                    output.write(block)
                    h.update(block)
                    size += len(block)
        output.flush()
        os.fsync(output.fileno())
    if size != row["normal_sra_bytes"] or h.hexdigest() != row["normal_sra_md5"]:
        # Re-fetch this run on the next invocation; no invalid final file is published.
        store.clear_owned_directory(pieces)
        raise RuntimeError(f"Assembled Normalized SRA failed size/MD5: {run}")
    require(not final.exists(), f"Final SRA appeared unexpectedly: {final}")
    os.replace(assembled, final)
    store.update(run, download="verified", sra_path=str(final), sra_verified_utc=now())
    store.clear_owned_directory(pieces)
    return final


# Serial export, compression and full FASTQ validation
def command(argv, environment, cwd, log_path, stdout_path=None):
    check_stop()
    with log_path.open("ab") as log:
        log.write(("\nCOMMAND " + json.dumps(argv) + "\n").encode())
        output = stdout_path.open("xb") if stdout_path else log
        try:
            process = subprocess.Popen(argv, cwd=cwd, env=environment, stdout=output,
                                       stderr=log, start_new_session=True)
            while process.poll() is None:
                if STOP.wait(1):
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    check_stop()
            require(process.returncode == 0, f"Command failed ({process.returncode}); see {log_path}")
        finally:
            if stdout_path:
                output.close()


def validate_pair(paths, row, baseline):
    run = row["run_accession"]
    hashes = [hashlib.sha256(), hashlib.sha256()]
    quality = [Counter(), Counter()]
    n_bases = [0, 0]
    records = 0
    expected_spots = row["total_spots"]
    with gzip.open(paths[0], "rb") as r1, gzip.open(paths[1], "rb") as r2:
        while True:
            headers = [r1.readline(), r2.readline()]
            if not any(headers):
                break
            require(all(headers), f"Unequal paired FASTQ record counts: {run}")
            records += 1
            if records % 10000 == 0:
                check_stop()
            require(records <= expected_spots, f"Too many spots: {run}")
            for index, handle in enumerate((r1, r2)):
                sequence, plus, scores = handle.readline(), handle.readline(), handle.readline()
                require(sequence and plus and scores, f"Truncated FASTQ record: {run}/{records}")
                header = headers[index].rstrip(b"\r\n")
                require(header == f"@{run}.{records}/{index + 1}".encode(),
                        f"Unexpected spot/mate identifier at {run}/{records}: {header!r}")
                sequence, scores = sequence.rstrip(b"\r\n"), scores.rstrip(b"\r\n")
                require(plus.rstrip(b"\r\n") == b"+", f"Invalid separator: {run}/{records}")
                require(len(sequence) == len(scores) == 75, f"Non-75bp read: {run}/{records}")
                require(re.fullmatch(b"[ACGTN]{75}", sequence) is not None,
                        f"Invalid sequence alphabet: {run}/{records}")
                require(all(33 <= q <= 126 for q in scores), f"Invalid Phred33 quality: {run}/{records}")
                hashes[index].update(sequence.upper() + b"\n")
                quality[index].update(scores)
                n_bases[index] += sequence.count(b"N")
    require(records == expected_spots, f"Missing spots: {run}: {records} != {expected_spots}")
    require(all(len(q) > 1 for q in quality), f"Constant quality remains in at least one mate: {run}")
    mates = []
    for index, path in enumerate(paths):
        expected = baseline[(run, index + 1)]
        require(records == expected["records"] and records * 75 == expected["bases"]
                and hashes[index].hexdigest() == expected["sequence_sha256"],
                f"Full-quality sequences differ from the ENA baseline: {run} mate {index + 1}")
        histogram = {str(q - 33): count for q, count in sorted(quality[index].items())}
        count = records * 75
        mates.append({"mate": index + 1, "filename": path.name, "records": records,
                      "read_length": 75, "bases": count, "n_bases": n_bases[index],
                      "quality_histogram_phred33": histogram,
                      "quality_min": min(quality[index]) - 33,
                      "quality_max": max(quality[index]) - 33,
                      "quality_mean": sum((q - 33) * n for q, n in quality[index].items()) / count,
                      "q20_fraction": sum(n for q, n in quality[index].items() if q >= 53) / count,
                      "q30_fraction": sum(n for q, n in quality[index].items() if q >= 63) / count,
                      "sequence_sha256": hashes[index].hexdigest(),
                      "sequence_sha256_algorithm": "ordered_uppercase_sequence_plus_LF_v1",
                      "ena_sequence_match": True,
                      "gzip_bytes": path.stat().st_size, "gzip_md5": digest(path)})
    return mates


def qc_outputs(directory, run):
    return [directory / f"{run}_{mate}_fastqc.{extension}"
            for mate in (1, 2) for extension in ("html", "zip")]


def validate_fastqc(directory, run):
    output = []
    for path in qc_outputs(directory, run):
        require(path.is_file() and path.stat().st_size > 0, f"Missing FastQC output: {path}")
        item = {"filename": path.name, "bytes": path.stat().st_size, "md5": digest(path)}
        if path.suffix == ".zip":
            with zipfile.ZipFile(path) as archive:
                require(archive.testzip() is None, f"Corrupt FastQC ZIP: {path}")
                prefix = path.stem + "/"
                item["module_summary"] = archive.read(prefix + "summary.txt").decode().splitlines()
                data = archive.read(prefix + "fastqc_data.txt").decode()
                require(data.startswith("##FastQC\t0.12.1"), f"Unexpected FastQC version: {path}")
        output.append(item)
    return output


def postprocess(row, sra, work, store, executables, qc_command, environment, versions, baseline):
    check_stop()
    run = row["run_accession"]
    finals = [work / "fastq_full" / f"{run}_{mate}.fastq.gz" for mate in (1, 2)]
    prior = store.result(run)
    state = store.get(run)
    if prior and all(p.exists() for p in finals):
        require(all(p.stat().st_size == m["gzip_bytes"] and digest(p) == m["gzip_md5"]
                    for p, m in zip(finals, prior["mates"])),
                f"Existing verified FASTQ was changed; not overwriting: {run}")
        if prior.get("fastqc_complete"):
            current = validate_fastqc(work / "fastqc_full", run)
            require([(x["filename"], x["md5"]) for x in current] ==
                    [(x["filename"], x["md5"]) for x in prior["fastqc_outputs"]],
                    f"Existing FastQC output was changed: {run}")
            store.update(run, postprocess="complete")
            return
    else:
        # Publication intent is saved before atomic moves, allowing crash recovery.
        candidate = state.get("publication_candidate") or prior
        if candidate and all(p.exists() for p in finals):
            require(all(p.stat().st_size == m["gzip_bytes"] and digest(p) == m["gzip_md5"]
                        for p, m in zip(finals, candidate["mates"])),
                    f"Interrupted publication does not match validated FASTQ: {run}")
            prior = candidate
            store.put_result(run, prior)
        else:
            if any(p.exists() for p in finals):
                require(candidate is not None, f"Refusing unknown existing final FASTQ: {run}")
                for path, metric in zip(finals, candidate["mates"]):
                    if path.exists():
                        require(path.stat().st_size == metric["gzip_bytes"] and digest(path) == metric["gzip_md5"],
                                f"Refusing changed final FASTQ during recovery: {path}")
                        path.unlink()
            scratch = store.clear_owned_directory(work / "export_tmp" / run)
            temp = scratch / "scratch"
            temp.mkdir()
            require(shutil.disk_usage(work).free > max(2 * 1024 ** 3, 12 * row["normal_sra_bytes"]),
                    f"Insufficient working disk space for export: {run}")
            store.update(run, postprocess="exporting")
            log = work / "recovery_logs" / f"{run}.log"
            argv = [executables["fasterq-dump"], str(sra), "--split-files", "--threads", "4",
                    "--temp", str(temp), "--outdir", str(scratch), "--outfile", f"{run}.fastq",
                    "--seq-defline", f"@{run}.$si/$ri", "--qual-defline", "+"]
            command(argv, environment, scratch, log)
            raw = [scratch / f"{run}_{mate}.fastq" for mate in (1, 2)]
            require(set(scratch.glob("*.fastq")) == set(raw) and all(p.is_file() for p in raw),
                    f"Unexpected fasterq-dump output set: {run}")
            compressed = [p.with_suffix(p.suffix + ".gz") for p in raw]
            store.update(run, postprocess="compressing")
            for path, target in zip(raw, compressed):
                command([executables["pigz"], "-p", "4", "-n", "-c", str(path)],
                        environment, scratch, log, stdout_path=target)
            store.update(run, postprocess="validating")
            metrics = validate_pair(compressed, row, baseline)
            prior = {key: row[key] for key in IDENTITY}
            prior.update(normal_sra_url=row["normal_sra_url"], normal_sra_md5=row["normal_sra_md5"],
                         normal_sra_bytes=row["normal_sra_bytes"], expected_spots=row["total_spots"],
                         validated_utc=now(), paired_ids_verified=True, unfiltered=True,
                         sequence_sha256_convention="ordered_uppercase_sequence_plus_LF_v1",
                         mates=metrics, tool_versions=versions, fastqc_complete=False)
            store.update(run, publication_candidate=prior, postprocess="publishing")
            for source, final in zip(compressed, finals):
                require(not final.exists() and not final.is_symlink(), f"Final FASTQ already exists: {final}")
                os.replace(source, final)
            store.put_result(run, prior)
            store.clear_owned_directory(scratch)
    store.update(run, postprocess="fastqc")
    qc_scratch = store.clear_owned_directory(work / "export_tmp" / (run + "_fastqc"))
    log = work / "recovery_logs" / f"{run}.log"
    command(qc_command + ["--threads", "2", "--noextract", "--outdir", str(qc_scratch)] +
            [str(p) for p in finals], environment, qc_scratch, log)
    qc_metrics = validate_fastqc(qc_scratch, run)
    store.update(run, qc_publication_candidate=qc_metrics)
    old_qc = state.get("qc_publication_candidate")
    old_md5 = {m["filename"]: m["md5"] for m in old_qc or []}
    for source in qc_outputs(qc_scratch, run):
        destination = work / "fastqc_full" / source.name
        if destination.exists():
            require(source.name in old_md5 and digest(destination) == old_md5[source.name],
                    f"Refusing unknown existing FastQC output: {destination}")
        os.replace(source, destination)
    prior.update(fastqc_complete=True, fastqc_outputs=qc_metrics, completed_utc=now())
    store.put_result(run, prior)
    store.update(run, postprocess="complete", publication_candidate=None, qc_publication_candidate=None)
    store.clear_owned_directory(qc_scratch)


# Orchestration and resumable termination
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # Set the local work/input/tool location, or supply the matching command-line option.
    parser.add_argument("--work", type=Path, default=ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--env-bin", type=Path)
    parser.add_argument("--fastqc-dir", type=Path)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--download-only", action="store_true", help="Download and verify 73 SRA files, then stop.")
    modes.add_argument("--process-only", action="store_true", help="Require local verified SRA files; never download.")
    arguments = parser.parse_args()
    require(sys.platform.startswith("linux"), "This recovery script must run on Linux.")
    work = arguments.work.resolve()
    require(work.is_dir(), f"WORK does not exist: {work}")
    # Set the local work/input/tool location, or supply the matching command-line option.
    manifest = arguments.manifest or ROOT
    rows, manifest_sha = load_manifest(manifest)
    # Set the local work/input/tool location, or supply the matching command-line option.
    baseline, baseline_sha = load_baseline(arguments.baseline or ROOT, rows)
    # All subsequently created threads and subprocesses inherit at most four CPUs.
    allowed_cpus = sorted(os.sched_getaffinity(0))
    require(allowed_cpus, "No CPUs are available in the job affinity mask.")
    os.sched_setaffinity(0, allowed_cpus[:4])
    import fcntl
    lockfile = (work / "quality_recovery.lock").open("a+")
    try:
        fcntl.flock(lockfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("Another recovery process already holds this WORK lock.")
    for event in (signal.SIGINT, signal.SIGTERM):
        signal.signal(event, lambda *_: STOP.set())
    store = Store(work, manifest_sha, baseline_sha)
    directories = ["normal_sra_parts", "normal_sra"]
    if not arguments.download_only:
        directories += ["export_tmp", "fastq_full", "fastqc_full", "recovery_logs"]
    for name in directories:
        store.owned_directory(work / name)
    executables, qc, environment, versions = {}, [], {}, {}
    if arguments.process_only:
        for row in rows:
            path = work / "normal_sra" / (row["run_accession"] + ".sra")
            require(path.is_file() and not path.is_symlink() and path.stat().st_size == row["normal_sra_bytes"],
                    f"--process-only requires all 73 correctly sized local SRA files: {path}")
    if not arguments.download_only:
        executables, qc, environment, versions = tools_config(
            # Set the local work/input/tool location, or supply the matching command-line option.
            work, (arguments.env_bin or ROOT).resolve(),
            # Set the local work/input/tool location, or supply the matching command-line option.
            (arguments.fastqc_dir or ROOT).resolve())
    mode = "download_only" if arguments.download_only else "process_only" if arguments.process_only else "all"
    store.data["configuration"] = {"network_workers": NETWORK_WORKERS, "active_runs": ACTIVE_RUNS,
                                    "chunks_per_run": PER_RUN_CHUNKS, "chunk_bytes": CHUNK,
                                    "request_timeout_seconds": TIMEOUT, "attempts_per_chunk": ATTEMPTS,
                                    "postprocess_workers": 1, "export_threads": 4, "pigz_threads": 4,
                                    "fastqc_threads": 2, "tool_versions": versions, "mode": mode}
    store.save()
    print(f"Recovering 73 runs; mode={mode}; network limit=16; serial export/QC worker=1.", flush=True)
    with cf.ThreadPoolExecutor(max_workers=NETWORK_WORKERS) as chunks, \
            cf.ThreadPoolExecutor(max_workers=ACTIVE_RUNS) as downloads, \
            cf.ThreadPoolExecutor(max_workers=1) as processor:
        waiting = {downloads.submit(download_run, row, work, store, chunks,
                                    not arguments.process_only): ("download", row)
                   for row in rows}
        last_report = 0.0
        try:
            while waiting:
                check_stop()
                finished, _ = cf.wait(waiting, timeout=1, return_when=cf.FIRST_COMPLETED)
                for future in finished:
                    phase, row = waiting.pop(future)
                    result = future.result()
                    if phase == "download" and not arguments.download_only:
                        task = processor.submit(postprocess, row, result, work, store,
                                                executables, qc, environment, versions, baseline)
                        waiting[task] = ("postprocess", row)
                if time.monotonic() - last_report >= 30:
                    with store.lock:
                        verified = sum(r.get("download") == "verified" for r in store.data["runs"].values())
                        complete = sum(r.get("postprocess") == "complete" for r in store.data["runs"].values())
                    partial_bytes = 0
                    for path in (work / "normal_sra_parts").glob("*/*.part"):
                        try:
                            partial_bytes += path.stat().st_size
                        except FileNotFoundError:
                            pass
                    print(f"SRA verified {verified}/73; FASTQ + QC complete {complete}/73; "
                          f"active partial data {partial_bytes / 1024 ** 2:.1f} MiB.", flush=True)
                    with store.lock:
                        store.data["progress"] = {"sra_verified": verified, "fastq_qc_complete": complete,
                                                   "partial_bytes": partial_bytes}
                        store.save()
                    last_report = time.monotonic()
        except BaseException as error:
            STOP.set()
            for future in waiting:
                future.cancel()
            with store.lock:
                store.data["last_error"] = {"utc": now(), "message": str(error)}
                store.save()
            raise
    if arguments.download_only:
        require(all(store.get(r["run_accession"]).get("download") == "verified" for r in rows),
                "Download phase ended without 73 verified Normalized SRA files.")
        store.data.pop("last_error", None)
        store.data["downloads_completed_utc"] = now()
        store.save()
        print("Complete: 73 Normalized SRA files verified. Export and FastQC were not run.", flush=True)
        return
    require(len(store.results["runs"]) == 73 and all(r.get("fastqc_complete")
            for r in store.results["runs"].values()), "Recovery ended without 73 validated pairs.")
    store.data.pop("last_error", None)
    store.data["completed_utc"] = now()
    store.save()
    print("Complete: 73 full-quality pairs and FastQC reports. No trimming or alignment performed.", flush=True)


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        print(f"Recovery stopped: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        sys.exit(1)
