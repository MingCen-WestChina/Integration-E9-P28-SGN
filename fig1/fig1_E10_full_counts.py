#!/usr/bin/env python3
# Source: full_counts.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
"""Count the fixed 73-cell E10.5 cohort, reusing nine verified Q20 pilot results."""

# Fixed cohort, reference and processing parameters
import sys
sys.dont_write_bytecode = True
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

from fig1_E10_trim_pilot import (ADAPTER, BED, ENV_BIN, INDEX, MIN_LENGTH, PARAMETERS,
                        SOURCE, THREADS, command, digest, fastqc_metrics,
                        parse_star, require, save_json)

# Set each storage placeholder below to its own local file or directory.
ROOT = Path()
# Set the FastQC launcher and Java executable.
FASTQC_EXECUTABLE = ROOT
JAVA_EXECUTABLE = ROOT

# Set the file or directory for work.
WORK = ROOT
# Uses the verified nine-cell results from fig1_E10_trim_pilot.py.
PILOT = ROOT
# Set the manifest file; the released public manifests are in metadata/e10.
MANIFEST = ROOT
# Set normal_sra_manifest.tsv from metadata/e10.
NORMAL_MANIFEST = ROOT
# Set pilot_manifest.tsv from metadata/e10.
PILOT_MANIFEST = ROOT
# Uses full_quality_results.json from fig1_E10_recover_quality.py.
SOURCE_RESULTS = ROOT
# Set pilot_owner.json, written by fig1_E10_trim_pilot.py for the accepted pilot.
PILOT_OWNER_FILE = ROOT
# Set pilot_results.json, written by fig1_E10_trim_pilot.py.
PILOT_RESULTS_FILE = ROOT
# Set full_counts_results.json, written here and read by fig1_E10_Counts.py.
RESULTS_FILE = ROOT
# Set counts_unstranded.tsv.gz, written here and read by fig1_E10_Counts.py and fig1_E10_clustering.ipynb.
COUNTS_FILE = ROOT
BRANCH = "adapter_q20"
EXPECTED_CELLS, EXPECTED_REUSE, EXPECTED_NEW, EXPECTED_GENES = 73, 9, 64, 55487
MANIFEST_SHA256 = "239da5e99f2a1c9c82b906e06d39fc998096bcbde49b50431b4e2bc03058f59b"
NORMAL_SHA256 = "e556e9ccba2871a4bbe7b9cdff79e8b7354249a078b659436f147e9a3076ad85"
IDENTITY = ("cell_id_global", "original_cell_id", "plate_id", "GSM", "experiment_accession", "run_accession")
COUNT_COLUMNS = ("unstranded", "R1_sense", "R2_sense")
CATEGORIES = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]


def read_json(path):
    return json.loads(path.read_text())


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_path(path):
    """Reject symlinks in data/output paths below the existing work root."""
    relative = path.relative_to(SOURCE)
    cursor = SOURCE
    for part in relative.parts:
        cursor = cursor / part
        require(not cursor.is_symlink(), f"Refusing a symlink data/output path: {cursor}")
    return path


def check_file(path, expected, full_hash):
    safe_path(path)
    require(path.is_file() and path.stat().st_size == int(expected["bytes"]), f"Missing or resized file: {path}")
    require(re.fullmatch(r"[0-9a-f]{32}", expected["md5"]), f"Invalid recorded MD5: {path}")
    if full_hash:
        require(digest(path) == expected["md5"], f"File MD5 changed: {path}")


# Validate fixed identities, source metadata and source FASTQ files
def load_inputs():
    require(sha256(MANIFEST) == MANIFEST_SHA256, "The selected full-cohort manifest changed.")
    require(sha256(NORMAL_MANIFEST) == NORMAL_SHA256, "The original source manifest changed.")
    rows, normal_rows, pilot_rows = read_tsv(MANIFEST), read_tsv(NORMAL_MANIFEST), read_tsv(PILOT_MANIFEST)
    require(len(rows) == len(normal_rows) == EXPECTED_CELLS, "Expected exactly 73 selected source rows.")
    for field in IDENTITY:
        if field != "plate_id":
            require(len({r[field] for r in rows}) == EXPECTED_CELLS, f"Nonunique cohort identity: {field}")
    require(Counter(r["stagewise_annotation_code"] for r in rows) == {"O-2": 44, "O-3": 29}, "Annotation counts changed.")
    require({r["selection_reason"] for r in rows} == {"full_selected_cohort"}, "Unexpected cohort selection reason.")
    normal = {r["run_accession"]: r for r in normal_rows}
    require(set(normal) == {r["run_accession"] for r in rows}, "Full/source run sets differ.")
    by_run = {r["run_accession"]: r for r in rows}
    for row in rows:
        run = row["run_accession"]
        require(re.fullmatch(r"SRR[0-9]+", run), f"Invalid run name: {run}")
        require(all(row[k] == v for k, v in normal[run].items()), f"Source identity/metadata changed: {run}")
    require(len(pilot_rows) == len({r["run_accession"] for r in pilot_rows}) == EXPECTED_REUSE, "Expected nine unique pilot runs.")
    for row in pilot_rows:
        run = row["run_accession"]
        require(run in by_run and all(row[k] == v for k, v in normal[run].items()), f"Pilot/source mapping differs: {run}")
        require(row["stagewise_annotation_code"] == by_run[run]["stagewise_annotation_code"], f"Pilot annotation changed: {run}")
    sources = read_json(SOURCE_RESULTS)["runs"]
    require(set(sources) == set(by_run), "Full-quality FASTQ metadata must cover exactly the selected 73 runs.")
    for row in rows:
        run, source = row["run_accession"], sources[row["run_accession"]]
        require(all(source[k] == row[k] for k in IDENTITY), f"FASTQ cell mapping differs: {run}")
        require(source["fastqc_complete"] and source["paired_ids_verified"] and source["unfiltered"], f"Incomplete raw-read validation: {run}")
        require(source["expected_spots"] == int(row["total_spots"]), f"Original read count differs: {run}")
        require(source["normal_sra_url"] == row["normal_sra_url"] and source["normal_sra_md5"] == row["normal_sra_md5"]
                and source["normal_sra_bytes"] == int(row["normal_sra_bytes"]), f"Normalized SRA source differs: {run}")
        require([m["mate"] for m in source["mates"]] == [1, 2], f"Unexpected mate order: {run}")
        for mate in source["mates"]:
            require(mate["filename"] == f"{run}_{mate['mate']}.fastq.gz" and mate["records"] == source["expected_spots"]
                    and mate["ena_sequence_match"], f"Incomplete mate evidence: {run}")
    pilot_token = {"schema": "e10-trimming-pilot-v1", "manifest_md5": digest(PILOT_MANIFEST), "parameters": PARAMETERS}
    require(read_json(PILOT_OWNER_FILE) == pilot_token, "Pilot configuration no longer matches the validated parameters.")
    pilot_results = read_json(PILOT_RESULTS_FILE)
    require(all(pilot_results[k] == v for k, v in pilot_token.items()), "Pilot result configuration differs.")
    require(pilot_results["completed_comparisons"] == pilot_results["expected_comparisons"] == 18, "Pilot is incomplete.")
    require(len(pilot_results["cells"]) == 18 and len({(r["run"], r["branch"]) for r in pilot_results["cells"]}) == 18,
            "Pilot records are missing or duplicated.")
    token = {"schema": "e10-full-counts-v1", "manifest_sha256": MANIFEST_SHA256, "normal_manifest_sha256": NORMAL_SHA256,
             "source_results_md5": digest(SOURCE_RESULTS), "pilot_configuration": pilot_token,
             "parameters": PARAMETERS, "selected_branch": BRANCH, "STAR_count_column": 2,
             "expected_cells": EXPECTED_CELLS, "expected_genes": EXPECTED_GENES}
    return rows, {r["run_accession"]: r for r in pilot_rows}, sources, pilot_results, token


def check_sources(rows, sources, full_hash):
    for index, row in enumerate(rows, 1):
        run = row["run_accession"]
        for mate in sources[run]["mates"]:
            check_file(SOURCE / "fastq_full" / mate["filename"],
                       {"bytes": mate["gzip_bytes"], "md5": mate["gzip_md5"]}, full_hash)
        if full_hash and (index % 10 == 0 or index == EXPECTED_CELLS):
            print(f"Source FASTQ size/MD5 verified: {index}/{EXPECTED_CELLS} cells", flush=True)


# Parse STAR counts and validate completed cells before any reuse
def gene_counts(directory, retained):
    with (directory / "ReadsPerGene.out.tab").open() as handle:
        rows = [line.rstrip("\r\n").split("\t") for line in handle]
    require(len(rows) == EXPECTED_GENES + 4 and all(len(r) == 4 for r in rows), f"Unexpected STAR count dimensions: {directory}")
    require([r[0] for r in rows[:4]] == CATEGORIES, f"Unexpected STAR count categories: {directory}")
    genes = [r[0] for r in rows[4:]]
    require(len(set(genes)) == EXPECTED_GENES and all(genes), f"Duplicate or empty gene IDs: {directory}")
    require(all(re.fullmatch(r"[0-9]+", value) for row in rows for value in row[1:]), f"Noninteger/negative counts: {directory}")
    values = [[int(v) for v in r[1:]] for r in rows]
    require(all(sum(v[j] for v in values) == retained for j in range(3)), f"STAR category totals differ from retained pairs: {directory}")
    return genes, values


def cell_metrics(row, source, directory, with_fastqc=True):
    run = row["run_accession"]
    cutadapt = read_json(directory / "cutadapt.json")
    original, retained = source["expected_spots"], cutadapt["read_counts"]["output"]
    require(cutadapt["read_counts"]["input"] == original and 0 < retained <= original, f"Cutadapt pair-count mismatch: {run}")
    star = parse_star(directory / "Log.final.out")
    require(int(star["Number of input reads"]) == retained, f"STAR/Cutadapt count mismatch: {run}")
    genes, counts = gene_counts(directory, retained)
    strand_text = (directory / "strand.log").read_text()
    require("This is PairEnd Data" in strand_text, f"RSeQC did not identify paired-end data: {run}")
    fractions = [float(v) for v in re.findall(r"Fraction of reads.*?:\s*([0-9.]+)", strand_text)]
    sampled = re.search(r"Total (\d+) usable reads were sampled", strand_text)
    require(len(fractions) == 3 and all(0 <= v <= 1 for v in fractions) and abs(sum(fractions) - 1) <= 0.001,
            f"Invalid RSeQC strand fractions: {run}")
    require(sampled is not None and 0 < int(sampled.group(1)) <= PARAMETERS["strandedness_sample_alignments"], f"Invalid RSeQC sample size: {run}")
    metrics = {"original_pairs": original, "retained_pairs": retained, "retained_pair_fraction": retained / original,
               "retained_bases": cutadapt["basepair_counts"]["output"], "unique_mapped_pairs": int(star["Uniquely mapped reads number"]),
               "multimapped_pairs": int(star["Number of reads mapped to multiple loci"]), "star_metrics": star,
               "gene_counts_by_column": {name: {"assigned_pairs": sum(v[j] for v in counts[4:]),
                   "assigned_per_original_pair": sum(v[j] for v in counts[4:]) / original,
                   "assigned_per_retained_pair": sum(v[j] for v in counts[4:]) / retained,
                   "detected_genes": sum(v[j] > 0 for v in counts[4:]), "no_feature_pairs": counts[2][j],
                   "ambiguous_pairs": counts[3][j]} for j, name in enumerate(COUNT_COLUMNS)},
               "rseqc_strand_fractions": dict(zip(("undetermined", "R1_sense", "R2_sense"), fractions)),
               "rseqc_sampled_alignments": int(sampled.group(1))}
    if with_fastqc:
        metrics["fastqc"] = fastqc_metrics(directory / "fastqc", run, retained)
    return metrics, genes, [v[0] for v in counts[4:]]


def required_outputs(run):
    return {f"{run}_1.fastq.gz", f"{run}_2.fastq.gz", "Aligned.sortedByCoord.out.bam", "Aligned.out.bam",
            "Aligned.sortedByCoord.out.bam.bai", "ReadsPerGene.out.tab", "Log.final.out", "cutadapt.json", "strand.log"}


def validate_record(row, source, directory, token, full_hash):
    safe_path(directory)
    record = read_json(directory / "complete.json")
    require(record["configuration"] == token and record["branch"] == BRANCH, f"Completed configuration differs: {directory}")
    for key, column in (("run", "run_accession"), ("cell", "original_cell_id"), ("plate", "plate_id"),
                        ("annotation", "stagewise_annotation_code"), ("selection_reason", "selection_reason")):
        require(record[key] == row[column], f"Completed identity differs: {directory}/{key}")
    require(required_outputs(row["run_accession"]) <= set(record["output_files"]), f"Missing output checksum records: {directory}")
    for filename, expected in record["output_files"].items():
        relative = Path(filename)
        require(not relative.is_absolute() and ".." not in relative.parts, f"Unsafe output filename: {filename}")
        check_file(directory / relative, expected, full_hash)
    for mate in (1, 2):
        for suffix in ("zip", "html"):
            path = safe_path(directory / "fastqc" / f"{row['run_accession']}_{mate}_fastqc.{suffix}")
            require(path.is_file() and path.stat().st_size > 0, f"Missing trimmed FastQC output: {path}")
    metrics, genes, counts = cell_metrics(row, source, directory, with_fastqc=full_hash)
    require(all(record[key] == value for key, value in metrics.items()), f"Completed metrics differ from their files: {directory}")
    return record, genes, counts


# Check tools and output ownership without writing during preflight
def environment_and_tools():
    environment = os.environ.copy()
    environment.update(PATH=str(ENV_BIN) + os.pathsep + environment.get("PATH", ""), TMPDIR=str(WORK),
                       OMP_NUM_THREADS=str(THREADS), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    for tool, expected in (("cutadapt", PARAMETERS["Cutadapt_version"]), ("STAR", PARAMETERS["STAR_version"])):
        version = subprocess.check_output([str(ENV_BIN / tool), "--version"], env=environment, text=True).strip()
        require(version == expected, f"Unexpected {tool} version: {version}")
    for executable in (ENV_BIN / "samtools", ENV_BIN / "infer_experiment.py", ENV_BIN / "perl", JAVA_EXECUTABLE, FASTQC_EXECUTABLE):
        require(executable.is_file(), f"Missing processing tool: {executable}")
    require((INDEX / "complete.txt").is_file() and (INDEX / "Genome").stat().st_size > 0, "The original STAR index is incomplete.")
    require(BED.is_file() and BED.stat().st_size > 0, "The original pilot strand-reference BED is missing.")
    return environment


def output_safety(rows, pilot_rows, token):
    safe_path(WORK)
    owner = WORK / "full_counts_owner.json"
    if WORK.exists():
        require(WORK.is_dir() and owner.is_file() and read_json(owner) == token, "Unrecognized or changed full-count output directory.")
        permitted = {"full_counts_owner.json", "full_counts.lock", "full_counts_results.json", "counts_unstranded.tsv.gz", BRANCH}
        require(all(p.name in permitted for p in WORK.iterdir()), "Unexpected files or interrupted publication in full-count directory.")
    branch = safe_path(WORK / BRANCH)
    new_runs = {r["run_accession"] for r in rows} - set(pilot_rows)
    require(len(new_runs) == EXPECTED_NEW, "Expected exactly 64 nonpilot runs.")
    if branch.exists():
        require(branch.is_dir() and all(p.name in new_runs for p in branch.iterdir()), "Unexpected runs in the full-count branch directory.")
    completed = []
    for run in sorted(new_runs):
        directory = safe_path(branch / run)
        if directory.exists():
            require(directory.is_dir() and (directory / "complete.json").is_file(), f"Inspect incomplete cell directory before resuming; no overwrite: {directory}")
            completed.append(run)
    prior = RESULTS_FILE
    if prior.exists():
        require(read_json(prior)["configuration"] == token, "Existing full-count summary belongs to another configuration.")
    return completed


def existing_cells(rows, pilot_rows, sources, pilot_results, token, full_hash):
    entries, reference_genes = {}, None
    pilot_records = {(r["run"], r["branch"]): r for r in pilot_results["cells"]}
    for row in rows:
        run = row["run_accession"]
        reuse = run in pilot_rows
        directory = (PILOT if reuse else WORK) / BRANCH / run
        if not reuse and not directory.exists():
            continue
        expected_row = pilot_rows[run] if reuse else row
        config = token["pilot_configuration"] if reuse else token
        record, genes, _ = validate_record(expected_row, sources[run], directory, config, full_hash)
        if reuse:
            require(record == pilot_records[(run, BRANCH)], f"Pilot completion and pilot summary differ: {run}")
        if reference_genes is None:
            reference_genes = genes
        require(genes == reference_genes, f"Gene IDs/order differ: {run}")
        entries[run] = {"cell_id_global": row["cell_id_global"], "origin": "pilot_reused" if reuse else "full_counts",
                        "directory": str(directory), "completion_record": str(directory / "complete.json"), "record": record}
        if full_hash:
            print(f"Completed cell validated: {row['original_cell_id']} ({entries[run]['origin']})", flush=True)
    require(sum(v["origin"] == "pilot_reused" for v in entries.values()) == EXPECTED_REUSE, "Nine pilot cells must be reused.")
    return entries, reference_genes


# Run the validated Q20 workflow for each of the remaining 64 cells
def process_cell(row, source, token, environment):
    run = row["run_accession"]
    directory = safe_path(WORK / BRANCH / run)
    require(not directory.exists(), f"Refusing to overwrite existing cell directory: {directory}")
    directory.mkdir(parents=True)
    reads = [SOURCE / "fastq_full" / f"{run}_{mate}.fastq.gz" for mate in (1, 2)]
    trimmed = [directory / f"{run}_{mate}.fastq.gz" for mate in (1, 2)]
    command([ENV_BIN / "cutadapt", "-j", str(THREADS), "-a", ADAPTER, "-A", ADAPTER, "-e", "0.1", "-O", "3",
             "-m", str(MIN_LENGTH), "--pair-filter=any", "--json", directory / "cutadapt.json", "-q", "20",
             "-o", trimmed[0], "-p", trimmed[1]] + reads, directory / "cutadapt.log", environment)
    cutadapt = read_json(directory / "cutadapt.json")
    retained = cutadapt["read_counts"]["output"]
    require(cutadapt["read_counts"]["input"] == source["expected_spots"] and 0 < retained <= source["expected_spots"],
            f"Cutadapt read-pair count mismatch: {run}")
    qc_dir = directory / "fastqc"
    qc_dir.mkdir()
    command([ENV_BIN / "perl", FASTQC_EXECUTABLE, "--java", "/usr/bin/java", "--threads", "2",
             "--dir", directory, "--nogroup", "--quiet", "--noextract", "--outdir", qc_dir] + trimmed,
            directory / "fastqc.log", environment)
    fastqc_metrics(qc_dir, run, retained)
    command([ENV_BIN / "STAR", "--runThreadN", str(THREADS), "--genomeDir", INDEX, "--readFilesIn"] + trimmed +
            ["--readFilesCommand", "zcat", "--outFileNamePrefix", str(directory) + "/", "--outSAMtype", "BAM", "Unsorted",
             "SortedByCoordinate", "--limitBAMsortRAM", "8000000000", "--outSAMattributes", "NH", "HI", "AS", "nM", "MD",
             "--outSAMunmapped", "Within", "KeepPairs", "--quantMode", "GeneCounts"], directory / "star.log", environment)
    bam, unsorted = directory / "Aligned.sortedByCoord.out.bam", directory / "Aligned.out.bam"
    command([ENV_BIN / "samtools", "quickcheck", "-v", bam], directory / "bam_check.log", environment)
    command([ENV_BIN / "samtools", "index", "-@", "4", bam], directory / "bam_check.log", environment)
    command([ENV_BIN / "samtools", "quickcheck", "-v", unsorted], directory / "bam_check.log", environment)
    command([ENV_BIN / "infer_experiment.py", "-i", unsorted, "-r", BED, "-s", "200000", "-q", "30"], directory / "strand.log", environment)
    metrics, genes, _ = cell_metrics(row, source, directory)
    record = {"configuration": token, "run": run, "cell": row["original_cell_id"], "cell_id_global": row["cell_id_global"],
              "plate": row["plate_id"], "annotation": row["stagewise_annotation_code"], "selection_reason": row["selection_reason"],
              "branch": BRANCH, **metrics, "output_files": {}}
    paths = [directory / name for name in sorted(required_outputs(run))]
    paths += [qc_dir / f"{run}_{mate}_fastqc.{suffix}" for mate in (1, 2) for suffix in ("zip", "html")]
    for path in paths:
        record["output_files"][str(path.relative_to(directory))] = {"bytes": path.stat().st_size, "md5": digest(path)}
    save_json(directory / "complete.json", record)
    return {"cell_id_global": row["cell_id_global"], "origin": "full_counts", "directory": str(directory),
            "completion_record": str(directory / "complete.json"), "record": record}, genes


# Publish and verify a 73 by 55,487 integer unstranded count matrix
def aggregate(rows, entries, reference_genes):
    path = COUNTS_FILE
    require(len(entries) == EXPECTED_CELLS and len(reference_genes) == EXPECTED_GENES, "Cannot aggregate incomplete cells/genes.")
    if not path.exists():
        temporary = path.with_suffix(path.suffix + ".writing")
        require(not temporary.exists(), "Inspect interrupted count-matrix publication before resuming.")
        with gzip.open(temporary, "wt", encoding="utf-8", newline="") as output:
            output.write("cell_id_global\t" + "\t".join(reference_genes) + "\n")
            for row in rows:
                entry = entries[row["run_accession"]]
                genes, counts = gene_counts(Path(entry["directory"]), entry["record"]["retained_pairs"])
                require(genes == reference_genes, f"Gene IDs/order differ during aggregation: {row['run_accession']}")
                output.write(row["cell_id_global"] + "\t" + "\t".join(str(c[0]) for c in counts[4:]) + "\n")
        os.replace(temporary, path)
    total = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        require(handle.readline().rstrip("\n").split("\t") == ["cell_id_global"] + reference_genes, "Published matrix gene header differs.")
        for row in rows:
            fields = handle.readline().rstrip("\n").split("\t")
            require(len(fields) == EXPECTED_GENES + 1 and fields[0] == row["cell_id_global"], "Published matrix cell order/dimensions differ.")
            require(all(re.fullmatch(r"[0-9]+", v) for v in fields[1:]), "Published matrix contains noninteger counts.")
            values = [int(v) for v in fields[1:]]
            entry = entries[row["run_accession"]]
            genes, original = gene_counts(Path(entry["directory"]), entry["record"]["retained_pairs"])
            require(genes == reference_genes and values == [v[0] for v in original[4:]], "Published matrix differs from STAR column 2.")
            require(sum(values) == entry["record"]["gene_counts_by_column"]["unstranded"]["assigned_pairs"], "Published matrix row sum differs.")
            total += sum(values)
        require(handle.read() == "", "Published matrix contains extra cell rows.")
    return {"path": str(path), "format": "gzip TSV; cells are rows; GENCODE gene IDs are columns", "shape": [EXPECTED_CELLS, EXPECTED_GENES],
            "STAR_count_column": 2, "strandedness": "unstranded", "sum": total, "bytes": path.stat().st_size, "md5": digest(path)}


def result_summary(rows, entries, token, matrix=None):
    return {"configuration": token, "expected_cells": EXPECTED_CELLS, "completed_cells": len(entries),
            "reused_pilot_cells": sum(e["origin"] == "pilot_reused" for e in entries.values()),
            "completed_new_cells": sum(e["origin"] == "full_counts" for e in entries.values()),
            "cells": [entries[r["run_accession"]] for r in rows if r["run_accession"] in entries], "matrix": matrix}


# Read-only preflight or one locked production run
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true", help="Read-only metadata, sizes, pilot counts/configuration and tools; no large-file MD5 or processing.")
    args = parser.parse_args()
    rows, pilot_rows, sources, pilot_results, token = load_inputs()
    environment = environment_and_tools()
    completed = output_safety(rows, pilot_rows, token)
    check_sources(rows, sources, full_hash=False)
    if args.preflight:
        entries, genes = existing_cells(rows, pilot_rows, sources, pilot_results, token, full_hash=False)
        print(json.dumps({"preflight": "passed_read_only", "cells": len(rows), "annotations": {"O-2": 44, "O-3": 29},
                          "reuse": EXPECTED_REUSE, "completed_new": len(completed), "pending_new": EXPECTED_NEW - len(completed),
                          "source_FASTQs_size_checked": 146, "large_file_MD5": "required on compute node before processing",
                          "reference_genes": len(genes), "count_column": 2, "output": str(WORK)}, indent=2), flush=True)
        return
    require(os.environ.get("SLURM_JOB_ID"), "Run production through the Slurm job; large-file MD5 must run on a compute node.")
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:THREADS])
    WORK.mkdir(exist_ok=True)
    import fcntl
    with (WORK / "full_counts.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = WORK / "full_counts_owner.json"
        if not owner.exists():
            require(all(p.name == "full_counts.lock" for p in WORK.iterdir()), "Refusing unowned existing files.")
            save_json(owner, token)
        output_safety(rows, pilot_rows, token)
        check_sources(rows, sources, full_hash=True)
        entries, genes = existing_cells(rows, pilot_rows, sources, pilot_results, token, full_hash=True)
        save_json(RESULTS_FILE, result_summary(rows, entries, token))
        for row in rows:
            run = row["run_accession"]
            if run in entries:
                continue
            print(f"Processing {row['original_cell_id']} {run}; completed={len(entries)}/{EXPECTED_CELLS}", flush=True)
            entry, current_genes = process_cell(row, sources[run], token, environment)
            require(current_genes == genes, f"New cell gene IDs/order differ: {run}")
            entries[run] = entry
            save_json(RESULTS_FILE, result_summary(rows, entries, token))
            print(f"Counted {len(entries)}/{EXPECTED_CELLS}: {row['original_cell_id']}; unstranded assigned={entry['record']['gene_counts_by_column']['unstranded']['assigned_pairs']}", flush=True)
        require(sum(e["origin"] == "full_counts" for e in entries.values()) == EXPECTED_NEW, "Expected 64 completed new cells.")
        matrix = aggregate(rows, entries, genes)
        save_json(RESULTS_FILE, result_summary(rows, entries, token, matrix))
        print(json.dumps({"status": "complete", "cells": EXPECTED_CELLS, "pilot_reused": EXPECTED_REUSE,
                          "new_processed": EXPECTED_NEW, "matrix": matrix}, indent=2), flush=True)


if __name__ == "__main__":
    main()
