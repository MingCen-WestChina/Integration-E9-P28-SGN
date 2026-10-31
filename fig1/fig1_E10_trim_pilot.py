#!/usr/bin/env python3
# Source: trim_pilot.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
"""Compare two paired-read trimming settings and annotated STAR gene counts."""

# Fixed pilot parameters
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import zipfile

# Set each storage placeholder below to its own local file or directory.
ROOT = Path()
# Set the FastQC launcher and Java executable.
FASTQC_EXECUTABLE = ROOT
JAVA_EXECUTABLE = ROOT
# Uses the released fig1_E10_gtf_bed12.py helper.
BED12_SCRIPT = ROOT
# Set pilot_results.json, written here and read by fig1_E10_full_counts.py.
PILOT_RESULTS_FILE = ROOT

# Set the input or working directory used by this step.
SOURCE = ROOT
# Set the file or directory for work.
WORK = ROOT
# Set the directory containing the required command-line tools.
ENV_BIN = ROOT
# Set the file or directory for reference.
REFERENCE = ROOT
# Set the STAR GRCm38 / GENCODE M25 index directory.
INDEX = ROOT
# Set the uncompressed GENCODE M25 primary-assembly GTF file.
GTF = ROOT
# Set the transcript BED12 output made by fig1_E10_gtf_bed12.py.
BED = ROOT
THREADS = 8
ADAPTER = "CTGTCTCTTATACACATCT"
MIN_LENGTH = 30
BRANCHES = {"adapter": [], "adapter_q20": ["-q", "20"]}
PARAMETERS = {"threads": THREADS, "adapter_R1_R2": ADAPTER, "adapter_error_rate": 0.1,
              "adapter_min_overlap": 3, "min_length_both_mates": MIN_LENGTH, "pair_filter": "any",
              "quality_cutoffs": {"adapter": None, "adapter_q20": [0, 20]},
              "STAR_version": "2.7.11b", "Cutadapt_version": "5.1", "RSeQC_version": "5.0.4", "genome": "GRCm38 primary",
              "annotation": "GENCODE M25 PRI", "STAR_two_pass": False, "STAR_sort_RAM_bytes": 8000000000,
              "strandedness_sample_alignments": 200000, "strandedness_mapq": 30,
              "strandedness_input": "STAR unsorted BAM, first qualifying alignments"}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path, data):
    temporary = path.with_suffix(path.suffix + ".writing")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    os.replace(temporary, path)


def command(argv, log, environment):
    with log.open("ab") as output:
        output.write(("COMMAND " + json.dumps([str(a) for a in argv]) + "\n").encode())
        output.flush()
        subprocess.run([str(a) for a in argv], stdout=output, stderr=subprocess.STDOUT,
                       env=environment, cwd=log.parent, check=True)


def parse_star(path):
    return {key.strip(): value.strip() for line in path.read_text().splitlines()
            if "|" in line for key, value in [line.split("|", 1)]}


def fastqc_metrics(directory, run, expected_pairs):
    metrics = []
    for mate in (1, 2):
        path = directory / f"{run}_{mate}_fastqc.zip"
        with zipfile.ZipFile(path) as archive:
            require(archive.testzip() is None, f"Corrupt FastQC ZIP: {path}")
            data = archive.read(path.stem + "/fastqc_data.txt").decode().split("\n", 1)[1]
        modules = {}
        for block in data.split(">>END_MODULE"):
            lines = block.strip().splitlines()
            if lines and lines[0].startswith(">>"):
                modules[lines[0].split("\t")[0][2:]] = lines[1:]
        basic = dict(line.split("\t", 1) for line in modules["Basic Statistics"] if not line.startswith("#"))
        require(int(basic["Total Sequences"]) == expected_pairs, f"FastQC/Cutadapt count mismatch: {run}/{mate}")
        quality = [line.split("\t") for line in modules["Per base sequence quality"][1:]]
        adapters = [line.split("\t") for line in modules["Adapter Content"][1:]]
        adapter_index = modules["Adapter Content"][0][1:].split("\t").index("Nextera Transposase Sequence")
        metrics.append({"mate": mate, "reads": expected_pairs, "sequence_length": basic["Sequence length"],
                        "quality_positions": [int(x[0]) for x in quality], "quality_mean": [float(x[1]) for x in quality],
                        "adapter_positions": [int(x[0]) for x in adapters],
                        "adapter_percent": [float(x[adapter_index]) for x in adapters]})
    return metrics


# Verify the selected cells, source files and fixed tools
def main():
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:THREADS])
    WORK.mkdir(exist_ok=True)
    require(not WORK.is_symlink(), "The pilot work directory must not be a symlink.")
    import fcntl
    lock = (WORK / "pilot.lock").open("a+")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    # Set metadata/e10/pilot_manifest.tsv, the released fixed nine-cell pilot manifest.
    manifest = ROOT
    rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
    # Set full_quality_results.json, generated by fig1_E10_recover_quality.py.
    original = json.loads((ROOT).read_text())["runs"]
    require(len(rows) == len({r["run_accession"] for r in rows}) == 9, "Expected nine unique pilot cells.")
    require(len({r["plate_id"] for r in rows}) == 6, "The pilot must cover all six plates.")
    require({r["stagewise_annotation_code"] for r in rows} == {"O-2", "O-3"}, "Both existing annotations must be represented.")
    token = {"schema": "e10-trimming-pilot-v1", "manifest_md5": digest(manifest), "parameters": PARAMETERS}
    owner = WORK / "pilot_owner.json"
    if owner.exists():
        require(json.loads(owner.read_text()) == token, "Existing pilot parameters differ; do not overwrite.")
    else:
        require(not any(p.name != "pilot.lock" for p in WORK.iterdir()), "Unrecognized nonempty pilot directory.")
        save_json(owner, token)
    environment = os.environ.copy()
    environment.update(PATH=str(ENV_BIN) + os.pathsep + environment["PATH"], TMPDIR=str(WORK),
                       OMP_NUM_THREADS=str(THREADS), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for tool, expected in (("cutadapt", "5.1"), ("STAR", "2.7.11b")):
        version = subprocess.check_output([str(ENV_BIN / tool), "--version"], env=environment, text=True).strip()
        require(version == expected, f"Unexpected {tool} version: {version}")
    require((INDEX / "complete.txt").is_file(), "The STAR index is incomplete.")
    if not BED.exists():
        command([sys.executable, BED12_SCRIPT, "--gtf", GTF, "--bed", BED], WORK / "reference.log", environment)
    require(BED.stat().st_size > 0, "Empty strand reference BED.")
    results = dict(token, cells=[], expected_comparisons=18, completed_comparisons=0)
    for row in rows:
        run = row["run_accession"]
        source = original[run]
        require(source["fastqc_complete"] and source["cell_id_global"] == row["cell_id_global"], f"Incomplete or mismatched source: {run}")
        reads = [SOURCE / "fastq_full" / f"{run}_{mate}.fastq.gz" for mate in (1, 2)]
        for path, mate in zip(reads, source["mates"]):
            require(path.stat().st_size == mate["gzip_bytes"] and digest(path) == mate["gzip_md5"], f"Source FASTQ changed: {path}")
        original_pairs = source["expected_spots"]

        #3. =========== Paired trimming, full FastQC, annotated STAR alignment and strand check
        for branch, quality_arguments in BRANCHES.items():
            directory = WORK / branch / run
            complete = directory / "complete.json"
            if complete.exists():
                record = json.loads(complete.read_text())
                require(record["configuration"] == token, f"Changed configuration: {directory}")
                for filename, expected in record["output_files"].items():
                    path = directory / filename
                    require(path.is_file() and path.stat().st_size == expected["bytes"] and digest(path) == expected["md5"],
                            f"Completed pilot output changed: {path}")
            else:
                require(not directory.exists(), f"Inspect the incomplete pilot directory before resuming: {directory}")
                directory.mkdir(parents=True)
                trimmed = [directory / f"{run}_{mate}.fastq.gz" for mate in (1, 2)]
                command([ENV_BIN / "cutadapt", "-j", str(THREADS), "-a", ADAPTER, "-A", ADAPTER,
                         "-e", "0.1", "-O", "3", "-m", str(MIN_LENGTH), "--pair-filter=any",
                         "--json", directory / "cutadapt.json"] + quality_arguments +
                        ["-o", trimmed[0], "-p", trimmed[1]] + reads, directory / "cutadapt.log", environment)
                cutadapt = json.loads((directory / "cutadapt.json").read_text())
                require(cutadapt["read_counts"]["input"] == original_pairs, f"Cutadapt input count mismatch: {run}")
                retained = cutadapt["read_counts"]["output"]
                require(0 < retained <= original_pairs, f"Invalid retained read-pair count: {run}")
                qc_dir = directory / "fastqc"
                qc_dir.mkdir()
                command([ENV_BIN / "perl", FASTQC_EXECUTABLE, "--java", "/usr/bin/java",
                         "--threads", "2", "--dir", directory, "--nogroup", "--quiet", "--noextract",
                         "--outdir", qc_dir] + trimmed, directory / "fastqc.log", environment)
                qc = fastqc_metrics(qc_dir, run, retained)
                command([ENV_BIN / "STAR", "--runThreadN", str(THREADS), "--genomeDir", INDEX,
                         "--readFilesIn"] + trimmed + ["--readFilesCommand", "zcat", "--outFileNamePrefix", str(directory) + "/",
                         "--outSAMtype", "BAM", "Unsorted", "SortedByCoordinate", "--limitBAMsortRAM", "8000000000",
                         "--outSAMattributes", "NH", "HI", "AS", "nM", "MD", "--outSAMunmapped", "Within", "KeepPairs",
                         "--quantMode", "GeneCounts"], directory / "star.log", environment)
                bam = directory / "Aligned.sortedByCoord.out.bam"
                command([ENV_BIN / "samtools", "quickcheck", "-v", bam], directory / "bam_check.log", environment)
                command([ENV_BIN / "samtools", "index", "-@", "4", bam], directory / "bam_check.log", environment)
                unsorted_bam = directory / "Aligned.out.bam"
                command([ENV_BIN / "samtools", "quickcheck", "-v", unsorted_bam], directory / "bam_check.log", environment)
                command([ENV_BIN / "infer_experiment.py", "-i", unsorted_bam, "-r", BED, "-s", "200000", "-q", "30"],
                        directory / "strand.log", environment)

                #4. =========== Compare absolute gene-assigned yields against the original pairs
                star = parse_star(directory / "Log.final.out")
                require(int(star["Number of input reads"]) == retained, f"STAR/Cutadapt count mismatch: {run}")
                with (directory / "ReadsPerGene.out.tab").open() as handle:
                    count_rows = [line.rstrip().split("\t") for line in handle]
                require([r[0] for r in count_rows[:4]] == ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"],
                        f"Unexpected STAR gene-count categories: {run}")
                counts = [[int(v) for v in r[1:]] for r in count_rows]
                require(all(len(c) == 3 and min(c) >= 0 for c in counts), f"Invalid gene counts: {run}")
                require(all(sum(c[j] for c in counts) == retained for j in range(3)), f"STAR gene-count totals differ from input: {run}")
                strand_text = (directory / "strand.log").read_text()
                require("This is PairEnd Data" in strand_text, f"RSeQC did not identify paired-end data: {run}")
                fractions = [float(v) for v in re.findall(r"Fraction of reads.*?:\s*([0-9.]+)", strand_text)]
                require(len(fractions) == 3, f"Incomplete RSeQC strand summary: {run}")
                sampled = re.search(r"Total (\d+) usable reads were sampled", strand_text)
                require(sampled is not None and int(sampled.group(1)) > 0, f"No usable RSeQC alignments: {run}")
                record = {"configuration": token, "run": run, "cell": row["original_cell_id"], "plate": row["plate_id"],
                          "annotation": row["stagewise_annotation_code"], "selection_reason": row["selection_reason"],
                          "branch": branch, "original_pairs": original_pairs, "retained_pairs": retained,
                          "retained_pair_fraction": retained / original_pairs,
                          "retained_bases": cutadapt["basepair_counts"]["output"],
                          "unique_mapped_pairs": int(star["Uniquely mapped reads number"]),
                          "multimapped_pairs": int(star["Number of reads mapped to multiple loci"]),
                          "star_metrics": star, "fastqc": qc,
                          "gene_counts_by_column": {name: {"assigned_pairs": sum(c[j] for c in counts[4:]),
                              "assigned_per_original_pair": sum(c[j] for c in counts[4:]) / original_pairs,
                              "assigned_per_retained_pair": sum(c[j] for c in counts[4:]) / retained,
                              "detected_genes": sum(c[j] > 0 for c in counts[4:]),
                              "no_feature_pairs": counts[2][j], "ambiguous_pairs": counts[3][j]}
                              for j, name in enumerate(("unstranded", "R1_sense", "R2_sense"))},
                          "rseqc_strand_fractions": {"undetermined": fractions[0], "R1_sense": fractions[1], "R2_sense": fractions[2]},
                          "rseqc_sampled_alignments": int(sampled.group(1)),
                          "output_files": {}}
                for path in trimmed + [bam, unsorted_bam, bam.with_suffix(".bam.bai"), directory / "ReadsPerGene.out.tab",
                                       directory / "Log.final.out", directory / "cutadapt.json", directory / "strand.log"]:
                    record["output_files"][path.name] = {"bytes": path.stat().st_size, "md5": digest(path)}
                save_json(complete, record)
            results["cells"].append(record)
            results["completed_comparisons"] += 1
            save_json(PILOT_RESULTS_FILE, results)
            print(f"Pilot complete {results['completed_comparisons']}/18: {row['original_cell_id']} {branch}; "
                  f"retained={record['retained_pair_fraction']:.4f}; "
                  f"unstranded_gene_yield={record['gene_counts_by_column']['unstranded']['assigned_per_original_pair']:.4f}", flush=True)
    require(results["completed_comparisons"] == 18, "Pilot comparisons are incomplete.")
    print("Pilot complete: nine cells, two trimming settings. No all-cell counting or integration performed.", flush=True)


if __name__ == "__main__":
    main()
