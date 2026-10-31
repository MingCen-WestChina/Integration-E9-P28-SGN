#!/usr/bin/env python3
# Source: gtf_bed12.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
"""Convert GENCODE GTF to transcript BED12 for RSeQC strandedness inference.

The thick interval spans the whole transcript; it does not represent true CDS.
Do not interpret CDS/UTR categories from read_distribution.py using this BED.
All chromosomes and biotypes are retained. Only Python's standard library is used.
Format references: https://rseqc.sourceforge.net/#input-format
https://genome.ucsc.edu/FAQ/FAQformat.html#format1
"""

# Imports and transcript records
import argparse
from dataclasses import dataclass, field
import gzip
import math
import os
from pathlib import Path
import re
import tempfile

# Set each storage placeholder below to its own local file or directory.
ROOT = Path()


@dataclass
class Transcript:
    chrom: str
    strand: str
    gene_id: str
    exons: list = field(default_factory=list)
    bounds: tuple | None = None


ATTRIBUTE = re.compile(r'\s*([A-Za-z_]\w*)\s+(?:"([^"]*)"|([^;\s]+))\s*;')


# GTF parsing and validation
def parse_ids(text, line_number):
    identifiers = {}
    position = 0
    while position < len(text):
        if not text[position:].strip():
            break
        match = ATTRIBUTE.match(text, position)
        if match is None:
            raise ValueError(f"Line {line_number}: malformed GTF attributes")
        key = match.group(1)
        value = match.group(2) if match.group(2) is not None else match.group(3)
        if key in {"gene_id", "transcript_id"}:
            if key in identifiers or not value or re.search(r"\s", value):
                raise ValueError(f"Line {line_number}: invalid or repeated {key}")
            identifiers[key] = value
        position = match.end()
    if "gene_id" not in identifiers:
        raise ValueError(f"Line {line_number}: missing gene_id")
    return identifiers


def read_gtf(gtf):
    transcripts = {}
    exon_count = 0
    opener = gzip.open if gtf.suffix.lower() == ".gz" else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Line {line_number}: expected nine GTF columns")
            chrom, source, feature, first, last, score, strand, phase, attributes = fields
            if not chrom or chrom == "." or re.search(r"\s", chrom):
                raise ValueError(f"Line {line_number}: invalid chromosome")
            if not source or not feature or strand not in {"+", "-"}:
                raise ValueError(f"Line {line_number}: invalid source, feature, or strand")
            if phase not in {".", "0", "1", "2"}:
                raise ValueError(f"Line {line_number}: invalid GTF phase")
            try:
                if not first.isdecimal() or not last.isdecimal():
                    raise ValueError
                start, end = int(first) - 1, int(last)
                if start < 0 or end <= start:
                    raise ValueError
                if score != "." and not math.isfinite(float(score)):
                    raise ValueError
            except ValueError as error:
                raise ValueError(f"Line {line_number}: invalid coordinates or score") from error
            identifiers = parse_ids(attributes, line_number)
            transcript_id = identifiers.get("transcript_id")
            if transcript_id is None:
                if feature != "gene":
                    raise ValueError(f"Line {line_number}: missing transcript_id")
                continue
            identity = (chrom, strand, identifiers["gene_id"])
            model = transcripts.get(transcript_id)
            if model is None:
                model = transcripts[transcript_id] = Transcript(*identity)
            elif (model.chrom, model.strand, model.gene_id) != identity:
                raise ValueError(
                    f"Line {line_number}: inconsistent chromosome, strand, or gene for {transcript_id}"
                )
            if feature == "transcript":
                if model.bounds is not None:
                    raise ValueError(f"Line {line_number}: duplicate transcript record {transcript_id}")
                model.bounds = (start, end)
            elif feature == "exon":
                model.exons.append((start, end))
                exon_count += 1
    if not transcripts or not exon_count:
        raise ValueError("No transcript exons found")
    return transcripts, exon_count


# Exon merging and BED12 construction
def bed_records(transcripts):
    records = []
    block_count = 0
    for transcript_id, model in transcripts.items():
        if not model.exons:
            raise ValueError(f"Transcript {transcript_id}: no exons")
        merged = []
        for start, end in sorted(model.exons):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        start, end = merged[0][0], merged[-1][1]
        if model.bounds is not None and model.bounds != (start, end):
            raise ValueError(f"Transcript {transcript_id}: exon span differs from transcript bounds")
        sizes = [right - left for left, right in merged]
        offsets = [left - start for left, _ in merged]
        if any(size <= 0 or offset < 0 or offset + size > end - start
               for size, offset in zip(sizes, offsets)):
            raise ValueError(f"Transcript {transcript_id}: invalid BED block bounds")
        fields = [model.chrom, start, end, transcript_id, 0, model.strand,
                  start, end, 0, len(merged),
                  ",".join(map(str, sizes)), ",".join(map(str, offsets))]
        records.append((model.chrom, start, end, transcript_id, "\t".join(map(str, fields))))
        block_count += len(merged)
    records.sort(key=lambda row: row[:4])
    return records, block_count


# Atomic output and command line
def write_bed(bed, records):
    if os.path.lexists(bed):
        raise FileExistsError(f"Refusing existing output: {bed}")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                         prefix=bed.name + ".", suffix=".tmp",
                                         dir=bed.parent, delete=False) as handle:
            temporary = Path(handle.name)
            for record in records:
                handle.write(record[4] + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        # Same-directory hard link publishes the complete file atomically, without replacement.
        os.link(temporary, bed)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, type=Path, help="GENCODE GTF or GTF.gz")
    parser.add_argument("--bed", required=True, type=Path, help="New BED12 output path")
    args = parser.parse_args()
    if not args.gtf.is_file():
        parser.error(f"GTF file not found: {args.gtf}")
    if os.path.lexists(args.bed):
        parser.error(f"Refusing existing output: {args.bed}")
    if not args.bed.parent.is_dir():
        parser.error(f"Output directory not found: {args.bed.parent}")
    try:
        transcripts, exon_count = read_gtf(args.gtf)
        records, block_count = bed_records(transcripts)
        write_bed(args.bed, records)
    except (OSError, ValueError, EOFError) as error:
        parser.exit(1, f"Error: {error}\n")
    gene_count = len({model.gene_id for model in transcripts.values()})
    print(f"genes={gene_count} transcripts={len(records)} exon_records={exon_count} "
          f"bed_blocks={block_count} merged_exon_records={exon_count - block_count}")


if __name__ == "__main__":
    main()
