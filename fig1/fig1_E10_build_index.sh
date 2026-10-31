#!/usr/bin/env bash
set -euo pipefail

# Fixed reference and observed read length
# Set the working directory.
ROOT=.
WORK="$ROOT"
# Set the reference directory.
REF="$ROOT"
# Set GRCm38.primary_assembly.genome.fa.gz, the downloaded GRCm38 primary-assembly genome.
GENOME_GZIP="$ROOT"
# Set gencode.vM25.primary_assembly.annotation.gtf.gz, the downloaded GENCODE M25 primary-assembly annotation.
GTF_GZIP="$ROOT"
# Set the expanded GRCm38.primary_assembly.genome.fa output file.
GENOME_FASTA="$ROOT"
# Set the expanded gencode.vM25.primary_assembly.annotation.gtf output file.
GTF_FILE="$ROOT"
# Set the STAR index output directory.
INDEX="$ROOT"
# Set the directory containing STAR and the reference utilities.
ENV_BIN="$ROOT"
export PATH="$ENV_BIN:$PATH"
export TMPDIR="$WORK"
cd "$REF"
printf '3bc591be24b77f710b6ba5d41022fc5a  %s\nc5125258a0a2c5250ddb4c192abbf4e8  %s\n' "$GENOME_GZIP" "$GTF_GZIP" | md5sum -c -

# Expand matching primary-assembly reference files
if ! test -f "$GENOME_FASTA"; then gzip -dc "$GENOME_GZIP" > "${GENOME_FASTA}.writing"; mv "${GENOME_FASTA}.writing" "$GENOME_FASTA"; fi
if ! test -f "$GTF_FILE"; then gzip -dc "$GTF_GZIP" > "${GTF_FILE}.writing"; mv "${GTF_FILE}.writing" "$GTF_FILE"; fi
if test -e "$INDEX"; then printf 'Index path already exists; inspect it before rerunning.\n' >&2; exit 1; fi
mkdir "$INDEX"

# Build a splice-aware STAR index; no read alignment
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir "$INDEX" --genomeFastaFiles "$GENOME_FASTA" --sjdbGTFfile "$GTF_FILE" --sjdbOverhang 74 --genomeSAindexNbases 14 --limitGenomeGenerateRAM 60000000000 --outFileNamePrefix "$INDEX/"
test -s "$INDEX/Genome" && test -s "$INDEX/SA" && test -s "$INDEX/SAindex" && test -s "$INDEX/geneInfo.tab"
printf 'STAR 2.7.11b; GRCm38 primary; GENCODE M25 PRI; sjdbOverhang=74; threads=8\n' > "$INDEX/complete.txt"
