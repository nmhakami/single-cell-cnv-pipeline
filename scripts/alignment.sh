#!/bin/bash

# Align paired-end sequencing reads to the reference genome using BWA-MEM.S
set -euo pipefail

# Check required arguments.
if [ "$#" -lt 4 ]; then
    echo "Usage: bash alignment.sh <reference> <read1.fastq.gz> <read2.fastq.gz> <output_prefix> [threads]"
    exit 1
fi
# Usage:
# bash alignment.sh <reference> <read1.fastq.gz> <read2.fastq.gz> <output_prefix> [threads]

REFERENCE=$1
READ1=$2
READ2=$3
OUTPUT=$4
THREADS=${5:-8}

bwa mem -t "$THREADS" "$REFERENCE" "$READ1" "$READ2" \
    > "${OUTPUT}.sam" \
    2> "${OUTPUT}_bwa.log"

# Convert SAM to compressed BAM.
samtools view -@ "$THREADS" -b \
    -o "${OUTPUT}.bam" \
    "${OUTPUT}.sam"

# Remove PCR duplicates.

# Group paired reads together.
samtools collate -@ "$THREADS" \
    -o "${OUTPUT}_collate.bam" \
    "${OUTPUT}.bam"

# Add/fix mate-pair information required for duplicate detection.
samtools fixmate -@ "$THREADS" -m \
    "${OUTPUT}_collate.bam" \
    "${OUTPUT}_fixmate.bam"

# Sort alignments by genomic coordinate.
samtools sort -@ "$THREADS" \
    -o "${OUTPUT}_sorted.bam" \
    "${OUTPUT}_fixmate.bam"

# Remove PCR duplicates.
samtools markdup -@ "$THREADS" -r \
    "${OUTPUT}_sorted.bam" \
    "${OUTPUT}_rmdup.bam"

# Keep only uniquely mapped reads with mapping quality >= 30.
samtools view -@ "$THREADS" \
    -q 30 \
    -F 0x800 \
    -b \
    -o "${OUTPUT}_unique.bam" \
    "${OUTPUT}_rmdup.bam"

# Keep only the first read of each paired-end fragment.
samtools view -@ "$THREADS" \
    -f 0x40 \
    -h \
    -b \
    -o "${OUTPUT}_fwd.bam" \
    "${OUTPUT}_unique.bam"

# Summarize alignment statistics for quality control.
samtools flagstat -@ "$THREADS" \
    "${OUTPUT}_fwd.bam" \
    > "${OUTPUT}_flagstat.txt"

# Calculate insert-size metrics.
java -jar picard.jar CollectInsertSizeMetrics \
    I="${OUTPUT}_unique.bam" \
    O="${OUTPUT}_insert_size.txt" \
    H="${OUTPUT}_insert_size.pdf"ok