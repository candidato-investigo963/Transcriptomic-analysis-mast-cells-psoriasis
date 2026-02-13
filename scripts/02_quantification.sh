#!/bin/bash
set -euo pipefail

# ======================================
# RNA-seq Quantification
# Project: PRJNA896634
# ======================================

# 1. Directories

RAW_DIR="../data/raw"
REF_DIR="../reference"
INDEX_DIR="$REF_DIR/salmon_index"
OUT_DIR="../results/salmon_quants"

TRANSCRIPTOME="$REF_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz"

THREADS=8

mkdir -p "$OUT_DIR"
mkdir -p "$REF_DIR"


# 1. Create Salmon index

if [ ! -d "$INDEX_DIR" ]; then
    echo "Salmon index not found. Creating index..."

    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$INDEX_DIR" \
        -k 31

else
    echo "Salmon index already exists. Skipping indexing."
fi


# 2. Show FASTQ files

echo "FASTQ files detected:"
ls "$RAW_DIR"/*_1.fastq.gz
echo "-----------------------------------"


# 3. Quantification loop

for R1 in "$RAW_DIR"/*_1.fastq.gz
do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    echo "Processing sample: $SAMPLE"

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p $THREADS \
        --validateMappings \
        -o "$OUT_DIR/$SAMPLE"

done

echo "Quantification finished."
