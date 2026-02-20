#!/bin/bash
set -euo pipefail

# ======================================
# RNA-seq Quantification
# Project: PRJNA896634
# ======================================


# ======================================
# 1. Define directories and parameters
# ======================================

TRIM_DIR="results/tmp_trimmed"
REF_DIR="data/reference"
INDEX_DIR="$REF_DIR/salmon_index"
OUT_DIR="results/salmon_quants"

TRANSCRIPTOME="$REF_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz"
THREADS=8

mkdir -p "$OUT_DIR"
mkdir -p "$REF_DIR"


# ======================================
# 2. Download transcriptome if missing
# ======================================

if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "Downloading human transcriptome (Ensembl 109)..."
    curl -L ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o "$TRANSCRIPTOME"
fi


# ======================================
# 3. Create Salmon index (if not exists)
# ======================================

if [ ! -d "$INDEX_DIR" ] || [ -z "$(ls -A "$INDEX_DIR")" ]; then
    echo "Creating Salmon index..."
    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$INDEX_DIR" \
        -k 31 \
        --threads $THREADS
else
    echo "Salmon index already exists. Skipping."
fi


# ======================================
# 4. Quantification
# ======================================

echo "Starting quantification..."

# Loop through paired reads
for R1 in "$TRIM_DIR"/*_1_paired.fastq.gz
do
    R2=${R1/_1_paired.fastq.gz/_2_paired.fastq.gz}
    SAMPLE=$(basename "$R1" _1_paired.fastq.gz)

    if [ -f "$R2" ]; then
        echo "--------------------------------------------"
        echo "Processing sample: $SAMPLE"
        echo "--------------------------------------------"

       salmon quant \
            -i "$INDEX_DIR" \
            -l A \
            -1 "$R1" \
            -2 "$R2" \
            -p $THREADS \
            --validateMappings \
            --gcBias \
            --seqBias \
            -o "$OUT_DIR/$SAMPLE"

    else
        echo "ERROR: Missing R2 pair for $R1"
    fi
done

echo "Quantification finished."