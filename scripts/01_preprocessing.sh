#!/bin/bash

set -euo pipefail


exec > >(tee -i preprocessing.log)
exec 2>&1

# ======================================
#           01 - PREPROCESSING
# ======================================

echo "Starting preprocessing pipeline..."

# ======================================
# PROJECT ROOT
# ======================================

PROJECT_DIR=$(pwd)

# ======================================
# 1. RUN ACCESSIONS
# ======================================

RUNS=(
SRR22134160 SRR22134161 SRR22134162 SRR22134163
SRR22134164 SRR22134165 SRR22134166 SRR22134167
SRR22134168 SRR22134169 SRR22134170 SRR22134171
SRR22134172
)

# ======================================
# 2. DIRECTORIES
# ======================================

RAW_DIR="$PROJECT_DIR/data/raw"
QC_RAW="$PROJECT_DIR/results/qc/raw"
QC_TRIM="$PROJECT_DIR/results/qc/trimmed"
TRIM_DIR="$PROJECT_DIR/results/tmp_trimmed"
ADAPTERS="$PROJECT_DIR/data/reference/TruSeq3-PE.fa"

mkdir -p "$RAW_DIR" "$QC_RAW" "$QC_TRIM" "$TRIM_DIR"

# ======================================
# 3. DOWNLOAD FASTQ FILES
# ======================================

echo "Downloading FASTQ files..."

cd "$RAW_DIR"

for SRR in "${RUNS[@]}"; do
    if [ ! -f "${SRR}_1.fastq.gz" ]; then
        echo "Downloading $SRR..."

        fasterq-dump "$SRR" --split-files --threads 8 --progress

        echo "Compressing FASTQ..."
        gzip ${SRR}_*.fastq

    else
        echo "$SRR already exists, skipping download."
    fi
done

cd "$PROJECT_DIR"

# ======================================
# 4. QC RAW READS
# ======================================

echo "Running FastQC (raw)..."

fastqc "$RAW_DIR"/*.fastq.gz -o "$QC_RAW" -t 8
multiqc "$QC_RAW" -o "$QC_RAW"

# ======================================
# 5. TRIMMING
# ======================================

echo "Trimming reads..."

for R1 in "$RAW_DIR"/*_1.fastq.gz; do

    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    # Verificar que existe el par
    if [ ! -f "$R2" ]; then
        echo "ERROR: Missing pair for $R1"
        exit 1
    fi

    echo "Processing sample: $SAMPLE"

    trimmomatic PE -threads 8 -phred33 \
        "$R1" "$R2" \
        "$TRIM_DIR/${SAMPLE}_1_paired.fastq.gz" "$TRIM_DIR/${SAMPLE}_1_unpaired.fastq.gz" \
        "$TRIM_DIR/${SAMPLE}_2_paired.fastq.gz" "$TRIM_DIR/${SAMPLE}_2_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done

# ======================================
# 6. QC TRIMMED READS
# ======================================

echo "Running FastQC (trimmed)..."

fastqc "$TRIM_DIR"/*_paired.fastq.gz -o "$QC_TRIM" -t 8
multiqc "$QC_TRIM" -o "$QC_TRIM"

echo "Generating global MultiQC report..."
multiqc "$PROJECT_DIR/results" -o "$PROJECT_DIR/results"

echo "Preprocessing completed successfully."
