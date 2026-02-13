#!/bin/bash

set -euo pipefail

# ======================================
#           01 - PREPROCESSING
# ======================================

# 1. RUN ACCESSIONS

RUNS=(
SRR22134160 SRR22134161 SRR22134162 SRR22134163
SRR22134164 SRR22134165 SRR22134166 SRR22134167
SRR22134168 SRR22134169 SRR22134170 SRR22134171
SRR22134172
)

# 2. DIRECTORIES

RAW_DIR="../data/raw"
QC_RAW="../results/qc/raw"
QC_TRIM="../results/qc/trimmed"
TRIM_DIR="../results/tmp_trimmed"

ADAPTERS="/path/to/TruSeq3-PE.fa"

mkdir -p "$RAW_DIR" "$QC_RAW" "$QC_TRIM" "$TRIM_DIR"


# 3. DOWNLOAD FASTQ FILES

echo "Downloading FASTQ files..."

cd "$RAW_DIR"

for SRR in "${RUNS[@]}"; do
    echo "Processing $SRR"
    
    # Descarga en formato FASTQ
    fasterq-dump $SRR --split-files --gzip
    
done

cd -


# 4. QC RAW READS

echo "Running FastQC (raw)..."
fastqc "$RAW_DIR"/*.fastq.gz -o "$QC_RAW"
multiqc "$QC_RAW" -o "$QC_RAW"


# 5. TRIMMING

echo "Trimming reads..."

for R1 in "$RAW_DIR"/*_1.fastq.gz
do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    trimmomatic PE \
        "$R1" "$R2" \
        "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" /dev/null \
        "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" /dev/null \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done


# 6. QC TRIMMED READS

echo "Running FastQC (trimmed)..."
fastqc "$TRIM_DIR"/*.fastq.gz -o "$QC_TRIM"
multiqc "$QC_TRIM" -o "$QC_TRIM"

echo "Preprocessing completed."
