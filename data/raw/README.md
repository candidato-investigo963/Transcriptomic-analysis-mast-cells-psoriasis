# Raw RNA-Seq Data

This folder stores the raw RNA-Seq FASTQ files for the Mast Cell Transcriptomics project (psoriasis vs healthy controls).

## Study Design

- 13 total samples
  - 6 psoriasis
  - 7 healthy controls
- Sequencing type: Illumina paired-end
- Expected read length: 150 bp

## Content

- Paired-end FASTQ files:
  - `*_1.fastq.gz` → Forward reads
  - `*_2.fastq.gz` → Reverse reads
- `metadata.tsv` → Sample metadata (condition, accession ID, replicate, etc.)

## Data Source

Raw data available from:

- ENA Project: PRJNA896634
- GEO Series: GSE217060

## Download Instructions

Example using SRA Toolkit:

```bash
prefetch SRRXXXXXXX
fasterq-dump SRRXXXXXXX --split-files --gzip
```

Place the resulting files in: data/raw/


Ensure files follow the naming convention:

SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz


