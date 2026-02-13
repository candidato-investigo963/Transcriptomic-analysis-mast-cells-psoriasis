# Raw RNA-Seq Data

This folder is intended to store the raw RNA-Seq FASTQ files for the Mast Cell Transcriptomics project (psoriasis vs healthy controls).

## Content

- FASTQ files for each sample in paired-end format:
  - `*_1.fastq.gz` → Forward reads
  - `*_2.fastq.gz` → Reverse reads
- 13 samples in total:
  - 6 psoriasis
  - 7 healthy controls

## Instructions

1. Download the raw FASTQ files from the public repository (ENA / GEO):
   - ENA Project: [PRJNA896634](https://www.ebi.ac.uk/ena/browser/view/PRJNA896634)
   - GEO Series: [GSE217060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217060)

2. Place the downloaded files in this folder (`data/raw/`).

3. Ensure that files follow the paired-end naming convention:

SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz

