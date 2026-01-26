# Transcriptomic-Analysis-Mast-Cells-Psoriasis

Comparative transcriptomic analysis of mast cells isolated from psoriasis patients and healthy donors using RNA-Seq data (GSE217060).

## 1. Project Overview

This project presents a complete RNA-seq analysis pipeline to identify differentially expressed genes in mast cells isolated from psoriasis patients and healthy donors. The study focuses on uncovering potential biomarker candidates associated with inflammatory and immune-related pathways involved in psoriasis.

The workflow follows a standard end-to-end bioinformatics approach, starting from raw sequencing data and progressing through quality control, preprocessing, transcript quantification, statistical analysis, and biological interpretation.

## 2. Dataset Description

The dataset used in this study corresponds to GEO accession GSE217060 and consists of paired-end RNA-seq data obtained from isolated human mast cells, with a total of 13 samples:

- Healthy donors: 7
- Psoriasis patients: 6

Sequencing type: Paired-end RNA-seq
Total FASTQ files: 26 (2 per biological sample)


## 3. Project Objectives

- To perform a complete RNA-seq processing pipeline starting from raw sequencing data
- To assess and improve read quality using standardized quality control methods such as adapter trimming
- To quantify transcript-level expression efficiently using a pseudo-alignment strategy
- To explore global transcriptional patterns between healthy and diseased conditions
- To identify significantly differentially expressed genes
- To highlight potential biomarker candidates associated with psoriasis-related pathways

## 4. Workflow Overview

The analysis follows a modular and reproducible workflow:

- Data acquisition from public repositories
- Quality control of raw sequencing reads
- Adapter trimming and preprocessing
- Transcript quantification using pseudo-alignment (Salmon)
- Exploratory data analysis using PCA
- Differential gene expression analysis
- Functional and pathway-level interpretation (Pathway analysis and GSEA)

## 5. Analysis Steps


### 5.1 Data Acquisition

Raw sequencing data were downloaded from the European Nucleotide Archive (ENA) using the GEO accession number as a reference. Each sample contains two FASTQ files corresponding to forward and reverse reads. All paired-end FASTQ files were organized by sample and condition (healthy vs. psoriasis).


### 5.2 Quality Control

Initial quality assessment was performed on all 26 FASTQ files using FastQC, followed by a global summary using MultiQC.

Key metrics evaluated included:

- Per-base sequence quality
- Adapter content
- GC content distribution
- Sequence duplication levels
- Overrepresented sequences

Most samples showed high and consistent quality across read positions, with minor adapter presence and some typical RNA-seq duplication patterns.


### 5.3 Preprocessing

Adapter trimming was performed to remove residual sequencing adapters and low-quality bases. The trimming process removed only a small fraction of reads, confirming that the original data quality was already high.


### 5.4 Alignment and Gene Expression Quantification

Instead of full genome alignment, a pseudo-alignment strategy was used for efficiency computational. Salmon was employed to quantify transcript-level expression directly against a reference transcriptome, enabling fast and accurate estimation of transcript abundances without generating alignment files.

- Reference: Reference transcriptome downloaded from Ensembl.

- Output: One **quant.sf** file per sample, containing expression matrices where rows correspond to genes/transcripts and columns correspond to samples.


### 5.5 Exploratory Data Analysis (PCA)

Principal Component Analysis (PCA) was performed to visualize global transcriptional differences between conditions, identify clustering patterns, and detect potential outliers.

The PCA revealed a clear separation between healthy and psoriasis samples, with one potential outlier corresponding to a sample with lower read depth and altered GC content.


### 5.6 Differential Gene Expression Analysis

Differential expression analysis was conducted using DESeq2.

Results were visualized using:

- MA plots (mean expression vs. log fold change)
- Volcano plots (log2 fold change vs. adjusted p-value)

Genes were considered significantly differentially expressed and marked in **green** based on the following criteria:

- Adjusted p-value < **0.05**
- |log2 fold change| ≥ **1**

Genes that exceeded a more stringent second threshold were marked in **red** as potential biomarker candidates:

- Adjusted p-value < **1e-6**
- |log2 fold change| ≥ **2**


### 5.7 Functional Enrichment Analysis
## 6. Results
## 7. Conclusion
## 8. References

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217060

https://www.ebi.ac.uk/ena/browser/view/PRJNA896634

West PW, Tontini C, Atmoko H, Kiss O, Garner T, Bahri R, Warren RB, Griffiths CEM, Stevens A, Bulfone-Paus S. Human Mast Cells Upregulate Cathepsin B, a Novel Marker of Itch in Psoriasis. Cells. 2023 Aug 30;12(17):2177. doi: 10.3390/cells12172177. PMID: 37681909; PMCID: PMC10486964.

