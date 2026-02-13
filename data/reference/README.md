# Reference Data

This directory contains the reference transcriptome (Ensembl release 110, GRCh38 assembly) and the corresponding Salmon index used for RNA-seq quantification.

## Content

- `Homo_sapiens.GRCh38.cdna.all.fa.gz` → Reference transcriptome (Ensembl release 110)
- `Homo_sapiens.GRCh38.110.gtf.gz` → Gene annotation file (Ensembl release 110)
- `salmon_index/` → Salmon index generated from the transcriptome

## Notes

- The Salmon index is created automatically by the `02_quantification_salmon.sh` script if it does not exist.
- Large reference files are excluded from version control via `.gitignore`.

## Download Instructions

To download the reference files:

```bash
cd data/reference

wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
