# Reference Data

This folder contains the reference transcriptome and Salmon index used for quantification.

## Content

- `Homo_sapiens.GRCh38.cdna.all.fa.gz` → Reference transcriptome (Ensembl release 110)
- `salmon_index/` → Salmon index generated from the transcriptome

## Note

- The Salmon index is created automatically by the `02_quantification_salmon.sh` script if it does not exist.
