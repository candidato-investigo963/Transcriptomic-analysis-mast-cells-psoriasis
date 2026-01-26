###
##### Transcriptomic Analysis of Mast Cell Gene Expression in Psoriasis
###

# =========================
# --- Libraries ---
# =========================
library(tximport)
library(DESeq2)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library(org.Hs.eg.db)
library(EnhancedVolcano)

# =========================
# --- Sample Metadata ---
# =========================
samples <- data.frame(
  row.names = c(
    "SRR22134160", "SRR22134161", "SRR22134162",
    "SRR22134163", "SRR22134164", "SRR22134165",
    "SRR22134166", "SRR22134167", "SRR22134168",
    "SRR22134169", "SRR22134170", "SRR22134171",
    "SRR22134172"
  ),
  condition = c(rep("Psoriasis", 6), rep("Healthy", 7))
)

# =========================
# --- Salmon quant files ---
# =========================
# Base directory where Salmon output folders are stored
base_dir <- "Path_hidden_for_privacy"

files <- file.path(base_dir, rownames(samples), "quant.sf")
names(files) <- rownames(samples)

# Check if reference FASTA files are available
list.files(pattern = "\\.fa$")

# =========================
# --- Build tx2gene from FASTA ---
# =========================
# Load the transcriptome used by Salmon
fa <- readDNAStringSet("Homo_sapiens.GRCh38.cdna.all.fa")
tx_ids <- names(fa)

# Create transcript-to-gene mapping table
tx2gene <- data.frame(
  tx = tx_ids,
  gene = sub(".*gene:([^ ]+).*", "\\1", tx_ids),
  stringsAsFactors = FALSE
)

# Remove transcripts without a gene ID
tx2gene <- tx2gene[!is.na(tx2gene$gene), ]

# Clean transcript and gene IDs
tx2gene_clean <- tx2gene %>%
  mutate(
    # Keep only the transcript ID (before any space)
    tx = str_split(tx, " ", simplify = TRUE)[,1],
    # Remove version numbers (e.g. .1, .2)
    tx = gsub("\\..*", "", tx),
    gene = gsub("\\..*", "", gene)
  )

# Quick check to confirm matching works
any(tx2gene_clean$tx %in% c("ENST00000415118", "ENST00000448914"))

# =========================
# --- Import Salmon data ---
# =========================
txi_gene <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene_clean,
  ignoreTxVersion = TRUE
)

# =========================
# --- Create DESeq2 object ---
# =========================
dds <- DESeqDataSetFromTximport(
  txi_gene,
  colData = samples,
  design = ~ condition
)

# =========================
# --- Low expression filter ---
# =========================
# Keep genes with at least 10 counts in 3 or more samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# =========================
# --- VST Normalization ---
# =========================
vsd <- vst(dds, blind = TRUE)

# =========================
# --- PCA Analysis ---
# =========================
pca <- prcomp(t(assay(vsd)))

# Percentage of variance explained
var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
round(var_exp[1:3], 2)

# Build PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  condition = samples$condition
)

# =========================
# --- PCA 2D Plot ---
# =========================
pca_plot_2d <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA of mast cell RNA-seq data (genes, VST)",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  )

# Save PCA 2D figure
ggsave(
  filename = "results/figures/PCA_2D.png",
  plot = pca_plot_2d,
  width = 8,
  height = 6,
  dpi = 300
)

# =========================
# --- PCA 3D Plot ---
# =========================
pca_plot_3d <- plot_ly(
  pca_df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~condition,
  type = "scatter3d",
  mode = "markers"
)

# Save interactive PCA
saveWidget(
  pca_plot_3d,
  "results/figures/PCA_3D.html",
  selfcontained = TRUE
)

# =========================
# --- Differential Expression Analysis ---
# =========================
dds <- DESeq(dds)

# Psoriasis vs Healthy contrast
res <- results(dds, contrast = c("condition", "Psoriasis", "Healthy"))

# Shrink log2 fold changes
res_shrink <- lfcShrink(
  dds,
  coef = "condition_Psoriasis_vs_Healthy",
  type = "apeglm"
)

# =========================
# --- Annotation ---
# =========================
res_df <- as.data.frame(res_shrink) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db,
                    keys = gene_id,
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first"),
    genetype = mapIds(org.Hs.eg.db,
                      keys = gene_id,
                      column = "GENETYPE",
                      keytype = "ENSEMBL",
                      multiVals = "first")
  ) %>%
  dplyr::filter(!is.na(padj), !is.na(symbol))

# Keep only protein-coding genes
res_df_prot <- res_df %>%
  dplyr::filter(genetype == "protein-coding")

# =========================
# --- High-confidence gene selection ---
# =========================
ids_to_label <- res_df_prot %>%
  dplyr::filter(padj < 1e-6 & abs(log2FoldChange) > 2) %>%
  pull(gene_id)

# Create labels only for high-confidence genes
res_df_prot <- res_df_prot %>%
  mutate(label_volcano = ifelse(gene_id %in% ids_to_label, symbol, ""))

# =========================
# --- Volcano Plot ---
# =========================
# Keep only significant genes
res_clean <- res_df_prot %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Manual color assignment
colors_manual <- ifelse(
  res_clean$padj < 1e-6 & abs(res_clean$log2FoldChange) > 2,
  "#d62728",  # Red = high confidence
  "#7fbf7b"   # Green = significant
)

names(colors_manual)[colors_manual == "#d62728"] <- "p < 1e-6, |LFC| > 2"
names(colors_manual)[colors_manual == "#7fbf7b"] <- "p < 0.05, |LFC| > 1"

# Final volcano plot
EnhancedVolcano(
  res_clean,
  lab = res_clean$label_volcano,
  x = "log2FoldChange",
  y = "padj",
  colCustom = colors_manual,
  pCutoff = 1e-6,
  FCcutoff = 1.0,
  selectLab = res_df_prot$symbol[res_df_prot$gene_id %in% ids_to_label],
  pointSize = 2.5,
  labSize = 4.0,
  labCol = "black",
  labFace = "bold",
  colAlpha = 0.6,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = 50,
  title = "Psoriasis vs Healthy: Protein-Coding Genes",
  subtitle = "Significant (Green) vs High-Confidence Biomarkers (Red)",
  legendPosition = "top",
  legendLabSize = 12,
  legendIconSize = 5.0
)

# =========================
# --- Gene list for functional analysis ---
# =========================
genes_for_functional <- res_df_prot %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)


# =========================
# --- Biomarker Table ---
# =========================
biomarker_table <- res_df_prot %>%
  dplyr::filter(gene_id %in% ids_to_label) %>%
  dplyr::arrange(padj) %>%
  dplyr::select(symbol, log2FoldChange, padj) %>%
  dplyr::mutate(Direction = ifelse(
    log2FoldChange > 0,
    "Upregulated",
    "Downregulated"
  ))

# Print table
print(biomarker_table)

# Save table
write.csv(
  biomarker_table,
  "results/Tables/Biomarker_Table_Psoriasis.csv",
  row.names = FALSE
)

