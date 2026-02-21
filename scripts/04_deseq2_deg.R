# ============================================================
# 04_deseq2_deg.R
# Differential Expression
# ============================================================


# ============================================================
# --- Load Required Libraries ---
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("DESeq2", "apeglm", "EnhancedVolcano", "org.Hs.eg.db", "AnnotationDbi")
cran_pkgs <- c("tidyverse")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
  library(pkg, character.only = TRUE)
}
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}


# ============================================================
# --- Differential Expression Analysis ---
# ============================================================

dds <- DESeq(dds)

# Psoriasis vs Healthy contrast
res <- results(dds, contrast = c("condition", "Psoriasis", "Healthy"))

# Shrink log2 fold changes
res_shrink <- lfcShrink(
  dds,
  coef = "condition_Psoriasis_vs_Healthy",
  type = "apeglm"
)

# MA Plot
png("results/figures/MA_Plot_Psoriasis.png", width = 800, height = 600)
plotMA(res_shrink, ylim=c(-5,5), main="MA Plot: Psoriasis vs Healthy (shrunken)")
dev.off()

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
#res_clean <- res_df_prot %>%
  #dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Manual color assignment
colors_manual <- ifelse(
  res_df_prot$padj < 1e-6 & abs(res_df_prot$log2FoldChange) > 2,
  "#d62728",  # Red
  ifelse(
    res_df_prot$padj < 0.05 & abs(res_df_prot$log2FoldChange) > 1,
    "#7fbf7b",  # Green
    "grey70"    # Grey
  )
)

names(colors_manual)[colors_manual == "#d62728"] <- "High confidence"
names(colors_manual)[colors_manual == "#7fbf7b"] <- "Significant"
names(colors_manual)[colors_manual == "grey70"] <- "Not significant"

# Final volcano plot
EnhancedVolcano(
  res_df_prot,
  lab = res_df_prot$label_volcano,
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
  subtitle = "Red: p < 1e-6 & |LFC| > 2 | Green: p < 0.05 & |LFC| > 1",
  legendPosition = "top",
  legendLabSize = 12,
  legendIconSize = 5.0
)

ggsave("results/figures/Volcano_Psoriasis_vs_Healthy.png", volcano_plot, width = 10, height = 8)

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

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Save table
write.csv(biomarker_table, "results/tables/Biomarker_Table_Psoriasis.csv", row.names = FALSE)
