# ============================================================
# 04_differential_expression_volcano.R
# Differential Expression + Annotation + Volcano Plot
# ============================================================


# ============================================================
# --- Load Required Libraries (auto-install if missing) ---
# ============================================================

required_packages <- c(
  "DESeq2",
  "apeglm",
  "EnhancedVolcano",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "tidyverse"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
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


# ============================================================
# --- Annotation ---
# ============================================================

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


# ============================================================
# --- High-Confidence Gene Selection ---
# ============================================================

ids_to_label <- res_df_prot %>%
  dplyr::filter(padj < 1e-6 & abs(log2FoldChange) > 2) %>%
  pull(gene_id)

# Create labels only for high-confidence genes
res_df_prot <- res_df_prot %>%
  mutate(label_volcano = ifelse(gene_id %in% ids_to_label, symbol, ""))


# ============================================================
# --- Volcano Plot ---
# ============================================================

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


# ============================================================
# --- Gene List for Functional Analysis ---
# ============================================================

genes_for_functional <- res_df_prot %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)


# ============================================================
# --- Biomarker Table ---
# ============================================================

biomarker_table <- res_df_prot %>%
  dplyr::filter(gene_id %in% ids_to_label) %>%
  dplyr::arrange(padj) %>%
  dplyr::select(symbol, log2FoldChange, padj) %>%
  dplyr::mutate(Direction = ifelse(
    log2FoldChange > 0,
    "Upregulated",
    "Downregulated"
  ))

print(biomarker_table)

write.csv(
  biomarker_table,
  "results/Tables/Biomarker_Table_Psoriasis.csv",
  row.names = FALSE
)
