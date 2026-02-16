# ============================================================
# 03_tximport_pca
# Tximport + Normalization + PCA
# ============================================================

# =========================
# --- Package Management ---
# =========================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c("DESeq2", "tximport", "Biostrings")
cran_packages <- c("ggplot2", "dplyr", "plotly", "htmlwidgets", "stringr")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
  library(pkg, character.only = TRUE)
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =========================
# --- Sample Metadata ---
# =========================

metadata_path <- "data/metadata/metadata.tsv"

if (!file.exists(metadata_path)) {
  stop("Metadata file not found.")
}

samples <- read.delim(metadata_path, header = TRUE, stringsAsFactors = FALSE)

rownames(samples) <- samples$srr_id
samples$condition <- factor(samples$condition)

# =========================
# --- Salmon quant files ---
# =========================

base_dir <- "results/salmon_quants"

files <- file.path(base_dir, samples$srr_id, "quant.sf")
names(files) <- samples$srr_id

if (!all(file.exists(files))) {
  stop("Some quant.sf files are missing.")
}

# =========================
# --- Reference FASTA ---
# =========================

reference_fasta <- "data/reference/Homo_sapiens.GRCh38.cdna.all.fa"

if (!file.exists(reference_fasta)) {
  stop("Reference FASTA not found. Please download it as described in the README and place it in data/reference/")
}

# =========================
# --- Build tx2gene ---
# =========================

fa <- readDNAStringSet(reference_fasta)

tx_ids <- names(fa)

tx2gene <- data.frame(
  tx = tx_ids,
  gene = sub(".*gene:([^ ]+).*", "\\1", tx_ids),
  stringsAsFactors = FALSE
)

tx2gene_clean <- tx2gene %>%
  mutate(
    tx = str_split(tx, " ", simplify = TRUE)[,1],
    tx = gsub("\\..*", "", tx),
    gene = gsub("\\..*", "", gene)
  )

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
# --- DESeq2 object ---
# =========================

dds <- DESeqDataSetFromTximport(
  txi_gene,
  colData = samples,
  design = ~ condition
)

# =========================
# --- Filtering ---
# =========================

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# =========================
# --- VST ---
# =========================

vsd <- vst(dds, blind = TRUE)

# =========================
# --- PCA ---
# =========================

pca <- prcomp(t(assay(vsd)))
var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  condition = samples$condition
)

# =========================
# --- Output directory ---
# =========================

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# =========================
# --- PCA 2D ---
# =========================

pca_plot_2d <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA of Mast Cell RNA-seq Data (VST)",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  )

ggsave(
  "results/figures/PCA_2D.png",
  pca_plot_2d,
  width = 8,
  height = 6,
  dpi = 300
)

# =========================
# --- PCA 3D ---
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

saveWidget(
  pca_plot_3d,
  "results/figures/PCA_3D.html",
  selfcontained = TRUE
)

cat("tximport + PCA completed successfully.\n")
