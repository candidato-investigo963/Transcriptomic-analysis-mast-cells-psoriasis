# ============================================================
# 05_functional_enrichment.R
# Functional Enrichment
# ============================================================


# ============================================================
# --- Load Required Libraries ---
# ============================================================

packages <- c(
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "tidyverse"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(paste("Package not installed:", p))
  }
  library(p, character.only = TRUE)
}

# Create output folders (evita errores)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# ============================================================
# --- Gene Lists ---
# ============================================================

# Upregulated genes
genes_up <- genes_for_functional %>%
  dplyr::filter(log2FoldChange > 0) %>%
  pull(symbol)

# Downregulated genes
genes_down <- genes_for_functional %>%
  dplyr::filter(log2FoldChange < 0) %>%
  pull(symbol)

convert_to_entrez <- function(gene_symbols) {
  entrez <- mapIds(
    org.Hs.eg.db,
    keys = gene_symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  na.omit(entrez)
}

entrez_up <- convert_to_entrez(genes_up)
entrez_down <- convert_to_entrez(genes_down)


# ============================================================
# --- GO Enrichment Function ---
# ============================================================

run_GO <- function(gene_list) {
  enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
}


# ============================================================
# --- GO Enrichment (BP / CC / MF) ---
# ============================================================

ego_up_BP   <- enrichGO(entrez_up, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP")
ego_up_CC   <- enrichGO(entrez_up, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="CC")
ego_up_MF   <- enrichGO(entrez_up, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="MF")

ego_down_BP <- enrichGO(entrez_down, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP")
ego_down_CC <- enrichGO(entrez_down, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="CC")
ego_down_MF <- enrichGO(entrez_down, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="MF")


# Plot solo si hay resultados (evita crash)
if (!is.null(ego_up_BP) && nrow(as.data.frame(ego_up_BP)) > 0) {
  dotplot(ego_up_BP, showCategory = 15) +
    ggtitle("GO BP - Upregulated Genes")
}


# Save tables
write.csv(as.data.frame(ego_up_BP),   "results/tables/GO_UP_BP.csv", row.names = FALSE)
write.csv(as.data.frame(ego_up_CC),   "results/tables/GO_UP_CC.csv", row.names = FALSE)
write.csv(as.data.frame(ego_up_MF),   "results/tables/GO_UP_MF.csv", row.names = FALSE)

write.csv(as.data.frame(ego_down_BP), "results/tables/GO_DOWN_BP.csv", row.names = FALSE)
write.csv(as.data.frame(ego_down_CC), "results/tables/GO_DOWN_CC.csv", row.names = FALSE)
write.csv(as.data.frame(ego_down_MF), "results/tables/GO_DOWN_MF.csv", row.names = FALSE)


# ============================================================
# --- KEGG Enrichment ---
# ============================================================

ekegg_up <- enrichKEGG(
  gene         = entrez_up,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg_down <- enrichKEGG(
  gene         = entrez_down,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg_up   <- setReadable(ekegg_up, OrgDb = org.Hs.eg.db)
ekegg_down <- setReadable(ekegg_down, OrgDb = org.Hs.eg.db)

if (!is.null(ekegg_up) && nrow(as.data.frame(ekegg_up)) > 0) {
  dotplot(ekegg_up, showCategory = 15) +
    ggtitle("KEGG - Upregulated Genes")
}

write.csv(as.data.frame(ekegg_up),   "results/tables/KEGG_UP.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_down), "results/tables/KEGG_DOWN.csv", row.names = FALSE)
