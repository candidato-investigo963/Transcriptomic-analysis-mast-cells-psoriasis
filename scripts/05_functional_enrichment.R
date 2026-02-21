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
  "tidyverse",
  "pathview"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(paste("Package not installed:", p))
  }
  library(p, character.only = TRUE)
}

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(pathview)



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

ego_up_BP   <- enrichGO(entrez_up, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP", readable=TRUE)
ego_up_CC   <- enrichGO(entrez_up, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="CC", readable=TRUE)
ego_up_MF   <- enrichGO(entrez_up, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="MF", readable=TRUE)

ego_down_BP <- enrichGO(entrez_down, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP", readable=TRUE)
ego_down_CC <- enrichGO(entrez_down, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="CC", readable=TRUE)
ego_down_MF <- enrichGO(entrez_down, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="MF", readable=TRUE)

# Plot and save UP - BP
if (!is.null(ego_up_BP) && nrow(as.data.frame(ego_up_BP)) > 0) {
  # Simplify removes redundant GO terms
  ego_up_BP_sim <- clusterProfiler::simplify(ego_up_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  png("results/figures/GO_UP_BP_dotplot.png", width = 800, height = 600)
  print(dotplot(ego_up_BP_sim, showCategory = 15) + ggtitle("GO BP - Upregulated Genes"))
  dev.off()
}

# Plot and save DOWN - BP
if (!is.null(ego_down_BP) && nrow(as.data.frame(ego_down_BP)) > 0) {
  # Simplify removes redundant GO terms
  ego_down_BP_sim <- clusterProfiler::simplify(ego_down_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  png("results/figures/GO_DOWN_BP_dotplot.png", width = 800, height = 600)
  print(dotplot(ego_down_BP_sim, showCategory = 15) + ggtitle("GO BP - Downregulated Genes"))
  dev.off()
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

ekegg_up <- enrichKEGG(gene = entrez_up, organism = "hsa", pvalueCutoff = 0.05)
ekegg_down <- enrichKEGG(gene = entrez_down, organism = "hsa", pvalueCutoff = 0.05)

# Process UP results
if (!is.null(ekegg_up) && nrow(as.data.frame(ekegg_up)) > 0) {
  ekegg_up <- setReadable(ekegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  png("results/figures/KEGG_UP_dotplot.png", width = 800, height = 600)
  print(dotplot(ekegg_up, showCategory = 15) + ggtitle("KEGG - Upregulated Genes"))
  dev.off()
  write.csv(as.data.frame(ekegg_up), "results/tables/KEGG_UP_results.csv", row.names = FALSE)
}

# Process DOWN results
if (!is.null(ekegg_down) && nrow(as.data.frame(ekegg_down)) > 0) {
  ekegg_down <- setReadable(ekegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  png("results/figures/KEGG_DOWN_dotplot.png", width = 800, height = 600)
  print(dotplot(ekegg_down, showCategory = 15) + ggtitle("KEGG - Downregulated Genes"))
  dev.off()
  write.csv(as.data.frame(ekegg_down), "results/tables/KEGG_DOWN_results.csv", row.names = FALSE)
}

# KEGG PATHVIEW
genelist_pathview <- genes_for_functional$log2FoldChange
names(genelist_pathview) <- mapIds(org.Hs.eg.db, 
                                   keys = genes_for_functional$symbol, 
                                   column = "ENTREZID", 
                                   keytype = "SYMBOL")


genelist_pathview <- na.omit(genelist_pathview)


my_target_pathways <- c(
  "hsa04514", # Cell adhesion molecules (UP)
  "hsa04080", # Neuroactive ligand-receptor (UP)
  "hsa03050", # Proteasome (DOWN)
  "hsa04141", # Protein processing in ER (DOWN)
  "hsa05012"  # Parkinson disease 
)


dir.create("results/figures/kegg_maps", showWarnings = FALSE)
setwd("results/figures/kegg_maps")

for (pid in my_target_pathways) {
  pathview(gene.data  = genelist_pathview,
           pathway.id = pid,
           species    = "hsa",
           limit      = list(gene=2, cpd=1), 
           low        = "green",             
           mid        = "gray",              
           high       = "red",               
           kegg.native = TRUE)               
}


setwd("../../../")
