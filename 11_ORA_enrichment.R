library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# Read DEG files
all <- fread("DESeq2_C_vs_NC_ANNOTATED_ALL.csv")
up <- fread("DEG_UP_ANNOTATED.csv")
down <- fread("DEG_DOWN_ANNOTATED.csv")

cat("Background genes:", nrow(all), "\n")
cat("UP genes:", nrow(up), "\n")
cat("DOWN genes:", nrow(down), "\n")

# Remove Ensembl version numbers
strip_version <- function(x) sub("\\..*$", "", x)

background <- unique(strip_version(all$gene_id))
up_genes <- unique(strip_version(up$gene_id))
down_genes <- unique(strip_version(down$gene_id))

# Map to Entrez IDs
background_entrez <- unique(na.omit(mapIds(
    org.Hs.eg.db, keys = background, keytype = "ENSEMBL",
    column = "ENTREZID", multiVals = "first"
)))
up_entrez <- unique(na.omit(mapIds(
    org.Hs.eg.db, keys = up_genes, keytype = "ENSEMBL",
    column = "ENTREZID", multiVals = "first"
)))
down_entrez <- unique(na.omit(mapIds(
    org.Hs.eg.db, keys = down_genes, keytype = "ENSEMBL",
    column = "ENTREZID", multiVals = "first"
)))

# Create output directory
dir.create("enrichment_ORA", showWarnings = FALSE)

# GO Enrichment Function
run_GO <- function(gene_list, background_list, label) {
    ego_BP <- enrichGO(gene = gene_list, universe = background_list,
                       OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                       readable = TRUE)
    ego_MF <- enrichGO(gene = gene_list, universe = background_list,
                       OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                       ont = "MF", pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                       readable = TRUE)
    ego_CC <- enrichGO(gene = gene_list, universe = background_list,
                       OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                       ont = "CC", pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                       readable = TRUE)
    
    fwrite(as.data.table(ego_BP), paste0("enrichment_ORA/GO_BP_", label, ".csv"))
    fwrite(as.data.table(ego_MF), paste0("enrichment_ORA/GO_MF_", label, ".csv"))
    fwrite(as.data.table(ego_CC), paste0("enrichment_ORA/GO_CC_", label, ".csv"))
    
    return(list(BP = ego_BP, MF = ego_MF, CC = ego_CC))
}

# KEGG Enrichment Function
run_KEGG <- function(gene_list, background_list, label) {
    ekegg <- enrichKEGG(gene = gene_list, universe = background_list,
                        organism = "hsa", keyType = "ncbi-geneid",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
    fwrite(as.data.table(ekegg), paste0("enrichment_ORA/KEGG_", label, ".csv"))
    return(ekegg)
}

# Reactome Enrichment Function
run_Reactome <- function(gene_list, background_list, label) {
    er <- enrichPathway(gene = gene_list, universe = background_list,
                        organism = "human", pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                        readable = TRUE)
    fwrite(as.data.table(er), paste0("enrichment_ORA/Reactome_", label, ".csv"))
    return(er)
}

# Run enrichment for UP and DOWN genes
GO_UP <- run_GO(up_entrez, background_entrez, "UP")
KEGG_UP <- run_KEGG(up_entrez, background_entrez, "UP")
Reactome_UP <- run_Reactome(up_entrez, background_entrez, "UP")

GO_DOWN <- run_GO(down_entrez, background_entrez, "DOWN")
KEGG_DOWN <- run_KEGG(down_entrez, background_entrez, "DOWN")
Reactome_DOWN <- run_Reactome(down_entrez, background_entrez, "DOWN")

cat("\n✓ Enrichment analysis completed.\n")