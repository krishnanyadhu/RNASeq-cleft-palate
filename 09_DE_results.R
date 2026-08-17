library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)

# Load DESeq2 object
dds <- readRDS("dds_DESeq2_C_vs_NC.rds")
dds$condition <- relevel(dds$condition, ref = "NC")

# LFC shrinkage
res <- lfcShrink(dds, coef = "condition_C_vs_NC", type = "apeglm")
res <- res[order(res$padj), ]

# Convert to data.frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", 
                     "lfcSE", "pvalue", "padj")]

# Save all shrunken results
write.csv(res_df, "DESeq2_C_vs_NC_shrunken_all.csv", row.names = FALSE, quote = FALSE)

# Define significant genes
sig <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1
up <- sig & res_df$log2FoldChange >= 1
down <- sig & res_df$log2FoldChange <= -1

# Save DEG lists
write.csv(res_df[sig, ], "DEG_all_FDR0.05_LFC1.csv", row.names = FALSE, quote = FALSE)
write.csv(res_df[up, ], "DEG_up_FDR0.05_LFC1.csv", row.names = FALSE, quote = FALSE)
write.csv(res_df[down, ], "DEG_down_FDR0.05_LFCminus1.csv", row.names = FALSE, quote = FALSE)

# Summary
cat("\n==============================\n")
cat("DE RESULTS SUMMARY\n")
cat("==============================\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("FDR < 0.05:", sum(!is.na(res_df$padj) & res_df$padj < 0.05), "\n")
cat("FDR < 0.05 & |LFC| >= 1:", sum(sig), "\n")
cat("Upregulated in C:", sum(up), "\n")
cat("Downregulated in C:", sum(down), "\n")

# Volcano plot
volcano_df <- res_df
volcano_df$significance <- "Not significant"
volcano_df$significance[!is.na(volcano_df$padj) & 
                       volcano_df$padj < 0.05 & 
                       volcano_df$log2FoldChange >= 1] <- "Upregulated"
volcano_df$significance[!is.na(volcano_df$padj) & 
                       volcano_df$padj < 0.05 & 
                       volcano_df$log2FoldChange <= -1] <- "Downregulated"
volcano_df$plot_padj <- pmax(volcano_df$padj, 1e-300)

p_volcano <- ggplot(volcano_df, 
    aes(x = log2FoldChange, y = -log10(plot_padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    labs(title = "Cleft vs Non-cleft",
         x = "Shrunken log2 fold change",
         y = "-log10 adjusted p-value")
ggsave("DESeq2_volcano.png", p_volcano, width = 9, height = 7, dpi = 300)

# Heatmaps
vsd <- vst(dds, blind = FALSE)

# Top 50 DEGs
top50 <- head(res_df[!is.na(res_df$padj), ], 50)
top50_genes <- top50$gene_id

heatmap_mat <- assay(vsd)[top50_genes, , drop = FALSE]
heatmap_mat <- t(scale(t(heatmap_mat)))

annotation <- data.frame(
    Condition = colData(vsd)$condition,
    Sex = colData(vsd)$sex,
    Patient = colData(vsd)$patient
)
rownames(annotation) <- colnames(vsd)

png("DESeq2_top50_heatmap.png", width = 3000, height = 3500, res = 300)
pheatmap(
    heatmap_mat,
    annotation_col = annotation,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 6,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    main = "Top 50 DE genes"
)
dev.off()

# Save
saveRDS(res, "DESeq2_C_vs_NC_shrunken.rds")

cat("✓ All analyses completed.\n")