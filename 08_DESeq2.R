library(DESeq2)

# Load filtered DESeq2 object
dds <- readRDS("dds_filtered.rds")

# Ensure reference level is NC
dds$condition <- relevel(dds$condition, ref = "NC")

# Paired design
design(dds) <- ~ patient + condition

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "C", "NC"), alpha = 0.05)
res <- res[order(res$padj), ]

# Summary
cat("DESeq2 results summary:\n")
print(summary(res))

# Significant genes
sig <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1
cat("FDR < 0.05 & |log2FC| >= 1:", sum(sig), "\n")
cat("Upregulated in C:", sum(sig & res$log2FoldChange > 1), "\n")
cat("Downregulated in C:", sum(sig & res$log2FoldChange < -1), "\n")

# Save results
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", 
                     "lfcSE", "stat", "pvalue", "padj")]

write.csv(res_df, "DESeq2_C_vs_NC_all_genes.csv", row.names = FALSE)

# Save significant results
sig_df <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & 
                 abs(res_df$log2FoldChange) >= 1, ]
write.csv(sig_df, "DESeq2_C_vs_NC_significant.csv", row.names = FALSE)

# MA plot
png("DESeq2_MAplot_C_vs_NC.png", width = 2400, height = 2000, res = 300)
plotMA(res, alpha = 0.05, ylim = c(-5, 5), 
       main = "Cleft vs Non-cleft")
dev.off()

# Save DESeq2 object
saveRDS(dds, "dds_DESeq2_C_vs_NC.rds")
saveRDS(res, "DESeq2_results_C_vs_NC.rds")

cat("✓ DESeq2 completed successfully.\n")