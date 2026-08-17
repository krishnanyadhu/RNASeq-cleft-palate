library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load data
txi <- readRDS("txi_gencode_v49.rds")
metadata <- readRDS("sample_metadata.rds")

rownames(metadata) <- metadata$sample

# Verify sample order
stopifnot(identical(colnames(txi$counts), rownames(metadata)))

# Create DESeq2 object with paired design
dds <- DESeqDataSetFromTximport(
    txi,
    colData = metadata,
    design = ~ patient + condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 10
cat("Genes before filtering:", nrow(dds), "\n")
cat("Genes after filtering:", sum(keep), "\n")

dds <- dds[keep, ]

# Library sizes
library_size <- colSums(counts(dds))
write.csv(data.frame(sample = names(library_size), 
                     library_size = library_size), 
          "library_sizes.csv", row.names = FALSE)

# VST transformation
vsd <- vst(dds, blind = TRUE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = c("condition", "sex"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = sex)) +
    geom_point(size = 3) +
    labs(x = paste0("PC1 (", percent_var[1], "%)"),
         y = paste0("PC2 (", percent_var[2], "%)"),
         title = "PCA of RNA-seq samples") +
    theme_bw()
ggsave("PCA_condition_sex.png", p, width = 9, height = 7, dpi = 300)

# Sample correlation heatmap
cor_mat <- cor(assay(vsd), method = "pearson")

png("sample_correlation_heatmap.png", width = 2500, height = 2500, res = 300)
pheatmap(
    cor_mat,
    annotation_col = metadata[, c("condition", "sex")],
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    main = "Sample correlation"
)
dev.off()

# Save objects
saveRDS(dds, "dds_filtered.rds")
saveRDS(vsd, "vsd.rds")

cat("✓ QC completed successfully.\n")