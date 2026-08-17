library(data.table)

# Files
gtf_file <- "reference/gencode.v50.annotation.gtf"
de_file <- "DESeq2_C_vs_NC_shrunken_all.csv"

# Read DESeq2 results
de <- fread(de_file)
cat("DESeq2 genes:", nrow(de), "\n")

# Read GTF
gtf <- fread(gtf_file, sep = "\t", header = FALSE, 
             comment.char = "#", quote = "")
setnames(gtf, c("chr", "source", "feature", "start", "end",
                "score", "strand", "frame", "attribute"))

# Keep transcript records
tx <- gtf[feature == "transcript"]

# Extract gene information
tx[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attribute)]
tx[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attribute)]
tx[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', attribute)]

# Create annotation
annotation <- unique(tx[, .(gene_id, gene_name, gene_type)])
annotation <- annotation[gene_id != "" & gene_name != "" & gene_type != ""]

cat("Unique annotated genes:", uniqueN(annotation$gene_id), "\n")

# Merge annotation with DE results
de_annotated <- merge(de, annotation, by = "gene_id", all.x = TRUE)
setorder(de_annotated, padj)

# Define DEGs
sig <- !is.na(de_annotated$padj) & 
       de_annotated$padj < 0.05 & 
       abs(de_annotated$log2FoldChange) >= 1
up <- sig & de_annotated$log2FoldChange >= 1
down <- sig & de_annotated$log2FoldChange <= -1

# Save annotated results
fwrite(de_annotated, "DESeq2_C_vs_NC_ANNOTATED_ALL.csv")
fwrite(de_annotated[sig], "DEG_ANNOTATED.csv")
fwrite(de_annotated[up], "DEG_UP_ANNOTATED.csv")
fwrite(de_annotated[down], "DEG_DOWN_ANNOTATED.csv")

# Gene lists for enrichment
fwrite(de_annotated[sig, .(gene_id, gene_name)], "DEG_gene_list.csv")
fwrite(de_annotated[up, .(gene_id, gene_name)], "DEG_UP_gene_list.csv")
fwrite(de_annotated[down, .(gene_id, gene_name)], "DEG_DOWN_gene_list.csv")

cat("\nAnnotation completed.\n")
cat("Total genes:", nrow(de_annotated), "\n")
cat("DEGs:", sum(sig), "\n")
cat("Upregulated:", sum(up), "\n")
cat("Downregulated:", sum(down), "\n")
