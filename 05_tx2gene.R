library(data.table)

gtf <- "reference/gencode.v50.annotation.gtf"

# Read GTF
g <- fread(
    gtf,
    sep = "\t",
    header = FALSE,
    comment.char = "#"
)

setnames(g, c(
    "chr", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
))

# Keep transcript rows only
tx <- g[feature == "transcript"]

# Extract transcript and gene IDs
tx[, transcript_id := sub('.*transcript_id "([^"]+)".*', '\\1', attribute)]
tx[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attribute)]

# Create tx2gene mapping
tx2gene <- unique(tx[, .(transcript_id, gene_id)])
tx2gene <- tx2gene[!is.na(transcript_id) & !is.na(gene_id) & 
                   transcript_id != "" & gene_id != ""]

# Save
fwrite(tx2gene, "tx2gene_gencode_v50.csv")

cat("Unique transcripts:", uniqueN(tx2gene$transcript_id), "\n")
cat("Unique genes:", uniqueN(tx2gene$gene_id), "\n")
