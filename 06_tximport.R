library(tximport)
library(data.table)

# Paths
salmon_dir <- "salmon"
tx2gene_file <- "tx2gene_gencode_v50.csv"

# Read tx2gene mapping
tx2gene <- fread(tx2gene_file)

# Find all Salmon quant files
quant_files <- list.files(
    salmon_dir,
    pattern = "^quant\\.sf$",
    recursive = TRUE,
    full.names = TRUE
)
quant_files <- sort(quant_files)

# Extract sample names
sample_names <- basename(dirname(quant_files))
sample_names <- sub("_quant$", "", sample_names)

# Create metadata
metadata <- data.frame(sample = sample_names, stringsAsFactors = FALSE)

# Extract patient, sex, condition from sample names
# Expected format: PatientXX_M_C or PatientXX_F_NC
metadata$patient <- sub("^Patient([0-9]+)[FM](C|NC)$", "Patient\\1", metadata$sample)
metadata$sex <- sub("^Patient[0-9]+([FM])(C|NC)$", "\\1", metadata$sample)
metadata$condition <- sub("^Patient[0-9]+[FM](C|NC)$", "\\1", metadata$sample)

metadata$condition <- factor(metadata$condition, levels = c("NC", "C"))
metadata$patient <- factor(metadata$patient)
metadata$sex <- factor(metadata$sex)

# Name quant files
names(quant_files) <- metadata$sample

# Import Salmon quantification
txi <- tximport(
    quant_files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "lengthScaledTPM"
)

# Save outputs
write.csv(txi$counts, "gene_counts_lengthScaledTPM.csv", quote = FALSE)
write.csv(metadata, "sample_metadata.csv", row.names = FALSE, quote = FALSE)

saveRDS(txi, "txi_gencode_v50.rds")
saveRDS(metadata, "sample_metadata.rds")

cat("✓ tximport completed. Genes:", nrow(txi$counts), "Samples:", ncol(txi$counts), "\n")
