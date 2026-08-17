# RNA-seq Analysis Pipeline — Cleft Palate

An RNA-seq analysis workflow for paired cleft palate and non-cleft control samples, covering read preprocessing, genome alignment, transcript quantification, quality control, differential expression, functional enrichment, and visualization.
**Analysis flow**

`FASTQ → Read trimming/QC → STAR genome index → STAR alignment → Salmon quantification → tximport → QC → DESeq2 → LFC shrinkage → DEG annotation → ORA enrichment → Visualization`

The differential expression analysis uses a paired study design:

`~ patient + condition`

where `condition` represents cleft (`C`) versus non-cleft (`NC`) samples.

---

## Scripts

| Script | Description | Main output/location |
|---|---|---|
| `01_trimming.sh` | Performs paired-end adapter and quality trimming using Trim Galore and generates FastQC reports. | `trimmed_fastq/` |
| `02_STAR_index.sh` | Builds the STAR genome index from the GRCh38 reference genome and GENCODE v50 annotation. | `STAR_genomeIndices_50/` |
| `03_STAR_Alignment.sh` | Aligns trimmed reads to the STAR genome index and generates sorted BAM files, transcriptome alignments, gene counts, and unmapped reads. | `work/` |
| `04_Salmon_quant.sh` | Performs transcript-level quantification using Salmon with sequence- and GC-bias correction. | `salmon/` |
| `05_tx2gene.R` | Extracts transcript and gene identifiers from the GENCODE v50 GTF to create a transcript-to-gene mapping table. | `tx2gene_gencode_v50.csv` |
| `06_tximport.R` | Imports Salmon quantification, summarizes transcript-level estimates to gene-level counts, and generates sample metadata. | Root directory: `gene_counts_lengthScaledTPM.csv`, `sample_metadata.csv`, RDS objects |
| `07_QC.R` | Performs count filtering, library-size assessment, VST transformation, PCA, and sample correlation analysis. | Root directory: QC tables, plots, and RDS objects |
| `08_DESeq2.R` | Performs paired differential expression analysis comparing cleft and non-cleft samples using DESeq2. | Root directory: DESeq2 result tables, MA plot, and RDS objects |
| `09_DE_results.R` | Performs apeglm log2 fold-change shrinkage and generates DEG tables, volcano plot, and top-DEG heatmap. | Root directory: DEG tables, plots, and RDS objects |
| `10_annotate_DEGs.R` | Annotates DESeq2 results with gene names and gene types using GENCODE v50. | Root directory: annotated DEG tables and gene lists |
| `11_ORA_enrichment.R` | Performs GO, KEGG, and Reactome over-representation analysis for upregulated and downregulated genes. | `enrichment_ORA/` |
| `12_visualization.R` | Generates GO, KEGG, and Reactome enrichment plots. | Root directory: enrichment PDF plots |

---
## Tools and Package Versions


### Command-line tools

| Tool | Purpose | Version |
|---|---|---|
| Trim Galore | Adapter and quality trimming | **2.3.0** |
| FastQC | Read-level quality control | **0.12.1** |
| STAR | Genome indexing and read alignment | **2.7.11** |
| Salmon | Transcript quantification | **2.3.4** |


### R packages

| Package | Purpose | Version |
|---|---|---|
| `data.table` | Fast tabular data processing and GTF parsing | **1.18.4** |
| `tximport` | Import and summarize Salmon quantification | **1.30.0** |
| `DESeq2` | Differential expression analysis | **1.48.1** |
| `apeglm` | Log2 fold-change shrinkage | **1.30.0** |
| `ggplot2` | Data visualization | **4.0.2** |
| `pheatmap` | Heatmap generation | **1.0.13** |
| `clusterProfiler` | GO/KEGG enrichment analysis | **4.16.0** |
| `org.Hs.eg.db` | Human gene annotation | **3.21.0** |
| `ReactomePA` | Reactome pathway enrichment | **1.52.0** |
| `dplyr` | Data manipulation | **1.2.0** |
| `readr` | Tabular data import | **2.2.0** |
| `stringr` | String processing | **1.6.0** |

---

## Languages

- **Bash** — read preprocessing, STAR indexing/alignment, and Salmon quantification
- **R** — quantification import, QC, differential expression, annotation, enrichment analysis, and visualization

---

## Reference Data

The pipeline uses:

- **Genome:** GRCh38 primary assembly
- **Gene annotation:** GENCODE v50

---

## Statistical Analysis

Differential expression is performed using DESeq2 with a paired design:

```text
~ patient + condition
```

The reference level for `condition` is `NC`, and the primary contrast is:

```text
C vs NC
```

Significant DEGs are defined using:

- adjusted p-value (`FDR`) < 0.05
- absolute log2 fold change ≥ 1

Log2 fold-change shrinkage is subsequently performed using `apeglm`.

---

## Functional Enrichment

Over-representation analysis is performed separately for upregulated and downregulated genes.

The analysis includes:

- Gene Ontology Biological Process (GO-BP)
- Gene Ontology Molecular Function (GO-MF)
- Gene Ontology Cellular Component (GO-CC)
- KEGG pathways
- Reactome pathways

The background/universe is derived from the genes retained in the DESeq2 analysis.

Multiple-testing correction is performed using the Benjamini-Hochberg method.

---
