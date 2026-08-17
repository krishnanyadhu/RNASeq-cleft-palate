# RNA-seq Analysis Pipeline — Cleft Palate

An RNA-seq analysis workflow for paired cleft palate and non-cleft control samples, covering read preprocessing, genome alignment, transcript quantification, quality control, differential expression, functional enrichment, and visualization.
**Analysis flow**

`FASTQ → Read trimming/QC → STAR genome index → STAR alignment → Salmon quantification → tximport → QC → DESeq2 → LFC shrinkage → DEG annotation → ORA enrichment → Visualization`

The differential expression analysis uses a paired study design:

`~ patient + condition`

where `condition` represents cleft (`C`) versus non-cleft (`NC`) samples.

---

