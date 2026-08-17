#!/bin/bash
# Trim raw FASTQ reads with Trim Galore

SAMPLE_FILE="samples.txt"
OUTDIR="trimmed_fastq"

mkdir -p "$OUTDIR"

while read -r sample; do
    [[ -z "$sample" ]] && continue
    
    echo "=========================================="
    echo "Processing: $sample"
    echo "=========================================="

    trim_galore \
        --paired \
        --illumina \
        --phred33 \
        --quality 20 \
        --stringency 3 \
        --length 36 \
        --cores 4 \
        --gzip \
        --fastqc \
        --output_dir "$OUTDIR" \
        "${sample}_R1.fastq.gz" \
        "${sample}_R2.fastq.gz" \
        2>&1 | tee "${OUTDIR}/${sample}_trim_galore.log"

done < "$SAMPLE_FILE"

echo "All samples completed."