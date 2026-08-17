#!/bin/bash

SAMPLE_FILE="samples.txt"
SALMON_INDEX="salmon_index_v50"
OUTDIR="salmon"

mkdir -p "$OUTDIR"

while read -r sample; do
    [[ -z "$sample" ]] && continue
    
    echo "Processing: $sample"
    salmon quant \
        -i "$SALMON_INDEX" \
        -l A \
        -1 "trimmed_fastq/${sample}_R1_val_1.fq.gz" \
        -2 "trimmed_fastq/${sample}_R2_val_2.fq.gz" \
        -o "$OUTDIR/${sample}_quant" \
        --seqBias \
        --gcBias \
        --validateMappings \
        -p 8

done < "$SAMPLE_FILE"

echo "Salmon quantification completed."
