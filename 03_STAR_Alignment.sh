#!/bin/bash
# Align trimmed reads with STAR

set -euo pipefail
ulimit -n 65535

GENOME_DIR="STAR_genomeIndices_50"
WORKDIR="./work"
SAMPLE_FILE="samples.txt"

THREADS=32
MM_NUM=10

mkdir -p "$WORKDIR"

# Verify STAR index exists
if [ ! -f "$GENOME_DIR/SA" ]; then
    echo "ERROR: STAR index not found."
    exit 1
fi

while read -r sample; do
    [[ -z "$sample" ]] && continue
    
    echo "============================================================"
    echo "Processing sample: $sample"
    echo "============================================================"
    
    r1="${sample}_R1_val_1.fq.gz"
    r2="${sample}_R2_val_2.fq.gz"
    prefix="${sample}_mm${MM_NUM}_"
    out_bam="${prefix}Aligned.sortedByCoord.out.bam"
    
    # Check if alignment already exists
    if [ -f "$WORKDIR/$out_bam" ]; then
        echo "Skipping $sample — alignment already exists."
        continue
    fi
    
    # Check if trimmed files exist
    if [ ! -f "trimmed_fastq/$r1" ] || [ ! -f "trimmed_fastq/$r2" ]; then
        echo "ERROR: Trimmed FASTQ files not found."
        exit 1
    fi
    
    # Run STAR alignment
    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "trimmed_fastq/$r1" "trimmed_fastq/$r2" \
        --outFileNamePrefix "$WORKDIR/${prefix}" \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand zcat \
        --limitBAMsortRAM 100000000000 \
        --outFilterMultimapNmax "$MM_NUM"
    
    echo "Completed: $sample"
    
done < "$SAMPLE_FILE"

echo "All samples aligned successfully."
