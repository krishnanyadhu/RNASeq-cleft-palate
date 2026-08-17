#!/bin/bash

STAR \
    --runThreadN 40 \
    --runMode genomeGenerate \
    --genomeDir STAR_genomeIndices_50 \
    --genomeFastaFiles reference/GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile reference/gencode.v50.annotation.gtf \
    --sjdbOverhang 149

# Verify index
ls -lh STAR_genomeIndices_50/
