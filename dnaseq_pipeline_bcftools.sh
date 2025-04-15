#!/bin/bash

# Set paths
REF="reference.fa"
BAM="aligned_data/NGS0001_sorted.bam"
VCF_OUT="results/variants_bcftools.vcf"

# Step 1: Index the reference genome (only needed once)
samtools faidx $REF

# Step 2: Generate BCF using bcftools mpileup
bcftools mpileup -f $REF $BAM -Ou > results/raw.bcf

# Step 3: Call variants
bcftools call -mv -Ov -o $VCF_OUT results/raw.bcf

echo "Variant calling with bcftools completed!"
