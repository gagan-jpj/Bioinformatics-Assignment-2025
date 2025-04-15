#!/bin/bash

# ----------------------------------
# NGS Pipeline Script for Assignment
# ----------------------------------

# Set up directories
mkdir -p ~/ngs_assignement/dnaseq/{data/{untrimmed_fastq,trimmed_fastq,reference,aligned_data,qc_reports},annotation,vcf_files}

# ----------------------------------
# 1. Install Miniconda or Anaconda
# ----------------------------------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

# Optional: Install Anaconda
#wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
#chmod +x Anaconda3-2022.10-Linux-x86_64.sh
#bash Anaconda3-2022.10-Linux-x86_64.sh

# ----------------------------------
# 2. Install Tools
# ----------------------------------
conda create -n ngs_env -y
source ~/miniconda3/bin/activate ngs_env
conda install -y samtools bwa freebayes picard bedtools trimmomatic fastqc vcftools

# ----------------------------------
# 3. Download Input Files
# ----------------------------------
# FASTQ files
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz -P ~/ngs_assignement/dnaseq/data/untrimmed_fastq/
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz -P ~/ngs_assignement/dnaseq/data/untrimmed_fastq/

# Annotation and Reference
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed -P ~/ngs_assignement/dnaseq/data/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -P ~/ngs_assignement/dnaseq/data/reference/
gunzip ~/ngs_assignement/dnaseq/data/reference/hg19.fa.gz
bwa index ~/ngs_assignement/dnaseq/data/reference/hg19.fa

# ----------------------------------
# 4. Quality Control (Pre- and Post-Trimming)
# ----------------------------------
fastqc -o ~/ngs_assignement/dnaseq/data/qc_reports ~/ngs_assignement/dnaseq/data/untrimmed_fastq/*.fastq.qz

trimmomatic PE -threads 4 -phred33 \
~/ngs_assignement/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.qz \
~/ngs_assignement/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.qz \
-baseout ~/ngs_assignement/dnaseq/data/trimmed_fastq/NGS0001_trimmed_ \
ILLUMINACLIP:/path/to/TruSeq3-PE.fa:2:30:10 TRAILING:25 MINLEN:50

fastqc -o ~/ngs_assignement/dnaseq/data/qc_reports ~/ngs_assignement/dnaseq/data/trimmed_fastq/*1P ~/ngs_assignement/dnaseq/data/trimmed_fastq/*2P

# ----------------------------------
# 5. Alignment
# ----------------------------------
bwa mem -t 4 -R '@RG\tID:sample1\tSM:NGS0001\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
~/ngs_assignement/dnaseq/data/reference/hg19.fa \
~/ngs_assignement/dnaseq/data/trimmed_fastq/NGS0001_trimmed__1P \
~/ngs_assignement/dnaseq/data/trimmed_fastq/NGS0001_trimmed__2P \
> ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_aligned.sam

# Convert and sort BAM files
samtools view -S -b ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_aligned.sam > ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_unsorted.bam
samtools sort -n ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_unsorted.bam -o ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_namesorted.bam
samtools fixmate -m ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_namesorted.bam ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_fixmate.bam
samtools sort ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_fixmate.bam -o ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_sorted.bam
samtools markdup ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_sorted.bam ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_dedup.bam

# Filter BAM by quality
samtools view -b -q 30 ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_dedup.bam > ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_filtered.bam

# Generate statistics
samtools flagstat ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_filtered.bam > ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_flagstats.txt
samtools index ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_filtered.bam
samtools idxstats ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_filtered.bam > ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_idxstats.txt
samtools depth ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_filtered.bam > ~/ngs_assignement/dnaseq/data/aligned_data/coverage.txt
samtools stats ~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_filtered.bam > ~/ngs_assignement/dnaseq/data/aligned_data/stats.txt

# ----------------------------------
# 6. Variant Calling and Filtering
# ----------------------------------
freebayes -f ~/ngs_assignement/dnaseq/data/reference/hg19.fa \
-t ~/ngs_assignement/dnaseq/data/annotation.bed \
~/ngs_assignement/dnaseq/data/aligned_data/NGS0001_filtered.bam \
> ~/ngs_assignement/dnaseq/data/vcf_files/NGS0001_raw.vcf

vcftools --vcf ~/ngs_assignement/dnaseq/data/vcf_files/NGS0001_raw.vcf \
--minQ 30 --recode --out ~/ngs_assignement/dnaseq/data/vcf_files/NGS0001_filtered

# ----------------------------------
# 7. Variant Annotation with ANNOVAR
# ----------------------------------
# ./convert2annovar.pl -format vcf4 NGS0001_filtered.vcf -outfile NGS0001_filtered.avinput
# ./table_annovar.pl NGS0001_filtered.avinput humandb/ -buildver hg19 -out NGS0001_annotated \
# -remove -protocol refGene,ensGene,clinvar_20180603,exac03 -operation g,g,f,f -nastring . -vcfinput

# ----------------------------------
# 8. Annotation with snpEff
# ----------------------------------
# wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# unzip snpEff_latest_core.zip && cd snpEff
# java -Xmx4g -jar snpEff.jar download -v GRCh37.75
# java -Xmx4g -jar snpEff.jar -v GRCh37.75 NGS0001_filtered.vcf > snpeff_annotated.vcf

# ----------------------------------
# 9. Variant Prioritization
# ----------------------------------
grep "exonic" NGS0001_annotated.hg19_multianno.txt > exonic_variants.txt
awk '{ if ($7 == "novel") print $0 }' exonic_variants.txt > novel_exonic_variants.txt

# END OF SCRIPT

