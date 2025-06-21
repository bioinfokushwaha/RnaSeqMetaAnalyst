#!/bin/bash
# Set global variables
export THREADS=15
export MODE="SE"
export GENOME_DIR="data/genome"
export READ_DIR="data/Trim"
export GTF="data/genome"
export FASTA="data/genome"

set -e  # Stop if any script fails
set -x  # Print each command for debugging

mkdir -p DEG Indices/Bowtie Mapping/Bowtie Quantification/Bowtie/{HT,FC,RSEM}
mkdir -p DEG Indices/Hisat2 Mapping/Hisat2 Quantification/Hisat2/{HT,FC}
mkdir -p DEG Indices/STAR Mapping/STAR Quantification/STAR/{HT,FC,RSEM}
# === Load Config ===
source config.sh

# === Step 1: Quality Control ===
./quality_control.sh

# === Step 3: Alignment ===
./hisat2.sh
./bowtie.sh
./star.sh

# == === cleaning_scripts ===
./DEG/star_clean.sh
./DEG/bowtie_clean.sh
./DEG/hisat2_clean.sh

# === Step 6: Differential Expression Analysis ===

./run_all_edgeR.sh
./run_all_deseq2.sh

# === Step 7: R-based Analysis ===
#Rscript run_deseq2_analysis.R
#Rscript run_deseq222_analysis.R
#3Rscript run_edgeR_analysis.R
