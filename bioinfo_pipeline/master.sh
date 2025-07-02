#!/bin/bash


set -e  # Exit on any error
set -x  # Print each command before executing (for debugging)

# === Show Environment Variables (set via Docker `-e`) ===
echo "THREADS: $THREADS"
echo "MODE: $MODE"
echo "GENOME_DIR: $GENOME_DIR"
echo "READ_DIR: $READ_DIR"
echo "GTF: $GTF"
echo "FASTA: $FASTA"

# === Export Variables for Subscripts ===
export THREADS
export MODE
export GENOME_DIR
export READ_DIR
export GTF
export FASTA

mkdir -p DEG Indices/Bowtie Mapping/Bowtie Quantification/Bowtie/{HT,FC,RSEM}
mkdir -p DEG Indices/Hisat2 Mapping/Hisat2 Quantification/Hisat2/{HT,FC}
mkdir -p DEG Indices/STAR Mapping/STAR Quantification/STAR/{HT,FC,RSEM}


# === Step 1: Quality Control ===
./quality_control.sh

# === Step 3: Alignment ===
./hisat2.sh
./bowtie.sh
./star.sh

# == === cleaning_scripts ===
./star_clean.sh
./bowtie_clean.sh
./hisat2_clean.sh

# === Step 6: Differential Expression Analysis ===

./run_all_edgeR.sh
./run_all_deseq2.sh

# === Step 7: R-based Analysis ===
#Rscript run_deseq2_analysis.R
#Rscript run_deseq222_analysis.R
#3Rscript run_edgeR_analysis.R
