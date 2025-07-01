#!/bin/bash
# Set global variables
#!/bin/bash

set -e
set -x

# Use environment variables (passed via Docker)
echo "THREADS: $THREADS"
echo "MODE: $MODE"
echo "GENOME_DIR: $GENOME_DIR"
echo "READ_DIR: $READ_DIR"
echo "GTF: $GTF"
echo "FASTA: $FASTA"

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
./DEG/hisat2.sh
./DEG/bowtie.sh
./DEG/star.sh

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
