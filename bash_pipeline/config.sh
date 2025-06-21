#!/bin/bash

# === Global Configuration ===

# Number of CPU threads to use
export THREADS=15

# Read type: SE = Single-End, PE = Paired-End
export MODE="SE"

# Paths to genome and reads
export GENOME_DIR="data/genome"
export READ_DIR="data/Trim"

# Genome annotation and sequence files
export GTF="$GENOME_DIR"
export FASTA="$GENOME_DIR"

# Output directories (used in other scripts)
export BOWTIE_DIR="Mapping/Bowtie"
export HISAT2_DIR="Mapping/Hisat2"
export STAR_DIR="Mapping/STAR"

export QUANT_BOWTIE="Quantification/Bowtie"
export QUANT_HISAT2="Quantification/Hisat2"
export QUANT_STAR="Quantification/STAR"

export DEG_DIR="DEG"

# Logs (optional)
export LOG_DIR="logs"
mkdir -p "$LOG_DIR"
