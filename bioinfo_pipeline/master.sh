#!/bin/bash

set -e  # Exit on any error
set -x  # Print commands for debugging

# === Validate Required Environment Variables ===
: "${THREADS:?THREADS is not set}"
: "${MODE:?MODE is not set (should be SE or PE)}"
: "${PROJECT_NAME:?PROJECT_NAME is not set}"
: "${GENOME_DIR:?GENOME_DIR is not set}"
: "${READ_DIR:?READ_DIR is not set}"
: "${GTF:?GTF is not set}"
: "${FASTA:?FASTA is not set}"

# === Normalize MODE and Export Variables ===
MODE=$(echo "$MODE" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')  # Normalize to SE or PE
export THREADS MODE PROJECT_NAME GENOME_DIR READ_DIR GTF FASTA

# === Logging ===
echo "=== RNA-Seq Master Pipeline Started at $(date) ==="
echo "THREADS: $THREADS"
echo "MODE: $MODE"
echo "PROJECT_NAME: $PROJECT_NAME"
echo "GENOME_DIR: $GENOME_DIR"
echo "READ_DIR: $READ_DIR"
echo "GTF: $GTF"
echo "FASTA: $FASTA"

# === Output Directory Base ===
OUTPUT_BASE="/data/${PROJECT_NAME}/output"
export OUTPUT_BASE  # Export so subscripts can use it

# === Create Output Directory Structure ===
mkdir -p "$OUTPUT_BASE/DEG" \
         "$OUTPUT_BASE/Indices/Bowtie" \
         "$OUTPUT_BASE/Mapping/Bowtie" \
         "$OUTPUT_BASE/Quantification/Bowtie/HT" \
         "$OUTPUT_BASE/Quantification/Bowtie/FC" \
         "$OUTPUT_BASE/Quantification/Bowtie/RSEM" \
         "$OUTPUT_BASE/Indices/Hisat2" \
         "$OUTPUT_BASE/Mapping/Hisat2" \
         "$OUTPUT_BASE/Quantification/Hisat2/HT" \
         "$OUTPUT_BASE/Quantification/Hisat2/FC" \
         "$OUTPUT_BASE/Indices/STAR" \
         "$OUTPUT_BASE/Mapping/STAR" \
         "$OUTPUT_BASE/Quantification/STAR/HT" \
         "$OUTPUT_BASE/Quantification/STAR/FC" \
         "$OUTPUT_BASE/Quantification/STAR/RSEM"

# === Ensure Subscripts Are Executable ===
chmod +x *.sh

# === Run the Pipeline Steps ===
micromamba activate rnaseq-cli
# Step 1: Quality Control
./quality_control.sh

# Step 2: Alignment
./hisat2.sh
./bowtie.sh
./star.sh

# Step 3: Clean intermediate files (optional)
./star_clean.sh
./bowtie_clean.sh
./hisat2_clean.sh

micromamba activate rnaseq-r
# Step 4: Differential Expression Analysis
./run_all_edgeR.sh
./run_all_deseq2.sh

# Optional: R-based additional analysis (manually uncomment if needed)
# Rscript run_deseq2_analysis.R
# Rscript run_deseq222_analysis.R
# Rscript run_edgeR_analysis.R

# === Done ===
echo "=== RNA-Seq Master Pipeline Finished at $(date) ==="
