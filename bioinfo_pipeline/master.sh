#!/bin/bash

set -e  # Exit on any error
set -x  # Print each command

# === Step 0: Echo Environment Variables Passed from Docker ===
echo "[INFO] THREADS:       $THREADS"
echo "[INFO] MODE:          $MODE"
echo "[INFO] PROJECT_NAME:  $PROJECT_NAME"
echo "[INFO] GENOME_DIR:    $GENOME_DIR"
echo "[INFO] READ_DIR:      $READ_DIR"
echo "[INFO] GTF:           $GTF"
echo "[INFO] FASTA:         $FASTA"

# === Step 1: Export Variables for Subscripts ===
export THREADS MODE PROJECT_NAME GENOME_DIR READ_DIR GTF FASTA

# === Step 2: Define Output Directory Base ===
OUTPUT_BASE="/data/${PROJECT_NAME}/output"
export OUTPUT_BASE

# === Step 3: Create Directory Structure ===
mkdir -p \
    "$OUTPUT_BASE/DEG" \
    "$OUTPUT_BASE/Indices/Bowtie" \
    "$OUTPUT_BASE/Mapping/Bowtie" \
    "$OUTPUT_BASE/Quantification/Bowtie/HT" \
    "$OUTPUT_BASE/Quantification/Bowtie/FC" \
    "$OUTPUT_BASE/Quantification/Bowtie/RSEM" \
    "$OUTPUT_BASE/Indices/Hisat2" \
    "$OUTPUT_BASE/Mapping/Hisat2" \
    "$OUTPUT_BASE/Quantification/Hisat2/HT" \
    "$OUTPUT_BASE/Quantification/Hisat2/FC" \
    "$OUTPUT_BASE/Quantification/Hisat2/RSEM" \
    "$OUTPUT_BASE/Indices/STAR" \
    "$OUTPUT_BASE/Mapping/STAR" \
    "$OUTPUT_BASE/Quantification/STAR/HT" \
    "$OUTPUT_BASE/Quantification/STAR/FC" \
    "$OUTPUT_BASE/Quantification/STAR/RSEM"

# === Step 4: Quality Control ===
[ -x ./quality_control.sh ] && ./quality_control.sh || echo "[SKIP] quality_control.sh not found or not executable"

# === Step 5: Alignment Steps ===
[ -x ./hisat2.sh ] && ./hisat2.sh || echo "[SKIP] hisat2.sh not found or not executable"
[ -x ./bowtie.sh ] && ./bowtie.sh || echo "[SKIP] bowtie.sh not found or not executable"
[ -x ./star.sh ] && ./star.sh || echo "[SKIP] star.sh not found or not executable"

# === Step 6: Clean Intermediate Files (Optional) ===
[ -x ./hisat2_clean.sh ] && ./hisat2_clean.sh || echo "[SKIP] hisat2_clean.sh not found"
[ -x ./bowtie_clean.sh ] && ./bowtie_clean.sh || echo "[SKIP] bowtie_clean.sh not found"
[ -x ./star_clean.sh ] && ./star_clean.sh || echo "[SKIP] star_clean.sh not found"

# === Step 7: Differential Expression Analysis ===
[ -x ./run_all_edgeR.sh ] && ./run_all_edgeR.sh || echo "[SKIP] run_all_edgeR.sh not found"
[ -x ./run_all_deseq2.sh ] && ./run_all_deseq2.sh || echo "[SKIP] run_all_deseq2.sh not found"

# === Step 8: Optional R-based Analysis Scripts (Commented) ===
# Uncomment below if you want to run these separately
# Rscript run_deseq2_analysis.R
# Rscript run_deseq222_analysis.R
# Rscript run_edgeR_analysis.R

# === Done ===
echo "[DONE] Pipeline completed for project: $PROJECT_NAME"
echo "[OUTPUT] Results located at: $OUTPUT_BASE"
