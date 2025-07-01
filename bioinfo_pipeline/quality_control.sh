#!/bin/bash

set -e
set -x

# Enable nullglob so unmatched globs expand to nothing (avoids literal '*.fastq.gz')
shopt -s nullglob

# Docker-passed variables
RAW_DIR="$READ_DIR"
TRIM_DIR="data/Trim"
FASTQC_BEFORE="qc/fastqc_raw"
FASTQC_AFTER="qc/fastqc_clean"
FASTP_REPORTS="qc/fastp_reports"
MULTIQC_DIR="qc/multiqc"

# Create output directories
mkdir -p "$TRIM_DIR" "$FASTQC_BEFORE" "$FASTQC_AFTER" "$FASTP_REPORTS" "$MULTIQC_DIR"

echo "Processing mode: $MODE"
echo "Input directory: $RAW_DIR"
echo "Output (trimmed reads): $TRIM_DIR"

# === PE MODE ===
if [[ "$MODE" == "PE" ]]; then
    echo ">>> Running in Paired-End mode..."
    FILES=("$RAW_DIR"/*_R1.fastq.gz)
    if [[ ${#FILES[@]} -eq 0 ]]; then
        echo "❌ No paired-end (_R1.fastq.gz) files found in $RAW_DIR"
        exit 1
    fi

    for R1 in "${FILES[@]}"; do
        SAMPLE=$(basename "$R1" _R1.fastq.gz)
        R2="$RAW_DIR/${SAMPLE}_R2.fastq.gz"

        if [[ ! -f "$R2" ]]; then
            echo "❌ Missing R2 pair for sample $SAMPLE"
            continue
        fi

        echo "Processing sample: $SAMPLE"

        fastqc -t "$THREADS" -o "$FASTQC_BEFORE" "$R1" "$R2"

        fastp -w "$THREADS" \
            -i "$R1" -I "$R2" \
            -o "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            -O "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        fastqc -t "$THREADS" -o "$FASTQC_AFTER" \
            "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz"

        echo "--- Done with sample: $SAMPLE ---"
    done

# === SE MODE ===
else
    echo ">>> Running in Single-End mode..."
    FILES=("$RAW_DIR"/*.fastq.gz)
    if [[ ${#FILES[@]} -eq 0 ]]; then
        echo "❌ No single-end fastq.gz files found in $RAW_DIR"
        exit 1
    fi

    for R in "${FILES[@]}"; do
        BASENAME=$(basename "$R")
        if [[ "$BASENAME" == *_R1.fastq.gz || "$BASENAME" == *_R2.fastq.gz ]]; then
            continue  # Skip paired-end files
        fi

        SAMPLE=$(basename "$R" .fastq.gz)
        echo "Processing sample: $SAMPLE"

        fastqc -t "$THREADS" -o "$FASTQC_BEFORE" "$R"

        fastp -w "$THREADS" \
            -i "$R" \
            -o "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        fastqc -t "$THREADS" -o "$FASTQC_AFTER" "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz"

        echo "--- Done with sample: $SAMPLE ---"
    done
fi

# === Run MultiQC ===
echo "Running MultiQC..."
multiqc qc/ -o "$MULTIQC_DIR"
echo "✅ fastp, fastqc, and MultiQC preprocessing complete."
echo "📁 Final MultiQC report available at: $MULTIQC_DIR"

# === Rename & Decompress Final Files ===
echo "🧼 Finalizing files: renaming and decompressing..."

if [[ "$MODE" == "PE" ]]; then
    for FILE in "$TRIM_DIR"/*_R1_trimmed.fastq.gz; do
        SAMPLE=$(basename "$FILE" _R1_trimmed.fastq.gz)

        for READ in R1 R2; do
            OLD="$TRIM_DIR/${SAMPLE}_${READ}_trimmed.fastq.gz"
            NEW="$TRIM_DIR/${SAMPLE}_${READ}.fastq.gz"
            FINAL="$TRIM_DIR/${SAMPLE}_${READ}.fastq"

            mv "$OLD" "$NEW"

            if [[ ! -f "$FINAL" ]]; then
                echo "Decompressing: $NEW"
                gunzip "$NEW"
            else
                echo "⚠️  Skipping decompression for $NEW — $FINAL already exists."
            fi
        done
    done
else
    for FILE in "$TRIM_DIR"/*_trimmed.fastq.gz; do
        BASENAME=$(basename "$FILE" _trimmed.fastq.gz)
        NEWNAME="$TRIM_DIR/${BASENAME}.fastq.gz"
        FINAL="$TRIM_DIR/${BASENAME}.fastq"

        mv "$FILE" "$NEWNAME"

        if [[ ! -f "$FINAL" ]]; then
            echo "Decompressing: $NEWNAME"
            gunzip "$NEWNAME"
        else
            echo "⚠️  Skipping decompression for $NEWNAME — $FINAL already exists."
        fi
    done
fi
