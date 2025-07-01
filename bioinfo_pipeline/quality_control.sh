#!/bin/bash

# === INPUT FROM ENV VARIABLES PASSED BY DOCKER ===
RAW_DIR="$READ_DIR"                  # e.g., /data/raw
TRIM_DIR="/data/Trim"               # Where trimmed reads go

# === QC OUTPUT PATHS ===
FASTQC_BEFORE="/opt/bioinfo/qc/fastqc_raw"
FASTQC_AFTER="/opt/bioinfo/qc/fastqc_clean"
FASTP_REPORTS="/opt/bioinfo/qc/fastp_reports"
MULTIQC_DIR="/opt/bioinfo/qc/multiqc"

# === CREATE OUTPUT DIRECTORIES ===
mkdir -p "$TRIM_DIR" "$FASTQC_BEFORE" "$FASTQC_AFTER" "$FASTP_REPORTS" "$MULTIQC_DIR"

echo "📦 MODE: $MODE"
echo "📂 RAW input: $RAW_DIR"
echo "📂 Trimmed output: $TRIM_DIR"
echo "💻 Threads: $THREADS"

# === QC PIPELINE ===
if [[ "$MODE" == "PE" ]]; then
    echo "🔄 Running in Paired-End mode..."

    for R1 in "$RAW_DIR"/*_R1.fastq.gz; do
        [[ ! -e "$R1" ]] && echo "⚠️ No *_R1.fastq.gz files in $RAW_DIR" && continue

        SAMPLE=$(basename "$R1" _R1.fastq.gz)
        R2="$RAW_DIR/${SAMPLE}_R2.fastq.gz"

        if [[ ! -f "$R2" ]]; then
            echo "❌ Missing R2 file for $SAMPLE"
            continue
        fi

        echo "🧪 Processing sample: $SAMPLE"

        # Run FastQC BEFORE
        fastqc -t "$THREADS" -o "$FASTQC_BEFORE" "$R1" "$R2"

        # Run fastp
        fastp -w "$THREADS" \
            -i "$R1" -I "$R2" \
            -o "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            -O "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        # Run FastQC AFTER
        fastqc -t "$THREADS" -o "$FASTQC_AFTER" \
            "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz"

        echo "✅ Done with sample: $SAMPLE"
    done

else
    echo "🔄 Running in Single-End mode..."

    for R in "$RAW_DIR"/*.fastq.gz; do
        [[ ! -e "$R" ]] && echo "⚠️ No *.fastq.gz files in $RAW_DIR" && continue

        # Skip PE-style names just in case
        [[ "$R" == *_R1.fastq.gz || "$R" == *_R2.fastq.gz ]] && continue

        SAMPLE=$(basename "$R" .fastq.gz)

        echo "🧪 Processing sample: $SAMPLE"

        # Run FastQC BEFORE
        fastqc -t "$THREADS" -o "$FASTQC_BEFORE" "$R"

        # Run fastp
        fastp -w "$THREADS" \
            -i "$R" \
            -o "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        # Run FastQC AFTER
        fastqc -t "$THREADS" -o "$FASTQC_AFTER" \
            "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz"

        echo "✅ Done with sample: $SAMPLE"
    done
fi

# === Run MultiQC on all QC results ===
echo "📊 Running MultiQC..."
multiqc /opt/bioinfo/qc/ -o "$MULTIQC_DIR"
echo "📁 Final MultiQC report: $MULTIQC_DIR"

# === Final Rename and Decompress ===
echo "🧼 Finalizing files: renaming and decompressing..."

for FILE in "$TRIM_DIR"/*_trimmed.fastq.gz; do
    [[ ! -f "$FILE" ]] && continue

    BASENAME=$(basename "$FILE" _trimmed.fastq.gz)
    RENAMED="$TRIM_DIR/${BASENAME}.fastq.gz"
    FINAL="$TRIM_DIR/${BASENAME}.fastq"

    mv "$FILE" "$RENAMED"

    if [[ ! -f "$FINAL" ]]; then
        echo "🗜️ Decompressing: $RENAMED"
        gunzip "$RENAMED"
    else
        echo "⚠️ Already exists: $FINAL — skipping decompression"
    fi
done

echo "✅ All files renamed and decompressed successfully."
