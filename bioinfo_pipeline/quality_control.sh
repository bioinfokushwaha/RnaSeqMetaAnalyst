RAW_DIR="data/raw"
READ_DIR="data/Trim"           # Final trimmed reads for downstream
FASTQC_BEFORE="qc/fastqc_raw"
FASTQC_AFTER="qc/fastqc_clean"
FASTP_REPORTS="qc/fastp_reports"
MULTIQC_DIR="qc/multiqc"

# === CREATE OUTPUT DIRECTORIES ===
mkdir -p "$READ_DIR" "$FASTQC_BEFORE" "$FASTQC_AFTER" "$FASTP_REPORTS" "$MULTIQC_DIR"

echo "Processing mode: $MODE"
echo "Input directory: $RAW_DIR"
echo "Output (trimmed reads): $READ_DIR"

if [[ $MODE == "PE" ]]; then
    echo ">>> Running in Paired-End mode..."
    for R1 in "$RAW_DIR"/*_R1.fastq.gz; do
        SAMPLE=$(basename "$R1" _R1.fastq.gz)
        R2="$RAW_DIR/${SAMPLE}_R2.fastq.gz"

        echo "Processing sample: $SAMPLE"

        # FastQC BEFORE
        fastqc -t $THREADS -o "$FASTQC_BEFORE" "$R1" "$R2"

        # Run fastp
        fastp -w $THREADS \
            -i "$R1" -I "$R2" \
            -o "$READ_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            -O "$READ_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        # FastQC AFTER
        fastqc -t $THREADS -o "$FASTQC_AFTER" \
            "$READ_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            "$READ_DIR/${SAMPLE}_R2_trimmed.fastq.gz"

        echo "--- Done with sample: $SAMPLE ---"
    done

else
    echo ">>> Running in Single-End mode..."
    for R in "$RAW_DIR"/*.fastq.gz; do
        if [[ "$R" == *_R.fastq.gz || "$R" == *_R.fastq.gz ]]; then
            continue  # skip paired-end reads in SE mode
        fi

        SAMPLE=$(basename "$R" .fastq.gz)

        echo "Processing sample: $SAMPLE"

        # FastQC BEFORE
        fastqc -t $THREADS -o "$FASTQC_BEFORE" "$R"

        # Run fastp
        fastp -w $THREADS \
            -i "$R" \
            -o "$READ_DIR/${SAMPLE}_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        # FastQC AFTER
        fastqc -t $THREADS -o "$FASTQC_AFTER" \
            "$READ_DIR/${SAMPLE}_trimmed.fastq.gz"

        echo "--- Done with sample: $SAMPLE ---"
    done
fi

# === Run MultiQC ===
echo "Running MultiQC..."
multiqc qc/ -o "$MULTIQC_DIR"

echo "✅ fastp, fastqc, and MultiQC preprocessing complete."
echo "📁 Final MultiQC report available at: $MULTIQC_DIR"
####
# === Rename and Decompress Final Files ===
echo "🧼 Finalizing files: renaming and decompressing..."

for FILE in "$READ_DIR"/*_trimmed.fastq.gz; do
    BASENAME=$(basename "$FILE" _trimmed.fastq.gz)
    TEMP_NAME="$READ_DIR/${BASENAME}.fastq.gz"
    FINAL_NAME="$READ_DIR/${BASENAME}.fastq"

#    # Rename trimmed file to remove "_trimmed"
    mv "$FILE" "$TEMP_NAME"

    # Decompress only if final .fastq doesn't already exist
    if [[ ! -f "$FINAL_NAME" ]]; then
        echo "Decompressing: $TEMP_NAME"
        gunzip "$TEMP_NAME"
    else
        echo "⚠️  Skipping decompression for $TEMP_NAME — $FINAL_NAME already exists."
    fi
done

echo "✅ All files renamed and decompressed successfully."


