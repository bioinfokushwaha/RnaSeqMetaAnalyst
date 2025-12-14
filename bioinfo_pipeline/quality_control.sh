#!/bin/bash

# === CONFIGURATION ===

READ_DIR=$READ_DIR
TRIM_DIR=$TRIM_DIR          # Final trimmed reads for downstream
FASTQC_BEFORE="qc/fastqc_raw"
FASTQC_AFTER="qc/fastqc_clean"
FASTP_REPORTS="qc/fastp_reports"
MULTIQC_DIR="qc/multiqc"

# === CREATE OUTPUT DIRECTORIES ===
mkdir -p "$TRIM_DIR" "$FASTQC_BEFORE" "$FASTQC_AFTER" "$FASTP_REPORTS" "$MULTIQC_DIR"

echo "Processing mode: $MODE"
echo "Input directory: $READ_DIR"
echo "Output (trimmed reads): $TRIM_DIR"

if [[ $MODE == "PE" ]]; then
    echo ">>> Running in Paired-End mode..."
    for R1 in "$READ_DIR"/*_1.fastq.gz; do
        SAMPLE=$(basename "$R1" _1.fastq.gz)
        R2="$READ_DIR/${SAMPLE}_2.fastq.gz"  # FIXED: was $RAW_DIR

        echo "Processing sample: $SAMPLE"

        # FastQC BEFORE
        fastqc -t $THREADS -o "$FASTQC_BEFORE" "$R1" "$R2"

        # Run fastp
        fastp -w $THREADS \
            -i "$R1" -I "$R2" \
            -o "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            -O "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        # FastQC AFTER
        fastqc -t $THREADS -o "$FASTQC_AFTER" \
            "$TRIM_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
            "$TRIM_DIR/${SAMPLE}_R2_trimmed.fastq.gz"

        echo "--- Done with sample: $SAMPLE ---"
    done

else
    echo ">>> Running in Single-End mode..."
    for R in "$READ_DIR"/*.fastq.gz; do
        if [[ "$R" == *_1.fastq.gz || "$R" == *_2.fastq.gz ]]; then
            continue  # skip paired-end reads in SE mode
        fi

        SAMPLE=$(basename "$R" .fastq.gz)
        echo "Processing sample: $SAMPLE"

        # FastQC BEFORE
        fastqc -t $THREADS -o "$FASTQC_BEFORE" "$R"

        # Run fastp
        fastp -w $THREADS \
            -i "$R" \
            -o "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz" \
            -h "$FASTP_REPORTS/${SAMPLE}_fastp.html" \
            -j "$FASTP_REPORTS/${SAMPLE}_fastp.json"

        # FastQC AFTER
        fastqc -t $THREADS -o "$FASTQC_AFTER" \
            "$TRIM_DIR/${SAMPLE}_trimmed.fastq.gz"

        echo "--- Done with sample: $SAMPLE ---"
    done
fi


# === Run MultiQC ===
echo "Running MultiQC..."
multiqc qc/ -o "$MULTIQC_DIR"

echo "fastp, fastqc, and MultiQC preprocessing complete."
echo "Final MultiQC report available at: $MULTIQC_DIR"

# === Rename and Decompress Final Files ===
echo "Finalizing files: renaming and decompressing..."

if [[ $MODE == "PE" ]]; then
    for FILE in "$TRIM_DIR"/*_R1_trimmed.fastq.gz; do
        SAMPLE=$(basename "$FILE" _R1_trimmed.fastq.gz)

        for READ in R1 R2; do
            OLD="$TRIM_DIR/${SAMPLE}_${READ}_trimmed.fastq.gz"
            NEW="$TRIM_DIR/${SAMPLE}_${READ}.fastq.gz"
            FINAL="$TRIM_DIR/${SAMPLE}_${READ}.fastq"

            mv "$OLD" "$NEW"

            if [[ ! -f "$FINAL" ]]; then
                echo "Decompressing: $NEW"
                gunzip -k "$NEW"
            else
                echo "Skipping decompression for $NEW — $FINAL already exists."
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
            gunzip -d -k "$NEWNAME"
        else
            echo "Skipping decompression for $NEWNAME — $FINAL already exists."
        fi
    done
fi

echo "All files renamed and decompressed successfully."
