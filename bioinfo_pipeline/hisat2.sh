#!/bin/bash
# ===================================================================
# HISAT2 COMPLETE PIPELINE: ALIGNMENT + QUANTIFICATION
# ===================================================================
# This script performs:
# 1. HISAT2 genome indexing (if needed)
# 2. HISAT2 alignment (PE/SE mode)
# 3. SAM to BAM conversion and sorting
# 4. featureCounts quantification
# 5. HTSeq quantification (with optimized piping)
#
# NOTE: HISAT2 is NOT compatible with RSEM (genome aligner, not transcriptome)
# ===================================================================

set -e  # Exit on error
set -o pipefail  # Exit on pipe failure

# === CONFIGURATION ===
ALIGNER="HISAT2"
INDEX_BASE="Indices/Hisat2/genome"
MAPPING_DIR="Mapping/Hisat2"
FC_OUTPUT_DIR="Quantification/Hisat2/FC"
HTSEQ_OUTPUT_DIR="Quantification/Hisat2/HT"

# === CREATE DIRECTORIES ===
mkdir -p "$MAPPING_DIR" "$FC_OUTPUT_DIR" "$HTSEQ_OUTPUT_DIR"

echo ""
echo "=========================================="
echo "HISAT2 COMPLETE PIPELINE"
echo "=========================================="
echo "Mode: $MODE"
echo "Threads: $THREADS"
echo "Genome: $FASTA"
echo "GTF: $GTF"
echo ""
echo "NOTE: HISAT2 does NOT support RSEM quantification"
echo "      (RSEM requires transcriptome alignments)"
echo ""

# ===================================================================
# STEP 1: HISAT2 INDEX PREPARATION
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 1: HISAT2 Index Preparation"
echo "=========================================="

# Check if we're using pre-built indices
if [ -n "$INDEX_DIR" ] && [ -d "$INDEX_DIR/Hisat2" ]; then
    echo ">>> Using pre-built HISAT2 index from: $INDEX_DIR/Hisat2"
    INDEX_BASE="$INDEX_DIR/Hisat2/genome"
    if [ ! -f "${INDEX_BASE}.1.ht2" ]; then
        echo "ERROR: HISAT2 index not found at ${INDEX_BASE}.1.ht2" >&2
        exit 1
    fi
else
    echo ">>> Building HISAT2 index locally..."
    mkdir -p "$(dirname "$INDEX_BASE")"
    
    if [ ! -f "${INDEX_BASE}.1.ht2" ]; then
        echo "  Creating HISAT2 index from: $FASTA"
        hisat2-build -p "$THREADS" "$FASTA" "$INDEX_BASE"
        echo "  ✓ HISAT2 index created successfully"
    else
        echo "  ✓ HISAT2 index already exists"
    fi
fi

# ===================================================================
# STEP 2: HISAT2 ALIGNMENT
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 2: HISAT2 Alignment"
echo "=========================================="

if [[ $MODE == "PE" ]]; then
    echo ">>> Running in Paired-End mode..."
    
    for R1 in "$TRIM_DIR"/*_R1.fastq; do
        if [ ! -f "$R1" ]; then
            continue  # Skip if no files match
        fi
        
        SAMPLE=$(basename "$R1" _R1.fastq)
        R2="$TRIM_DIR/${SAMPLE}_R2.fastq"
        
        if [ ! -f "$R2" ]; then
            echo "WARNING: R2 file not found for $SAMPLE: $R2"
            continue
        fi
        
        echo ""
        echo "  Processing sample: $SAMPLE"
        
        SAM_FILE="$MAPPING_DIR/${SAMPLE}.sam"
        BAM_FILE="$MAPPING_DIR/${SAMPLE}.bam"
        SORTED_BAM="$MAPPING_DIR/${SAMPLE}_sorted.bam"
        
        # Run HISAT2 alignment
        echo "    - Running HISAT2 alignment..."
        hisat2 -p "$THREADS" \
            -x "$INDEX_BASE" \
            -1 "$R1" \
            -2 "$R2" \
            -S "$SAM_FILE" \
            2> "$MAPPING_DIR/${SAMPLE}.hisat2.log"
        
        if [ ! -f "$SAM_FILE" ]; then
            echo "ERROR: HISAT2 failed for $SAMPLE" >&2
            exit 1
        fi
        
        # Convert SAM to BAM
        echo "    - Converting SAM to BAM..."
        samtools view -@ "$THREADS" -b -h "$SAM_FILE" > "$BAM_FILE"
        rm "$SAM_FILE"
        
        # Sort BAM by coordinate
        echo "    - Sorting BAM file..."
        samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM_FILE"
        rm "$BAM_FILE"
        
        # Index BAM
        echo "    - Indexing BAM file..."
        samtools index -@ "$THREADS" "$SORTED_BAM"
        
        echo "    ✓ $SAMPLE completed"
    done

else
    echo ">>> Running in Single-End mode..."
    
    for R in "$TRIM_DIR"/*.fastq; do
        if [ ! -f "$R" ]; then
            continue  # Skip if no files match
        fi
        
        # Skip paired-end files in SE mode
        if [[ "$R" == *_R1.fastq || "$R" == *_R2.fastq ]]; then
            continue
        fi
        
        SAMPLE=$(basename "$R" .fastq)
        
        echo ""
        echo "  Processing sample: $SAMPLE"
        
        SAM_FILE="$MAPPING_DIR/${SAMPLE}.sam"
        BAM_FILE="$MAPPING_DIR/${SAMPLE}.bam"
        SORTED_BAM="$MAPPING_DIR/${SAMPLE}_sorted.bam"
        
        # Run HISAT2 alignment
        echo "    - Running HISAT2 alignment..."
        hisat2 -p "$THREADS" \
            -x "$INDEX_BASE" \
            -U "$R" \
            -S "$SAM_FILE" \
            2> "$MAPPING_DIR/${SAMPLE}.hisat2.log"
        
        if [ ! -f "$SAM_FILE" ]; then
            echo "ERROR: HISAT2 failed for $SAMPLE" >&2
            exit 1
        fi
        
        # Convert SAM to BAM
        echo "    - Converting SAM to BAM..."
        samtools view -@ "$THREADS" -b -h "$SAM_FILE" > "$BAM_FILE"
        rm "$SAM_FILE"
        
        # Sort BAM by coordinate
        echo "    - Sorting BAM file..."
        samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM_FILE"
        rm "$BAM_FILE"
        
        # Index BAM
        echo "    - Indexing BAM file..."
        samtools index -@ "$THREADS" "$SORTED_BAM"
        
        echo "    ✓ $SAMPLE completed"
    done
fi

echo ""
echo "✓ HISAT2 alignment completed"

# Verify BAM files were created
BAM_COUNT=$(find "$MAPPING_DIR" -name "*_sorted.bam" -type f | wc -l)
if [ "$BAM_COUNT" -eq 0 ]; then
    echo "ERROR: No BAM files were created!" >&2
    exit 1
fi
echo "  Found $BAM_COUNT BAM files"

# ===================================================================
# STEP 3: FEATURECOUNTS QUANTIFICATION
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 3: featureCounts Quantification"
echo "=========================================="

echo ">>> Running featureCounts..."

FC_OUTPUT="$FC_OUTPUT_DIR/FC_Count.txt"

if [[ $MODE == "PE" ]]; then
    echo "  Mode: Paired-End"
    featureCounts -p -g gene_id -T "$THREADS" -O -M \
        -a "$GTF" \
        -o "$FC_OUTPUT" \
        "$MAPPING_DIR"/*_sorted.bam
else
    echo "  Mode: Single-End"
    featureCounts -g gene_id -T "$THREADS" -O -M \
        -a "$GTF" \
        -o "$FC_OUTPUT" \
        "$MAPPING_DIR"/*_sorted.bam
fi

if [ ! -f "$FC_OUTPUT" ]; then
    echo "ERROR: featureCounts failed" >&2
    exit 1
fi

FC_LINES=$(wc -l < "$FC_OUTPUT")
echo "✓ featureCounts completed ($FC_LINES lines)"

# ===================================================================
# STEP 4: HTSEQ QUANTIFICATION (OPTIMIZED WITH PIPES)
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 4: HTSeq Quantification"
echo "=========================================="

echo ">>> Processing BAM files for HTSeq..."

# Process each BAM file
for BAM in "$MAPPING_DIR"/*_sorted.bam; do
    if [ ! -f "$BAM" ]; then
        continue
    fi
    
    SAMPLE=$(basename "$BAM" _sorted.bam)
    HTSEQ_OUTPUT="$HTSEQ_OUTPUT_DIR/${SAMPLE}.htseq.txt"
    
    echo ""
    echo "  Sample: $SAMPLE"
    echo "    - Running HTSeq-count with optimized piping..."
    
    # Optimized approach: use pipes instead of creating intermediate files
    # This saves disk space and I/O operations
    samtools sort -n -u -@ "$THREADS" "$BAM" | \
    samtools view -h -@ "$THREADS" - | \
    htseq-count -f sam -r name -m union - "$GTF" \
        > "$HTSEQ_OUTPUT" 2> "$HTSEQ_OUTPUT_DIR/${SAMPLE}.htseq.err"
    
    if [ ! -f "$HTSEQ_OUTPUT" ]; then
        echo "    ERROR: HTSeq failed for $SAMPLE" >&2
        exit 1
    fi
    
    echo "    ✓ HTSeq completed for $SAMPLE"
done

echo ""
echo ">>> Merging HTSeq results into count matrix..."

# Get list of HTSeq output files
HTSEQ_FILES=("$HTSEQ_OUTPUT_DIR"/*.htseq.txt)

if [ ${#HTSEQ_FILES[@]} -eq 0 ] || [ ! -f "${HTSEQ_FILES[0]}" ]; then
    echo "ERROR: No HTSeq output files found" >&2
    exit 1
fi

# FIXED: Extract gene IDs (no header in HTSeq files, exclude __ summary lines)
echo "  - Extracting gene IDs..."
awk '!/^__/ {print $1}' "${HTSEQ_FILES[0]}" > "$HTSEQ_OUTPUT_DIR/gene_ids.txt"

# Create header dynamically from sample names
echo "  - Creating header..."
HEADER="Gene"
for htseq_file in "${HTSEQ_FILES[@]}"; do
    SAMPLE=$(basename "$htseq_file" .htseq.txt)
    HEADER="${HEADER}	${SAMPLE}"
done

# Extract counts from each file and merge
echo "  - Extracting sample counts..."
TEMP_DIR="$HTSEQ_OUTPUT_DIR/temp_merge"
mkdir -p "$TEMP_DIR"

for htseq_file in "${HTSEQ_FILES[@]}"; do
    SAMPLE=$(basename "$htseq_file" .htseq.txt)
    # FIXED: Extract counts (column 2) excluding __ summary lines
    awk '!/^__/ {print $2}' "$htseq_file" > "$TEMP_DIR/${SAMPLE}_counts.txt"
done

# Merge all columns
echo "  - Merging columns..."
MERGED_FILE="$HTSEQ_OUTPUT_DIR/HTSeq_Count_union.txt"

# Write header
echo -e "$HEADER" > "$MERGED_FILE"

# Paste gene IDs with all sample counts
paste "$HTSEQ_OUTPUT_DIR/gene_ids.txt" "$TEMP_DIR"/*_counts.txt >> "$MERGED_FILE"

# Cleanup temporary files
rm -rf "$TEMP_DIR" "$HTSEQ_OUTPUT_DIR/gene_ids.txt"

HT_LINES=$(wc -l < "$MERGED_FILE")
echo "✓ HTSeq quantification completed ($HT_LINES lines)"

# ===================================================================
# FINAL VALIDATION
# ===================================================================
echo ""
echo "=========================================="
echo "HISAT2 PIPELINE COMPLETE"
echo "=========================================="

echo ""
echo "Output Summary:"
echo "  ✓ Mapping directory: $MAPPING_DIR"
echo "    - BAM files: $BAM_COUNT"
echo "  ✓ featureCounts output: $FC_OUTPUT"
echo "  ✓ HTSeq output: $MERGED_FILE"
echo ""
echo "⚠️  Note: HISAT2 does NOT support RSEM quantification"
echo "   (Use STAR or Bowtie2 for RSEM-based analysis)"
echo ""
echo "✅ HISAT2 pipeline completed successfully!"