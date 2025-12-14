#!/bin/bash
# ===================================================================
# STAR COMPLETE PIPELINE: ALIGNMENT + QUANTIFICATION + RSEM
# ===================================================================
# This script performs:
# 1. STAR genome indexing (if needed)
# 2. RSEM reference preparation (if needed)
# 3. STAR alignment (PE/SE mode)
# 4. SAM to BAM conversion and sorting
# 5. featureCounts quantification
# 6. HTSeq quantification (with optimized piping)
# 7. RSEM quantification
# ===================================================================

set -e  # Exit on error
set -o pipefail  # Exit on pipe failure

# === CONFIGURATION ===
ALIGNER="STAR"
INDEX_DIR_STAR="Indices/STAR/genome"
RSEM_INDEX_DIR="Indices/RSEM/STAR"
MAPPING_DIR="Mapping/STAR"
FC_OUTPUT_DIR="Quantification/STAR/FC"
HTSEQ_OUTPUT_DIR="Quantification/STAR/HT"
RSEM_OUTPUT_DIR="Quantification/STAR/RSEM"

# === CREATE DIRECTORIES ===
mkdir -p "$MAPPING_DIR" "$FC_OUTPUT_DIR" "$HTSEQ_OUTPUT_DIR" "$RSEM_OUTPUT_DIR"

echo ""
echo "=========================================="
echo "STAR COMPLETE PIPELINE (with RSEM)"
echo "=========================================="
echo "Mode: $MODE"
echo "Threads: $THREADS"
echo "Genome: $FASTA"
echo "GTF: $GTF"
echo ""

# ===================================================================
# STEP 1: STAR INDEX PREPARATION
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 1: STAR Index Preparation"
echo "=========================================="

# Check if we're using pre-built indices
if [ -n "$INDEX_DIR" ] && [ -d "$INDEX_DIR/STAR" ]; then
    echo ">>> Using pre-built STAR index from: $INDEX_DIR/STAR"
    INDEX_DIR_STAR="$INDEX_DIR/STAR/genome"
    if [ ! -f "$INDEX_DIR_STAR/SA" ]; then
        echo "ERROR: STAR index not found at $INDEX_DIR_STAR" >&2
        exit 1
    fi
else
    echo ">>> Building STAR index locally..."
    mkdir -p "$INDEX_DIR_STAR"
    
    if [ ! -f "$INDEX_DIR_STAR/SA" ]; then
        echo "  Creating STAR index..."
        echo "  This may take a while..."
        
        STAR --runThreadN "$THREADS" \
            --runMode genomeGenerate \
            --genomeDir "$INDEX_DIR_STAR" \
            --genomeFastaFiles "$FASTA" \
            --sjdbGTFfile "$GTF" \
            --sjdbOverhang 100
        
        echo "  ✓ STAR index created successfully"
    else
        echo "  ✓ STAR index already exists"
    fi
fi

# ===================================================================
# STEP 1.5: RSEM REFERENCE PREPARATION
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 1.5: RSEM Reference Preparation"
echo "=========================================="

if [ ! -f "$RSEM_INDEX_DIR/rsem_ref.idx.fa" ]; then
    echo ">>> Building RSEM reference (STAR-based)..."
    mkdir -p "$RSEM_INDEX_DIR"
    
    rsem-prepare-reference \
        --gtf "$GTF" \
        "$FASTA" \
        "$RSEM_INDEX_DIR/rsem_ref" \
        --star \
        -p "$THREADS"
    
    echo "  ✓ RSEM reference created successfully"
else
    echo "  ✓ RSEM reference already exists"
fi

# ===================================================================
# STEP 2: STAR ALIGNMENT
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 2: STAR Alignment"
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
        
        SAMPLE_DIR="$MAPPING_DIR/${SAMPLE}"
        mkdir -p "$SAMPLE_DIR"
        
        echo "    - Running STAR alignment..."
        STAR --runThreadN "$THREADS" \
            --genomeDir "$INDEX_DIR_STAR" \
            --readFilesIn "$R1" "$R2" \
            --outFileNamePrefix "$SAMPLE_DIR/" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            2>&1 | grep -E "(input|output|Completed)" || true
        
        # Verify output
        SAM_FILE="$SAMPLE_DIR/Aligned.sortedByCoord.out.bam"
        if [ ! -f "$SAM_FILE" ]; then
            echo "ERROR: STAR alignment failed for $SAMPLE" >&2
            exit 1
        fi
        
        # Rename to standardized name
        SORTED_BAM="$MAPPING_DIR/${SAMPLE}_sorted.bam"
        mv "$SAM_FILE" "$SORTED_BAM"
        
        # Index BAM
        echo "    - Indexing BAM file..."
        samtools index -@ "$THREADS" "$SORTED_BAM"
        
        # Cleanup STAR temp directory
        rm -rf "$SAMPLE_DIR"
        
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
        
        SAMPLE_DIR="$MAPPING_DIR/${SAMPLE}"
        mkdir -p "$SAMPLE_DIR"
        
        echo "    - Running STAR alignment..."
        STAR --runThreadN "$THREADS" \
            --genomeDir "$INDEX_DIR_STAR" \
            --readFilesIn "$R" \
            --outFileNamePrefix "$SAMPLE_DIR/" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            2>&1 | grep -E "(input|output|Completed)" || true
        
        # Verify output
        SAM_FILE="$SAMPLE_DIR/Aligned.sortedByCoord.out.bam"
        if [ ! -f "$SAM_FILE" ]; then
            echo "ERROR: STAR alignment failed for $SAMPLE" >&2
            exit 1
        fi
        
        # Rename to standardized name
        SORTED_BAM="$MAPPING_DIR/${SAMPLE}_sorted.bam"
        mv "$SAM_FILE" "$SORTED_BAM"
        
        # Index BAM
        echo "    - Indexing BAM file..."
        samtools index -@ "$THREADS" "$SORTED_BAM"
        
        # Cleanup STAR temp directory
        rm -rf "$SAMPLE_DIR"
        
        echo "    ✓ $SAMPLE completed"
    done
fi

echo ""
echo "✓ STAR alignment completed"

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
# STEP 5: RSEM QUANTIFICATION
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 5: RSEM Quantification"
echo "=========================================="

echo ">>> Running RSEM quantification..."

if [[ $MODE == "PE" ]]; then
    echo "  Mode: Paired-End"
    
    for R1 in "$TRIM_DIR"/*_R1.fastq; do
        if [ ! -f "$R1" ]; then
            continue
        fi
        
        SAMPLE=$(basename "$R1" _R1.fastq)
        R2="$TRIM_DIR/${SAMPLE}_R2.fastq"
        
        if [ ! -f "$R2" ]; then
            continue
        fi
        
        echo ""
        echo "  Quantifying sample: $SAMPLE"
        
        rsem-calculate-expression \
            --paired-end \
            --star \
            --estimate-rspd \
            --append-names \
            --num-threads "$THREADS" \
            "$R1" "$R2" \
            "$RSEM_INDEX_DIR/rsem_ref" \
            "$RSEM_OUTPUT_DIR/$SAMPLE" \
            2> "$RSEM_OUTPUT_DIR/${SAMPLE}.rsem.log"
        
        if [ ! -f "$RSEM_OUTPUT_DIR/${SAMPLE}.genes.results" ]; then
            echo "    ERROR: RSEM failed for $SAMPLE" >&2
            exit 1
        fi
        
        echo "    ✓ RSEM completed for $SAMPLE"
    done

else
    echo "  Mode: Single-End"
    
    for R in "$TRIM_DIR"/*.fastq; do
        if [ ! -f "$R" ]; then
            continue
        fi
        
        # Skip paired-end files in SE mode
        if [[ "$R" == *_R1.fastq || "$R" == *_R2.fastq ]]; then
            continue
        fi
        
        SAMPLE=$(basename "$R" .fastq)
        
        echo ""
        echo "  Quantifying sample: $SAMPLE"
        
        rsem-calculate-expression \
            --star \
            --estimate-rspd \
            --append-names \
            --num-threads "$THREADS" \
            "$R" \
            "$RSEM_INDEX_DIR/rsem_ref" \
            "$RSEM_OUTPUT_DIR/$SAMPLE" \
            2> "$RSEM_OUTPUT_DIR/${SAMPLE}.rsem.log"
        
        if [ ! -f "$RSEM_OUTPUT_DIR/${SAMPLE}.genes.results" ]; then
            echo "    ERROR: RSEM failed for $SAMPLE" >&2
            exit 1
        fi
        
        echo "    ✓ RSEM completed for $SAMPLE"
    done
fi

echo ""
echo ">>> Merging RSEM results into count matrix..."

# Get gene_id column from first file
FIRST_FILE=$(ls "$RSEM_OUTPUT_DIR"/*.genes.results | head -n1)
cut -f1 "$FIRST_FILE" | tail -n +2 > "$RSEM_OUTPUT_DIR/gene_ids.txt"

# Create temporary directory for counts
TEMP_DIR="$RSEM_OUTPUT_DIR/temp_counts"
mkdir -p "$TEMP_DIR"

# Extract raw counts (5th column) from each file
for FILE in "$RSEM_OUTPUT_DIR"/*.genes.results; do
    SAMPLE=$(basename "$FILE" .genes.results)
    echo "  - Processing $SAMPLE..."
    cut -f5 "$FILE" | tail -n +2 > "$TEMP_DIR/${SAMPLE}.txt"
done

# Merge gene_ids with all sample count files
paste "$RSEM_OUTPUT_DIR/gene_ids.txt" "$TEMP_DIR"/*.txt > "$RSEM_OUTPUT_DIR/rsem_gene_counts_matrix.txt"

# Prepend header to the final matrix
HEADER="Geneid"
for FILE in "$TEMP_DIR"/*.txt; do
    SAMPLE=$(basename "$FILE" .txt)
    HEADER="$HEADER	$SAMPLE"
done

(echo -e "$HEADER"; cat "$RSEM_OUTPUT_DIR/rsem_gene_counts_matrix.txt") > "$RSEM_OUTPUT_DIR/tmpfile" && \
mv "$RSEM_OUTPUT_DIR/tmpfile" "$RSEM_OUTPUT_DIR/rsem_gene_counts_matrix.txt"

# Cleanup temporary files
rm -rf "$TEMP_DIR" "$RSEM_OUTPUT_DIR/gene_ids.txt"

RSEM_LINES=$(wc -l < "$RSEM_OUTPUT_DIR/rsem_gene_counts_matrix.txt")
echo "✓ RSEM quantification completed ($RSEM_LINES lines)"

# ===================================================================
# FINAL VALIDATION
# ===================================================================
echo ""
echo "=========================================="
echo "STAR PIPELINE COMPLETE"
echo "=========================================="

echo ""
echo "Output Summary:"
echo "  ✓ Mapping directory: $MAPPING_DIR"
echo "    - BAM files: $BAM_COUNT"
echo "  ✓ featureCounts output: $FC_OUTPUT"
echo "  ✓ HTSeq output: $MERGED_FILE"
echo "  ✓ RSEM output: $RSEM_OUTPUT_DIR/rsem_gene_counts_matrix.txt"
echo ""
echo "✓ STAR pipeline completed successfully!"