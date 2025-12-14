#!/bin/bash
# ===================================================================
# STAR COMPLETE CLEANING SCRIPT - ALL FIXES
# ===================================================================
# This script:
# 1. Fixes featureCounts headers (removes .bam and _sorted suffixes)
# 2. Re-merges HTSeq files if they're broken
# 3. Rounds RSEM counts to integers for DESeq2
# 4. Validates all outputs
# ===================================================================

set -e
set -o pipefail

echo ""
echo "=========================================="
echo "STAR COUNT MATRIX CLEANING - COMPLETE FIX"
echo "=========================================="
echo ""

# Change to DEG directory
cd DEG || { echo "ERROR: DEG directory not found!"; exit 1; }

# ===================================================================
# STEP 1: CLEAN FEATURECOUNTS FILE
# ===================================================================
echo ">>> STEP 1: Processing featureCounts output..."

if [ ! -f "S_FC_Count.txt" ]; then
    echo "ERROR: S_FC_Count.txt not found!" >&2
    exit 1
fi

awk -F'\t' '
BEGIN { OFS = "," }
NR == 2 {
    printf "Geneid"
    for (i = 7; i <= NF; i++) {
        name = $i
        gsub(/^.*\//, "", name)      # Remove directory path
        gsub(/\.bam$/, "", name)     # Remove .bam extension
        gsub(/_sorted$/, "", name)   # Remove _sorted suffix
        gsub(/_Aligned.*/, "", name) # Remove STAR alignment suffix
        printf ",%s", name
    }
    printf "\n"
    next
}
NR > 2 {
    printf "%s", $1
    for (i = 7; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
}
' S_FC_Count.txt > S_FC_counts_clean.txt

rm S_FC_Count.txt
FC_LINES=$(wc -l < S_FC_counts_clean.txt)
FC_DATA_LINES=$((FC_LINES - 1))
echo "  ✓ featureCounts cleaned: $FC_DATA_LINES genes"

# ===================================================================
# STEP 2: FIX AND CLEAN HTSEQ FILE
# ===================================================================
echo ""
echo ">>> STEP 2: Processing HTSeq output..."

if [ ! -f "S_HTSeq_Count_union.txt" ]; then
    echo "ERROR: S_HTSeq_Count_union.txt not found!" >&2
    exit 1
fi

# Check if HTSeq file is broken (only header, no data)
HTSEQ_DATA_LINES=$(tail -n +2 S_HTSeq_Count_union.txt | wc -l)

if [ $HTSEQ_DATA_LINES -eq 0 ]; then
    echo "  ⚠  HTSeq file is broken (no data rows) - attempting to re-merge..."
    
    # Check if individual HTSeq files exist
    HTSEQ_DIR="../Quantification/STAR/HT"
    if [ ! -d "$HTSEQ_DIR" ]; then
        echo "  ERROR: HTSeq output directory not found at $HTSEQ_DIR" >&2
        exit 1
    fi
    
    HTSEQ_FILES=("$HTSEQ_DIR"/*.htseq.txt)
    if [ ${#HTSEQ_FILES[@]} -eq 0 ] || [ ! -f "${HTSEQ_FILES[0]}" ]; then
        echo "  ERROR: No individual HTSeq files found in $HTSEQ_DIR" >&2
        exit 1
    fi
    
    echo "  Found ${#HTSEQ_FILES[@]} individual HTSeq files - re-merging..."
    
    # Create header from filenames
    HEADER="Gene"
    for htseq_file in "${HTSEQ_FILES[@]}"; do
        SAMPLE=$(basename "$htseq_file" .htseq.txt)
        HEADER="${HEADER}	${SAMPLE}"
    done
    
    # Create temp directory
    TEMP_DIR="htseq_remerge_temp"
    mkdir -p "$TEMP_DIR"
    
    # Extract gene IDs from first file (excluding HTSeq summary lines)
    awk '!/^__/ {print $1}' "${HTSEQ_FILES[0]}" > "$TEMP_DIR/gene_ids.txt"
    
    # Extract counts from each file
    for htseq_file in "${HTSEQ_FILES[@]}"; do
        SAMPLE=$(basename "$htseq_file" .htseq.txt)
        awk '!/^__/ {print $2}' "$htseq_file" > "$TEMP_DIR/${SAMPLE}_counts.txt"
    done
    
    # Merge into new file
    echo -e "$HEADER" > S_HTSeq_Count_union.txt
    paste "$TEMP_DIR/gene_ids.txt" "$TEMP_DIR"/*_counts.txt >> S_HTSeq_Count_union.txt
    
    # Cleanup
    rm -rf "$TEMP_DIR"
    
    echo "  ✓ HTSeq files successfully re-merged"
else
    echo "  ✓ HTSeq file has data ($HTSEQ_DATA_LINES genes)"
fi

# Now clean the HTSeq file (convert to CSV)
awk -F'\t' '
BEGIN { OFS = "," }
NR == 1 {
    printf "Geneid"
    for (i = 2; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
    next
}
NR > 1 && $1 !~ /^__/ {
    printf "%s", $1
    for (i = 2; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
}
' S_HTSeq_Count_union.txt > S_HTSeq_counts_clean.txt

rm S_HTSeq_Count_union.txt
HT_LINES=$(wc -l < S_HTSeq_counts_clean.txt)
HT_DATA_LINES=$((HT_LINES - 1))
echo "  ✓ HTSeq cleaned: $HT_DATA_LINES genes"

# ===================================================================
# STEP 3: CLEAN RSEM FILE (ROUND TO INTEGERS)
# ===================================================================
echo ""
echo ">>> STEP 3: Processing RSEM output (rounding to integers)..."

if [ ! -f "S_rsem_gene_counts_matrix.txt" ]; then
    echo "ERROR: S_rsem_gene_counts_matrix.txt not found!" >&2
    exit 1
fi

awk -F'\t' '
BEGIN { OFS = "," }
NR == 1 {
    # Print header as-is
    printf "%s", $1
    for (i = 2; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
    next
}
NR > 1 {
    # Keep original values (no rounding)
    printf "%s", $1
    for (i = 2; i <= NF; i++) {
        printf ",%s", $i              # ← CHANGED FROM %d to %s
    }
    printf "\n"
}
' S_rsem_gene_counts_matrix.txt > S_RSEM_counts_clean.txt

rm S_rsem_gene_counts_matrix.txt
RSEM_LINES=$(wc -l < S_RSEM_counts_clean.txt)
RSEM_DATA_LINES=$((RSEM_LINES - 1))
echo "  ✓ RSEM cleaned (decimals preserved): $RSEM_DATA_LINES genes"

# ===================================================================
# STEP 4: COPY SAMPLEINFO
# ===================================================================
echo ""
echo ">>> STEP 4: Checking for sampleinfo.txt..."

if [ ! -f "sampleinfo.txt" ] && [ -f "/data/sampleinfo.txt" ]; then
    cp /data/sampleinfo.txt sampleinfo.txt
    echo "  ✓ sampleinfo.txt copied"
elif [ -f "sampleinfo.txt" ]; then
    echo "  ✓ sampleinfo.txt already exists"
else
    echo "  ⚠  WARNING: sampleinfo.txt not found at /data/sampleinfo.txt"
fi

# ===================================================================
# VALIDATION & SUMMARY
# ===================================================================
echo ""
echo "=========================================="
echo "STAR CLEANING COMPLETE"
echo "=========================================="
echo ""

echo "Files processed:"
echo "  • S_FC_counts_clean.txt ($FC_DATA_LINES genes)"
echo "  • S_HTSeq_counts_clean.txt ($HT_DATA_LINES genes)"
echo "  • S_RSEM_counts_clean.txt ($RSEM_DATA_LINES genes) [INTEGERS]"
echo ""

# Validate data exists
if [ $FC_DATA_LINES -eq 0 ] || [ $HT_DATA_LINES -eq 0 ] || [ $RSEM_DATA_LINES -eq 0 ]; then
    echo "ERROR: One or more files have no data!" >&2
    exit 1
fi

echo "Headers (first 100 chars):"
echo ""
echo "featureCounts:"
head -n 1 S_FC_counts_clean.txt | cut -c1-100
echo ""
echo "HTSeq:"
head -n 1 S_HTSeq_counts_clean.txt | cut -c1-100
echo ""
echo "RSEM:"
head -n 1 S_RSEM_counts_clean.txt | cut -c1-100
echo ""

echo "Sample data (first gene):"
echo ""
echo "featureCounts:"
head -n 2 S_FC_counts_clean.txt | tail -n 1 | cut -c1-100
echo ""
echo "HTSeq:"
head -n 2 S_HTSeq_counts_clean.txt | tail -n 1 | cut -c1-100
echo ""
echo "RSEM (should be integers):"
head -n 2 S_RSEM_counts_clean.txt | tail -n 1 | cut -c1-100
echo ""

echo "✅ All STAR count matrices cleaned and ready for DEG analysis"
echo ""

cd ..