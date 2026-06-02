#!/bin/bash
# ===================================================================
# RUN ALL DESeq2 ANALYSES
# ===================================================================
# Loops over every *_counts_clean.txt in DEG/ and runs DESeq2.
# Housekeeping genes (housekeeping_genes.xlsx) are located
# automatically from several candidate paths and copied to /data/
# so that run_deseq2_analysis.R (--project-dir /data) always finds
# them.  The BUSCO pipeline writes them to /data/ by default;
# this script handles the case where they live elsewhere.
# ===================================================================

set -e

# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------
RSCRIPT="/opt/project/scripts/run_deseq2_analysis.R"
PROJECT_DIR="/data"
DEG_DIR="DEG"

# -------------------------------------------------------------------
# VALIDATE SAMPLEINFO
# -------------------------------------------------------------------
# Accept SAMPLE_INFO_PATH env-var; fall back to DEG/sampleinfo.txt
if [ -z "$SAMPLE_INFO_PATH" ]; then
    SAMPLE_INFO_PATH="${DEG_DIR}/sampleinfo.txt"
fi

if [ ! -f "$SAMPLE_INFO_PATH" ]; then
    echo "ERROR: sampleinfo.txt not found at '$SAMPLE_INFO_PATH'" >&2
    echo "  Set the SAMPLE_INFO_PATH environment variable or place" >&2
    echo "  sampleinfo.txt inside the DEG/ directory." >&2
    exit 1
fi
echo "  ✓ Sample info: $SAMPLE_INFO_PATH"

# -------------------------------------------------------------------
# LOCATE AND COPY HOUSEKEEPING GENES FILE
# -------------------------------------------------------------------
echo ""
echo ">>> Locating housekeeping_genes.xlsx ..."

HKG_DEST="${PROJECT_DIR}/housekeeping_genes.xlsx"

# Candidate source locations (in order of preference)
HKG_CANDIDATES=(
    "${PROJECT_DIR}/housekeeping_genes.xlsx"                    # already in /data
    "${PROJECT_DIR}/busco_output/filtered/housekeeping_genes.xlsx"  # BUSCO filter output
    "/data/busco_output/filtered/housekeeping_genes.xlsx"
    "${DEG_DIR}/housekeeping_genes.xlsx"                        # manually placed in DEG/
    "./housekeeping_genes.xlsx"                                  # working directory
)

HKG_FOUND=""
for candidate in "${HKG_CANDIDATES[@]}"; do
    if [ -f "$candidate" ]; then
        HKG_FOUND="$candidate"
        break
    fi
done

if [ -n "$HKG_FOUND" ]; then
    if [ "$HKG_FOUND" != "$HKG_DEST" ]; then
        echo "  Found: $HKG_FOUND"
        cp "$HKG_FOUND" "$HKG_DEST"
        echo "  ✓ Copied to: $HKG_DEST"
    else
        echo "  ✓ Already at: $HKG_DEST"
    fi
else
    echo "  ⚠  WARNING: housekeeping_genes.xlsx not found in any candidate path."
    echo "     HKG analysis will be skipped by the R script."
    echo "     Run the BUSCO pipeline first, or place the file at:"
    echo "       $HKG_DEST"
    echo "     Candidate paths searched:"
    for candidate in "${HKG_CANDIDATES[@]}"; do
        echo "       $candidate"
    done
fi

# -------------------------------------------------------------------
# RUN DESeq2 FOR EACH COUNT MATRIX
# -------------------------------------------------------------------
echo ""
echo "=========================================="
echo "DESeq2 ANALYSIS"
echo "=========================================="

COUNT_FILES=( "${DEG_DIR}"/*_counts_clean.txt )
if [ ${#COUNT_FILES[@]} -eq 0 ] || [ ! -f "${COUNT_FILES[0]}" ]; then
    echo "ERROR: No *_counts_clean.txt files found in ${DEG_DIR}/" >&2
    exit 1
fi

echo "  Count files found: ${#COUNT_FILES[@]}"
echo ""

for count_file in "${COUNT_FILES[@]}"; do
    echo "-----------------------------"
    echo "Running DESeq2 for: $count_file"
    Rscript "$RSCRIPT" \
        "$count_file" \
        "$SAMPLE_INFO_PATH" \
        --project-dir "$PROJECT_DIR"
    echo "  ✓ Finished: $count_file"
done

echo ""
echo "=========================================="
echo "✅ All DESeq2 analyses complete!"
echo "   Results: ${PROJECT_DIR}/results/DEG/"
echo "=========================================="
