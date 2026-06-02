#!/bin/bash
# ===================================================================
# master_busco.sh — BUSCO Pipeline Orchestrator
# ===================================================================
# Called at runtime (NOT during Docker build).
# Runs inside the container using two isolated environments:
#
#   rnaseq-busco    → BUSCO 6.0.0 + SEPP  (Python 3.9, pandas 1.5.x)
#   rnaseq-busco-py → gffpandas filter     (Python 3.11, pandas 2.x)
#
# Required environment variables (set via -e in docker run):
#   PROTEIN_FILE   path to protein.faa inside container
#   BUSCO_LINEAGE  e.g. mammalia_odb12
#   GFF_FILE       path to genomic.gff inside container
#   GTF_FILE       path to genomic.gtf inside container
#   THREADS        number of CPU threads (default: 8)
#
# Output written to: /data/busco_output/
#
# Example docker run command:
#   docker run --rm \
#     -v $(pwd)/results:/data \
#     -v $(pwd)/data/genome:/data/genome:ro \
#     -e PROTEIN_FILE=/data/genome/protein.faa \
#     -e BUSCO_LINEAGE=mammalia_odb12 \
#     -e GFF_FILE=/data/genome/genomic.gff \
#     -e GTF_FILE=/data/genome/genomic.gtf \
#     -e THREADS=16 \
#     rnaseqmetaanalyst:latest \
#     /bin/bash -c "/opt/project/scripts/master_busco.sh"
# ===================================================================

set -o pipefail

echo ""
echo "=========================================="
echo "BUSCO PIPELINE STARTED: $(date)"
echo "=========================================="

# ===================================================================
# VALIDATE ENVIRONMENT VARIABLES
# ===================================================================
echo ""
echo ">>> Validating inputs..."

THREADS="${THREADS:-8}"
MISSING=0

for var in PROTEIN_FILE BUSCO_LINEAGE GFF_FILE GTF_FILE; do
    if [ -z "${!var}" ]; then
        echo "  ERROR: $var is not set" >&2
        MISSING=1
    fi
done

[ $MISSING -ne 0 ] && exit 1

for file in "$PROTEIN_FILE" "$GFF_FILE" "$GTF_FILE"; do
    if [ ! -f "$file" ]; then
        echo "  ERROR: File not found: $file" >&2
        exit 1
    fi
    echo "  [OK] $file"
done

echo "  [OK] Lineage : $BUSCO_LINEAGE"
echo "  [OK] Threads : $THREADS"

# ===================================================================
# OUTPUT DIRECTORIES
# ===================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUT_DIR="/data/busco_output"
BUSCO_RUN_DIR="$OUT_DIR/busco_run"
FILTER_DIR="$OUT_DIR/filtered"
BUSCO_LIST="$OUT_DIR/busco_core_genes1.txt"
LOG_DIR="/data/logs"

mkdir -p "$OUT_DIR" "$FILTER_DIR" "$LOG_DIR"

# ===================================================================
# STAGE 1: RUN BUSCO
# Environment: rnaseq-busco (Python 3.9, BUSCO 6.0.0, SEPP 4.5.5)
# ===================================================================
echo ""
echo "=========================================="
echo "STAGE 1: BUSCO Assessment"
echo "Environment: rnaseq-busco"
echo "=========================================="

BUSCO_LOG="$LOG_DIR/busco_${TIMESTAMP}.log"
echo "  Log: $BUSCO_LOG"
echo ""

if micromamba run -n rnaseq-busco \
    busco \
        -i "$PROTEIN_FILE" \
        -l "$BUSCO_LINEAGE" \
        -o busco_run \
        --out_path "$OUT_DIR" \
        -m proteins \
        -f \
        -c "$THREADS" \
    > "$BUSCO_LOG" 2>&1; then
    echo "  ✓ BUSCO completed successfully"
else
    echo "  ✗ BUSCO failed — check log: $BUSCO_LOG" >&2
    tail -20 "$BUSCO_LOG" >&2
    exit 1
fi

# ===================================================================
# STAGE 2: EXTRACT COMPLETE SINGLE-COPY BUSCO GENE IDs
# Uses awk — no environment needed
# ===================================================================
echo ""
echo "=========================================="
echo "STAGE 2: Extracting Complete BUSCO Gene IDs"
echo "=========================================="

FULL_TABLE="$BUSCO_RUN_DIR/run_${BUSCO_LINEAGE}/full_table.tsv"

if [ ! -f "$FULL_TABLE" ]; then
    echo "  ERROR: full_table.tsv not found at: $FULL_TABLE" >&2
    echo "  Check BUSCO log: $BUSCO_LOG" >&2
    exit 1
fi

awk -F'\t' '$2=="Complete" {
    printf $3 "\t"
    for (i=7; i<=NF; i++) { printf $i " " }
    print ""
}' "$FULL_TABLE" > "$BUSCO_LIST"

BUSCO_COUNT=$(wc -l < "$BUSCO_LIST")
echo "  ✓ Complete BUSCO genes extracted: $BUSCO_COUNT"
echo "  ✓ Saved: $BUSCO_LIST"

if [ "$BUSCO_COUNT" -eq 0 ]; then
    echo "  ERROR: No complete BUSCO genes found in full_table.tsv" >&2
    exit 1
fi

# ===================================================================
# STAGE 3: FILTER GFF3 AND GTF TO BUSCO GENES
# Environment: rnaseq-busco-py (Python 3.11, gffpandas, pandas 2.x)
# ===================================================================
echo ""
echo "=========================================="
echo "STAGE 3: GFF3 / GTF Filtering"
echo "Environment: rnaseq-busco-py"
echo "=========================================="

FILTER_LOG="$LOG_DIR/busco_filter_${TIMESTAMP}.log"
echo "  Log: $FILTER_LOG"
echo ""

if micromamba run -n rnaseq-busco-py \
    python /opt/project/scripts/busco_filter.py \
        --gff     "$GFF_FILE" \
        --gtf     "$GTF_FILE" \
        --busco   "$BUSCO_LIST" \
        --out-dir "$FILTER_DIR" \
    > "$FILTER_LOG" 2>&1; then
    echo "  ✓ GFF3/GTF filtering completed successfully"
    cat "$FILTER_LOG"
else
    echo "  ✗ Filtering failed — check log: $FILTER_LOG" >&2
    tail -20 "$FILTER_LOG" >&2
    exit 1
fi

# ===================================================================
# FINAL SUMMARY
# ===================================================================
echo ""
echo "=========================================="
echo "BUSCO PIPELINE COMPLETED: $(date)"
echo "=========================================="
echo ""
echo "Output Summary:"
echo "──────────────────────────────────────────"
echo ""
echo "Stage 1 — BUSCO Assessment:"
echo "  ✓ Full results   : $BUSCO_RUN_DIR/"
echo "  ✓ Full table     : $FULL_TABLE"
echo ""
echo "Stage 2 — Core Gene List:"
echo "  ✓ Complete genes : $BUSCO_COUNT"
echo "  ✓ Gene list file : $BUSCO_LIST"
echo ""
echo "Stage 3 — Filtered Annotations:"
echo "  ✓ Filtered GFF3  : $FILTER_DIR/busco_filtered.gff"
echo "  ✓ Filtered GTF   : $FILTER_DIR/busco_filtered.gtf"
echo "  ✓ Match table    : $FILTER_DIR/busco_matched.tsv"
echo ""
echo "Log Files:"
echo "  → BUSCO log      : $BUSCO_LOG"
echo "  → Filter log     : $FILTER_LOG"
echo ""
echo "──────────────────────────────────────────"
echo "✅ ALL BUSCO STAGES COMPLETED"
echo "──────────────────────────────────────────"
echo ""
