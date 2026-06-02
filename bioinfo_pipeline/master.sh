#!/bin/bash
# ===================================================================
# RNA-Seq Pipeline Master Script - FIXED VERSION
# ===================================================================
# FIXES:
# 1. BOWTIE_FAILED / HISAT2_FAILED / STAR_FAILED were set inside
#    subshells (parallel function) and never returned to parent shell.
#    All file-copy and clean steps now use direct file-existence checks
#    instead of those broken variables.
# 2. Sequential mode had same scoping issue - fixed same way.
# 3. set -e removed from top level so one aligner failure does not
#    abort the entire pipeline; each stage handles its own errors.
# 4. Minimum-aligner guard raised to 1 (was silently 2 in parallel).
# ===================================================================

set -o pipefail

echo "=========================================="
echo "RNA-Seq Pipeline Started: $(date)"
echo "=========================================="

# ===================================================================
# PARSE COMMAND LINE ARGUMENTS
# ===================================================================
RUN_ALIGNMENT="true"
RUN_DEG="true"
RUN_COMPARISON="true"
ALIGNMENT_MODE="parallel"

export LOG2FC_THRESHOLD="${LOG2FC_THRESHOLD:-0}"
export FDR_THRESHOLD="${FDR_THRESHOLD:-0.05}"
export PVALUE_THRESHOLD="${PVALUE_THRESHOLD:-0.05}"
export STAR_SJDB_OVERHANG="${STAR_SJDB_OVERHANG:-100}"
export INDEX_DIR="${INDEX_DIR:-}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-alignment)
            RUN_ALIGNMENT="false"
            shift ;;
        --skip-deg)
            RUN_DEG="false"
            RUN_COMPARISON="false"
            shift ;;
        --skip-comparison)
            RUN_COMPARISON="false"
            shift ;;
        --deg-only)
            RUN_ALIGNMENT="false"
            RUN_DEG="true"
            shift ;;
        --comparison-only)
            RUN_ALIGNMENT="false"
            RUN_DEG="false"
            RUN_COMPARISON="true"
            shift ;;
        --deg-and-comparison)
            RUN_ALIGNMENT="false"
            RUN_DEG="true"
            RUN_COMPARISON="true"
            shift ;;
        --sequential)
            ALIGNMENT_MODE="sequential"
            shift ;;
        --sequential-only)
            ALIGNMENT_MODE="sequential_only"
            shift ;;
        --log2fc)
            export LOG2FC_THRESHOLD="$2"
            shift 2 ;;
        --fdr)
            export FDR_THRESHOLD="$2"
            shift 2 ;;
        --pvalue)
            export PVALUE_THRESHOLD="$2"
            shift 2 ;;
        --sjdb-overhang)
            export STAR_SJDB_OVERHANG="$2"
            shift 2 ;;
        --index-dir)
            export INDEX_DIR="$2"
            shift 2 ;;
        --help)
            cat << 'EOF'
Usage: master.sh [OPTIONS]

Pipeline Control:
  --skip-alignment     Skip alignment & quantification
  --skip-deg           Skip DEG analysis
  --skip-comparison    Skip comparison analysis
  --deg-only           Only run DEG (requires existing count matrices)
  --comparison-only    Only run comparison (requires existing DEG results)
  --deg-and-comparison Run DEG and comparison (skip alignment)

Alignment Execution Mode:
  --sequential         Run aligners sequentially (one at a time)
  --sequential-only    Alias for --sequential
  Default: parallel    Run all 3 aligners in parallel

Reference Indices:
  --index-dir <path>   Path to pre-built indices directory

DEG Thresholds:
  --log2fc <value>     log2 Fold Change threshold (default: 0)
  --fdr <value>        FDR threshold (default: 0.05)
  --pvalue <value>     p-value threshold (default: 0.05)

Alignment Options:
  --sjdb-overhang <value>  STAR parameter (default: 100)

Environment Variables:
  INDEX_DIR    Alternative to --index-dir flag
  THREADS      Number of CPU threads
  MODE         SE or PE sequencing mode
  GENOME_DIR   Path to genome files
  READ_DIR     Path to FASTQ files
  TRIM_DIR     Path for trimmed reads
  GTF, FASTA, GFF  Reference file paths
EOF
            exit 0 ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            exit 1 ;;
    esac
done

# ===================================================================
# SETUP
# ===================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

mkdir -p /data/logs /data/DEG /data/results/DEG /data/results/pipeline_comparison

# Ensure sampleinfo.txt is present in /data/DEG/
if [ -f "/data/sampleinfo.txt" ]; then
    cp /data/sampleinfo.txt /data/DEG/sampleinfo.txt 2>/dev/null || true
    echo "✓ sampleinfo.txt copied to /data/DEG/"
fi

if [ ! -f "/data/DEG/sampleinfo.txt" ]; then
    echo "WARNING: sampleinfo.txt not found in /data/DEG/" >&2
fi

export SAMPLE_INFO_PATH="/data/DEG/sampleinfo.txt"

# ===================================================================
# HOUSEKEEPING GENES DISTRIBUTION
# ===================================================================
# Automatically finds housekeeping_genes.xlsx from any of these
# source locations (checked in priority order):
#
#   Priority 1: /data/busco_output/filtered/housekeeping_genes.xlsx
#               → auto-generated by busco_filter.py (rnaseq-busco-py)
#               → no user action needed if BUSCO pipeline was run first
#
#   Priority 2: /data/housekeeping_genes.xlsx
#               → user manually placed alongside sampleinfo.txt
#               → also written here by busco_filter.py as a second copy
#
#   Priority 3: /data/DEG/housekeeping_genes.xlsx
#               → already distributed from a previous pipeline run
#               → sync forward to /data/results/ only
#
# Distributed to:
#   /data/DEG/      → DESeq2 and edgeR R scripts read from here
#   /data/results/  → available parallel to DEG folder for reference
# ===================================================================
# ===================================================================
# SETUP
# ===================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
 
mkdir -p /data/logs /data/DEG /data/results/DEG /data/results/pipeline_comparison
 
# Ensure sampleinfo.txt is present in /data/DEG/
if [ -f "/data/sampleinfo.txt" ]; then
    cp /data/sampleinfo.txt /data/DEG/sampleinfo.txt 2>/dev/null || true
    echo "✓ sampleinfo.txt copied to /data/DEG/"
fi
 
if [ ! -f "/data/DEG/sampleinfo.txt" ]; then
    echo "WARNING: sampleinfo.txt not found in /data/DEG/" >&2
fi
 
export SAMPLE_INFO_PATH="/data/DEG/sampleinfo.txt"
 
# Ensure housekeeping_genes.xlsx is present in /data/DEG/
# Two independent if-blocks mirror the sampleinfo.txt pattern:
#   block 1 — BUSCO pipeline output (auto-generated, no user action)
#   block 2 — /data/ root           (user-placed, same dir as sampleinfo.txt)
# Both run independently so whichever source exists (or both) will copy.
# Last write wins; /data/ root takes precedence if both are present.
 
if [ -f "/data/busco_output/filtered/housekeeping_genes.xlsx" ]; then
    cp /data/busco_output/filtered/housekeeping_genes.xlsx \
       /data/DEG/housekeeping_genes.xlsx 2>/dev/null || true
    echo "✓ housekeeping_genes.xlsx copied to /data/DEG/ (source: busco_output)"
fi
 
if [ -f "/data/housekeeping_genes.xlsx" ]; then
    cp /data/housekeeping_genes.xlsx \
       /data/DEG/housekeeping_genes.xlsx 2>/dev/null || true
    echo "✓ housekeeping_genes.xlsx copied to /data/DEG/ (source: /data/)"
fi
 
if [ ! -f "/data/DEG/housekeeping_genes.xlsx" ]; then
    echo "WARNING: housekeeping_genes.xlsx not found in /data/DEG/" >&2
    echo "         HKG analysis will be skipped" >&2
    echo "         To enable: run BUSCO first, or place the file at" >&2
    echo "           /data/housekeeping_genes.xlsx  (next to sampleinfo.txt)" >&2
fi
 
export HKG_PATH="/data/DEG/housekeeping_genes.xlsx"
 
# ===================================================================
# VALIDATE ENVIRONMENT VARIABLES (only if running alignment)
# ===================================================================
if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    echo ""
    echo ">>> Validating environment variables..."
 
    MISSING_VARS=0
    for var in THREADS MODE GENOME_DIR READ_DIR TRIM_DIR GTF FASTA GFF; do
        if [ -z "${!var}" ]; then
            echo "ERROR: $var is not set" >&2
            MISSING_VARS=1
        fi
    done
    [ $MISSING_VARS -ne 0 ] && exit 1
 
    MODE=$(echo "$MODE" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')
    if [[ "$MODE" != "SE" && "$MODE" != "PE" ]]; then
        echo "ERROR: MODE must be 'SE' or 'PE', got: $MODE" >&2
        exit 1
    fi
 
    echo ">>> Checking for critical files..."
    for file in "$GTF" "$FASTA" "$GFF"; do
        if [ ! -f "$file" ]; then
            echo "ERROR: Required file not found: $file" >&2
            exit 1
        fi
        echo "  [OK] Found: $file"
    done
 
    if [ ! -d "$READ_DIR" ] || [ -z "$(ls -A "$READ_DIR")" ]; then
        echo "ERROR: READ_DIR is empty or doesn't exist: $READ_DIR" >&2
        exit 1
    fi
    echo "  [OK] Found reads in: $READ_DIR"
 
    echo ""
    echo ">>> Alignment execution mode: $ALIGNMENT_MODE"
fi
 


# ===================================================================
# VALIDATE ENVIRONMENT VARIABLES (only if running alignment)
# ===================================================================
if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    echo ""
    echo ">>> Validating environment variables..."

    MISSING_VARS=0
    for var in THREADS MODE GENOME_DIR READ_DIR TRIM_DIR GTF FASTA GFF; do
        if [ -z "${!var}" ]; then
            echo "ERROR: $var is not set" >&2
            MISSING_VARS=1
        fi
    done
    [ $MISSING_VARS -ne 0 ] && exit 1

    MODE=$(echo "$MODE" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')
    if [[ "$MODE" != "SE" && "$MODE" != "PE" ]]; then
        echo "ERROR: MODE must be 'SE' or 'PE', got: $MODE" >&2
        exit 1
    fi

    echo ">>> Checking for critical files..."
    for file in "$GTF" "$FASTA" "$GFF"; do
        if [ ! -f "$file" ]; then
            echo "ERROR: Required file not found: $file" >&2
            exit 1
        fi
        echo "  [OK] Found: $file"
    done

    if [ ! -d "$READ_DIR" ] || [ -z "$(ls -A "$READ_DIR")" ]; then
        echo "ERROR: READ_DIR is empty or doesn't exist: $READ_DIR" >&2
        exit 1
    fi
    echo "  [OK] Found reads in: $READ_DIR"

    echo ""
    echo ">>> Alignment execution mode: $ALIGNMENT_MODE"
fi

# ===================================================================
# HELPER FUNCTION: Run a single aligner (no subshell variable leaking)
# Returns 0 on success, 1 on failure.
# Usage: run_aligner <name> <script> <logfile>
# ===================================================================
run_aligner() {
    local NAME="$1"
    local SCRIPT="$2"
    local LOGFILE="$3"

    echo "    [$NAME] Starting... $(date)"
    if micromamba run -n rnaseq-cli bash -c "cd /data && $SCRIPT" > "$LOGFILE" 2>&1; then
        echo "    [$NAME] ✓ Completed successfully"
        return 0
    else
        echo "    [$NAME] ✗ Failed (see $LOGFILE)" >&2
        return 1
    fi
}

# ===================================================================
# HELPER FUNCTION: Copy count matrices for an aligner
# Uses file-existence checks — no dependency on *_FAILED variables.
# Usage: copy_matrices_bowtie / copy_matrices_hisat2 / copy_matrices_star
# ===================================================================
copy_matrices_bowtie() {
    local COPIED=0
    [ -f "/data/Quantification/Bowtie/FC/FC_Count.txt" ] && \
        cp "/data/Quantification/Bowtie/FC/FC_Count.txt" "/data/DEG/B_FC_Count.txt" && ((COPIED++))
    [ -f "/data/Quantification/Bowtie/HT/HTSeq_Count_union.txt" ] && \
        cp "/data/Quantification/Bowtie/HT/HTSeq_Count_union.txt" "/data/DEG/B_HTSeq_Count_union.txt" && ((COPIED++))
    [ -f "/data/Quantification/Bowtie/RSEM/rsem_gene_counts_matrix.txt" ] && \
        cp "/data/Quantification/Bowtie/RSEM/rsem_gene_counts_matrix.txt" "/data/DEG/B_rsem_gene_counts_matrix.txt" && ((COPIED++))
    if [ $COPIED -gt 0 ]; then
        echo "      ✓ Bowtie2: $COPIED matrix file(s) copied"
    else
        echo "      ⚠ Bowtie2: no quantification output found to copy"
    fi
}

copy_matrices_hisat2() {
    local COPIED=0
    [ -f "/data/Quantification/Hisat2/FC/FC_Count.txt" ] && \
        cp "/data/Quantification/Hisat2/FC/FC_Count.txt" "/data/DEG/H_FC_Count.txt" && ((COPIED++))
    [ -f "/data/Quantification/Hisat2/HT/HTSeq_Count_union.txt" ] && \
        cp "/data/Quantification/Hisat2/HT/HTSeq_Count_union.txt" "/data/DEG/H_HTSeq_Count_union.txt" && ((COPIED++))
    if [ $COPIED -gt 0 ]; then
        echo "      ✓ HISAT2: $COPIED matrix file(s) copied"
    else
        echo "      ⚠ HISAT2: no quantification output found to copy"
    fi
}

copy_matrices_star() {
    local COPIED=0
    [ -f "/data/Quantification/STAR/FC/FC_Count.txt" ] && \
        cp "/data/Quantification/STAR/FC/FC_Count.txt" "/data/DEG/S_FC_Count.txt" && ((COPIED++))
    [ -f "/data/Quantification/STAR/HT/HTSeq_Count_union.txt" ] && \
        cp "/data/Quantification/STAR/HT/HTSeq_Count_union.txt" "/data/DEG/S_HTSeq_Count_union.txt" && ((COPIED++))
    [ -f "/data/Quantification/STAR/RSEM/rsem_gene_counts_matrix.txt" ] && \
        cp "/data/Quantification/STAR/RSEM/rsem_gene_counts_matrix.txt" "/data/DEG/S_rsem_gene_counts_matrix.txt" && ((COPIED++))
    if [ $COPIED -gt 0 ]; then
        echo "      ✓ STAR: $COPIED matrix file(s) copied"
    else
        echo "      ⚠ STAR: no quantification output found to copy"
    fi
}

# ===================================================================
# HELPER FUNCTION: Clean count matrices for an aligner
# Uses file-existence checks — no dependency on *_FAILED variables.
# ===================================================================
clean_matrices() {
    local NAME="$1"
    local TRIGGER_FILE="$2"
    local CLEAN_SCRIPT="$3"
    local LOGFILE="$4"

    if [ -f "$TRIGGER_FILE" ]; then
        echo "    - Cleaning $NAME matrices..."
        if micromamba run -n rnaseq-cli bash "$CLEAN_SCRIPT" > "$LOGFILE" 2>&1; then
            echo "      ✓ $NAME cleaning completed"
        else
            echo "      ✗ $NAME cleaning failed! (see $LOGFILE)" >&2
            tail -20 "$LOGFILE" >&2
        fi
    else
        echo "    - Skipping $NAME cleaning (no trigger file: $TRIGGER_FILE)"
    fi
}

# ===================================================================
# HELPER FUNCTION: Distribute housekeeping_genes.xlsx into every
# pipeline subfolder that already exists under /data/results/DEG/
# Called: (1) after Stage 1 matrix cleaning
#         (2) after Stage 2 DEG analysis completes
# This ensures housekeeping_genes.xlsx sits alongside every
# pipeline's xlsx outputs for easy reference:
#   results/DEG/B_FC/edgeR/housekeeping_genes.xlsx
#   results/DEG/B_FC/DESeq2/housekeeping_genes.xlsx
#   results/DEG/S_RSEM/edgeR/housekeeping_genes.xlsx  ... etc.
# ===================================================================
distribute_hkg_to_subfolders() {
    local SRC=""

    # Find source — same priority order as SETUP section
    if   [ -f "/data/DEG/housekeeping_genes.xlsx" ];                            then SRC="/data/DEG/housekeeping_genes.xlsx"
    elif [ -f "/data/results/DEG/housekeeping_genes.xlsx" ];                    then SRC="/data/results/DEG/housekeeping_genes.xlsx"
    elif [ -f "/data/busco_output/filtered/housekeeping_genes.xlsx" ];          then SRC="/data/busco_output/filtered/housekeeping_genes.xlsx"
    elif [ -f "/data/housekeeping_genes.xlsx" ];                                then SRC="/data/housekeeping_genes.xlsx"
    fi

    if [ -z "$SRC" ]; then
        echo "    ⚠ housekeeping_genes.xlsx not found — skipping subfolder distribution"
        return 0
    fi

    # Walk every subdirectory under /data/results/DEG/ and copy
    local COUNT=0
    while IFS= read -r -d '' subdir; do
        cp "$SRC" "$subdir/housekeeping_genes.xlsx" 2>/dev/null && ((COUNT++)) || true
    done < <(find /data/results/DEG -mindepth 2 -maxdepth 2 -type d -print0 2>/dev/null)

    if [ $COUNT -gt 0 ]; then
        echo "    ✓ housekeeping_genes.xlsx copied into $COUNT pipeline subfolder(s)"
        echo "      under /data/results/DEG/"
    else
        echo "    ℹ No pipeline subfolders found yet under /data/results/DEG/"
        echo "      (will be distributed again after DEG analysis completes)"
    fi
}

# ===================================================================
# STAGE 1: ALIGNMENT & QUANTIFICATION
# ===================================================================
if [[ "$RUN_ALIGNMENT" == "true" ]]; then

    echo ""
    echo "=========================================="
    echo "STAGE 1: ALIGNMENT & QUANTIFICATION"
    echo "Environment: rnaseq-cli"
    echo "Execution: ${ALIGNMENT_MODE^^}"
    echo "=========================================="
    echo ""

    # ------------------------------------------------------------------
    # STEP 1: Quality Control
    # ------------------------------------------------------------------
    echo ">>> STEP 1: Quality Control (fastp + FastQC)"
    echo "  Started: $(date)"
    if ! micromamba run -n rnaseq-cli bash -c "cd /data && /opt/project/scripts/quality_control.sh" \
            > /data/logs/qc_${TIMESTAMP}.log 2>&1; then
        echo "ERROR: Quality control failed!" >&2
        tail -30 /data/logs/qc_${TIMESTAMP}.log >&2
        exit 1
    fi
    echo "  Status: ✓ Completed"

    # ------------------------------------------------------------------
    # STEP 2: Alignment
    # ------------------------------------------------------------------
    echo ""
    echo ">>> STEP 2: Running Aligners (mode: $ALIGNMENT_MODE)"
    echo ""

    BOWTIE_OK=0
    HISAT2_OK=0
    STAR_OK=0

    if [[ "$ALIGNMENT_MODE" == "sequential" || "$ALIGNMENT_MODE" == "sequential_only" ]]; then
        # ---- SEQUENTIAL ----
        run_aligner "Bowtie2" "/opt/project/scripts/bowtie.sh" \
            "/data/logs/bowtie_${TIMESTAMP}.log"  && BOWTIE_OK=1 || true

        echo ""
        run_aligner "HISAT2"  "/opt/project/scripts/hisat2.sh" \
            "/data/logs/hisat2_${TIMESTAMP}.log"  && HISAT2_OK=1 || true

        echo ""
        run_aligner "STAR"    "/opt/project/scripts/star.sh"   \
            "/data/logs/star_${TIMESTAMP}.log"    && STAR_OK=1   || true

    else
        # ---- PARALLEL ----
        run_aligner "Bowtie2" "/opt/project/scripts/bowtie.sh" \
            "/data/logs/bowtie_${TIMESTAMP}.log"  &
        BOWTIE_PID=$!

        run_aligner "HISAT2"  "/opt/project/scripts/hisat2.sh" \
            "/data/logs/hisat2_${TIMESTAMP}.log"  &
        HISAT2_PID=$!

        run_aligner "STAR"    "/opt/project/scripts/star.sh"   \
            "/data/logs/star_${TIMESTAMP}.log"    &
        STAR_PID=$!

        echo "    Waiting for all aligners to finish..."
        wait $BOWTIE_PID && BOWTIE_OK=1 || true
        wait $HISAT2_PID && HISAT2_OK=1 || true
        wait $STAR_PID   && STAR_OK=1   || true
    fi

    # Tally results
    SUCCESSFUL_ALIGNERS=$((BOWTIE_OK + HISAT2_OK + STAR_OK))
    echo ""
    echo "    Alignment summary: $SUCCESSFUL_ALIGNERS/3 aligners succeeded"
    [ $BOWTIE_OK -eq 1 ] && echo "      ✓ Bowtie2" || echo "      ✗ Bowtie2"
    [ $HISAT2_OK -eq 1 ] && echo "      ✓ HISAT2"  || echo "      ✗ HISAT2"
    [ $STAR_OK   -eq 1 ] && echo "      ✓ STAR"    || echo "      ✗ STAR"

    if [ $SUCCESSFUL_ALIGNERS -eq 0 ]; then
        echo "ERROR: All aligners failed — cannot continue." >&2
        exit 1
    fi

    # Verify BAM files
    BAM_COUNT=$(find /data/Mapping -name "*_sorted.bam" -type f 2>/dev/null | wc -l)
    if [ "$BAM_COUNT" -eq 0 ]; then
        echo "ERROR: No BAM files were created!" >&2
        exit 1
    fi
    echo "    Output: $BAM_COUNT BAM file(s) created"

    # ------------------------------------------------------------------
    # STEP 3: Copy count matrices to /data/DEG/
    # ------------------------------------------------------------------
    echo ""
    echo ">>> STEP 3: Copying count matrices to DEG directory..."

    # Copy sampleinfo
    if [ -f "/data/sampleinfo.txt" ]; then
        cp /data/sampleinfo.txt /data/DEG/sampleinfo.txt
        echo "    ✓ sampleinfo.txt copied"
    else
        echo "    ⚠ WARNING: /data/sampleinfo.txt not found"
    fi

    # Copy housekeeping_genes.xlsx — same priority logic as SETUP section
    # Priority 1: BUSCO pipeline output
    # Priority 2: /data/ root
    # Priority 3: already in DEG (skip re-copy, already distributed)
    if [ -f "/data/busco_output/filtered/housekeeping_genes.xlsx" ]; then
        cp /data/busco_output/filtered/housekeeping_genes.xlsx /data/DEG/housekeeping_genes.xlsx
        cp /data/busco_output/filtered/housekeeping_genes.xlsx /data/results/DEG/housekeeping_genes.xlsx
        echo "    ✓ housekeeping_genes.xlsx distributed (source: BUSCO output)"
    elif [ -f "/data/housekeeping_genes.xlsx" ]; then
        cp /data/housekeeping_genes.xlsx /data/DEG/housekeeping_genes.xlsx
        cp /data/housekeeping_genes.xlsx /data/results/DEG/housekeeping_genes.xlsx
        echo "    ✓ housekeeping_genes.xlsx distributed (source: /data/)"
    elif [ -f "/data/DEG/housekeeping_genes.xlsx" ]; then
        cp /data/DEG/housekeeping_genes.xlsx /data/results/DEG/housekeeping_genes.xlsx 2>/dev/null || true
        echo "    ✓ housekeeping_genes.xlsx synced to results/DEG/ (source: DEG dir)"
    else
        echo "    ⚠ housekeeping_genes.xlsx not found — HKG analysis will be skipped"
    fi

    # Copy per-aligner matrices — driven by file existence, not variables
    [ $BOWTIE_OK -eq 1 ] && copy_matrices_bowtie
    [ $HISAT2_OK -eq 1 ] && copy_matrices_hisat2
    [ $STAR_OK   -eq 1 ] && copy_matrices_star

    echo "    Status: ✓ Completed"

    # ------------------------------------------------------------------
    # STEP 4: Clean count matrices
    # ------------------------------------------------------------------
    echo ""
    echo ">>> STEP 4: Cleaning count matrices..."

    clean_matrices "Bowtie2" \
        "/data/DEG/B_FC_Count.txt" \
        "/opt/project/scripts/bowtie_clean.sh" \
        "/data/logs/bowtie_clean_${TIMESTAMP}.log"

    clean_matrices "HISAT2" \
        "/data/DEG/H_FC_Count.txt" \
        "/opt/project/scripts/hisat2_clean.sh" \
        "/data/logs/hisat2_clean_${TIMESTAMP}.log"

    clean_matrices "STAR" \
        "/data/DEG/S_FC_Count.txt" \
        "/opt/project/scripts/star_clean.sh" \
        "/data/logs/star_clean_${TIMESTAMP}.log"

    CLEAN_MATRICES=$(find /data/DEG -name "*_clean.txt" 2>/dev/null | wc -l)
    echo ""
    echo "    Output: $CLEAN_MATRICES cleaned count matrix/matrices"

    if [ $CLEAN_MATRICES -gt 0 ]; then
        echo "    Cleaned matrices:"
        find /data/DEG -name "*_clean.txt" -exec basename {} \; | sed 's/^/      - /'
    fi

    echo "    Status: ✓ Completed"

    # Distribute housekeeping_genes.xlsx into pipeline subfolders
    # (Stage 1 creates results/DEG/ structure via quantification)
    echo ""
    echo ">>> Distributing housekeeping_genes.xlsx to pipeline subfolders..."
    distribute_hkg_to_subfolders

else
    echo ""
    echo "=========================================="
    echo "STAGE 1: ALIGNMENT & QUANTIFICATION SKIPPED"
    echo "=========================================="
fi

# ===================================================================
# STAGE 2: DEG ANALYSIS
# ===================================================================
if [[ "$RUN_DEG" == "true" ]]; then

    echo ""
    echo "=========================================="
    echo "STAGE 2: DEG ANALYSIS"
    echo "Environment: rnaseq-r"
    echo "=========================================="
    echo "Thresholds: log2FC=$LOG2FC_THRESHOLD, FDR=$FDR_THRESHOLD, p-value=$PVALUE_THRESHOLD"
    echo ""

    if [ ! -f "/data/DEG/sampleinfo.txt" ]; then
        echo "ERROR: sampleinfo.txt not found in /data/DEG/" >&2
        exit 1
    fi
    echo "  [OK] sampleinfo.txt verified"

    # Check and redistribute housekeeping_genes.xlsx at DEG stage start
    # Handles --deg-only runs where Stage 1 was skipped
    # Same priority order as SETUP section
    if [ -f "/data/busco_output/filtered/housekeeping_genes.xlsx" ]; then
        cp /data/busco_output/filtered/housekeeping_genes.xlsx /data/DEG/housekeeping_genes.xlsx 2>/dev/null || true
        cp /data/busco_output/filtered/housekeeping_genes.xlsx /data/results/DEG/housekeeping_genes.xlsx 2>/dev/null || true
        echo "  [OK] housekeeping_genes.xlsx distributed (source: BUSCO output)"
    elif [ -f "/data/housekeeping_genes.xlsx" ]; then
        cp /data/housekeeping_genes.xlsx /data/DEG/housekeeping_genes.xlsx 2>/dev/null || true
        cp /data/housekeeping_genes.xlsx /data/results/DEG/housekeeping_genes.xlsx 2>/dev/null || true
        echo "  [OK] housekeeping_genes.xlsx distributed (source: /data/)"
    elif [ -f "/data/DEG/housekeeping_genes.xlsx" ]; then
        cp /data/DEG/housekeeping_genes.xlsx /data/results/DEG/housekeeping_genes.xlsx 2>/dev/null || true
        echo "  [OK] housekeeping_genes.xlsx found in DEG dir — synced to results/DEG/"
    else
        echo "  ⚠  housekeeping_genes.xlsx not found — HKG analysis will be skipped"
        echo "     Locations checked:"
        echo "       /data/busco_output/filtered/housekeeping_genes.xlsx"
        echo "       /data/housekeeping_genes.xlsx"
        echo "       /data/DEG/housekeeping_genes.xlsx"
    fi

    COUNT_FILES=$(find /data/DEG -name '*_clean.txt' 2>/dev/null | wc -l)
    if [ "$COUNT_FILES" -eq 0 ]; then
        echo "ERROR: No cleaned count matrices found in /data/DEG/" >&2
        echo "       Directory contents:" >&2
        ls -lh /data/DEG/ 2>/dev/null >&2 || echo "  (empty)" >&2
        exit 1
    fi
    echo "  [OK] Found $COUNT_FILES cleaned count matrix/matrices:"
    find /data/DEG -name '*_clean.txt' -exec basename {} \; | sed 's/^/    - /'

    # Run edgeR
    echo ""
    echo ">>> Running edgeR Analysis..."
    if ! micromamba run -n rnaseq-r bash -c "
        cd /data
        export SAMPLE_INFO_PATH='$SAMPLE_INFO_PATH'
        export LOG2FC_THRESHOLD='$LOG2FC_THRESHOLD'
        export FDR_THRESHOLD='$FDR_THRESHOLD'
        export PVALUE_THRESHOLD='$PVALUE_THRESHOLD'
        /opt/project/scripts/run_all_edgeR.sh
    " > /data/logs/edgeR_${TIMESTAMP}.log 2>&1; then
        echo "ERROR: edgeR failed!" >&2
        tail -30 /data/logs/edgeR_${TIMESTAMP}.log >&2
        exit 1
    fi
    echo "  ✓ edgeR completed successfully"

    # Run DESeq2
    echo ""
    echo ">>> Running DESeq2 Analysis..."
    if ! micromamba run -n rnaseq-r bash -c "
        cd /data
        export SAMPLE_INFO_PATH='$SAMPLE_INFO_PATH'
        export LOG2FC_THRESHOLD='$LOG2FC_THRESHOLD'
        export FDR_THRESHOLD='$FDR_THRESHOLD'
        export PVALUE_THRESHOLD='$PVALUE_THRESHOLD'
        /opt/project/scripts/run_all_deseq2.sh
    " > /data/logs/deseq2_${TIMESTAMP}.log 2>&1; then
        echo "ERROR: DESeq2 failed!" >&2
        tail -30 /data/logs/deseq2_${TIMESTAMP}.log >&2
        exit 1
    fi
    echo "  ✓ DESeq2 completed successfully"

    DEG_FILES=$(find /data/results/DEG -name "*.xlsx" 2>/dev/null | wc -l)
    echo ""
    echo "  Generated $DEG_FILES DEG result file(s)"

    # Distribute housekeeping_genes.xlsx into every pipeline subfolder
    # now that all DEG output directories have been created
    echo ""
    echo ">>> Distributing housekeeping_genes.xlsx to all pipeline subfolders..."
    distribute_hkg_to_subfolders

else
    echo ""
    echo "=========================================="
    echo "STAGE 2: DEG ANALYSIS SKIPPED"
    echo "=========================================="
fi

# ===================================================================
# STAGE 3: COMPARISON ANALYSIS
# ===================================================================
if [[ "$RUN_COMPARISON" == "true" ]]; then

    echo ""
    echo "=========================================="
    echo "STAGE 3: COMPARISON ANALYSIS"
    echo "Environment: rnaseq-comparison (Python)"
    echo "=========================================="

    DEG_COUNT=$(find /data/results/DEG \( -name "*_UP.xlsx" -o -name "*_DOWN.xlsx" \) 2>/dev/null | wc -l)

    if [ "$DEG_COUNT" -eq 0 ]; then
        echo "WARNING: No DEG results found — skipping comparison" >&2
    else
        echo "  [OK] Found $DEG_COUNT DEG result file(s)"
        echo ""
        echo ">>> Running Pipeline Comparison Analysis..."

        if ! micromamba run -n rnaseq-comparison bash -c "
            cd /data
            python /opt/project/scripts/run_pipeline_comparison.py
        " > /data/logs/comparison_${TIMESTAMP}.log 2>&1; then
            echo "ERROR: Comparison analysis failed!" >&2
            tail -30 /data/logs/comparison_${TIMESTAMP}.log >&2
            exit 1
        fi

        echo "  ✓ Comparison analysis completed successfully"
        COMP_FILES=$(find /data/results/pipeline_comparison -type f 2>/dev/null | wc -l)
        echo "  Generated $COMP_FILES comparison output file(s)"
    fi

else
    echo ""
    echo "=========================================="
    echo "STAGE 3: COMPARISON ANALYSIS SKIPPED"
    echo "=========================================="
fi

# ===================================================================
# FINAL SUMMARY
# ===================================================================
echo ""
echo "=========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "Completion Time: $(date)"
echo "=========================================="
echo ""
echo "Output Summary:"
echo "──────────────────────────────────────────"

if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    BAM_COUNT=$(find /data/Mapping -name "*.bam" 2>/dev/null | wc -l)
    CLEAN_COUNT=$(find /data/DEG -name "*_clean.txt" 2>/dev/null | wc -l)
    echo ""
    echo "Stage 1 - Alignment & Quantification:"
    echo "  ✓ Execution mode: ${ALIGNMENT_MODE^^}"
    if [ -n "$INDEX_DIR" ]; then
        echo "  ✓ Index source: PRE-BUILT ($INDEX_DIR)"
    else
        echo "  ✓ Index source: BUILT LOCALLY"
    fi
    echo "  ✓ BAM files created: $BAM_COUNT"
    echo "  ✓ Cleaned count matrices: $CLEAN_COUNT"
    echo "  ✓ QC reports: /data/qc/multiqc/"
fi

if [[ "$RUN_DEG" == "true" ]]; then
    DEG_UP=$(find /data/results/DEG -name "*_UP.xlsx" 2>/dev/null | wc -l)
    DEG_DOWN=$(find /data/results/DEG -name "*_DOWN.xlsx" 2>/dev/null | wc -l)
    DEG_TOTAL=$((DEG_UP + DEG_DOWN))
    echo ""
    echo "Stage 2 - Differential Expression:"
    echo "  ✓ UP-regulated gene lists: $DEG_UP"
    echo "  ✓ DOWN-regulated gene lists: $DEG_DOWN"
    echo "  ✓ Total DEG files: $DEG_TOTAL"
    echo "  ✓ Results location: /data/results/DEG/"
fi

if [[ "$RUN_COMPARISON" == "true" ]]; then
    COMP_FILES=$(find /data/results/pipeline_comparison -type f 2>/dev/null | wc -l)
    if [ "$COMP_FILES" -gt 0 ]; then
        echo ""
        echo "Stage 3 - Pipeline Comparison:"
        echo "  ✓ Comparison files: $COMP_FILES"
        echo "  ✓ Results location: /data/results/pipeline_comparison/"
    fi
fi

QC_FILES=$(find /data/results/QC_Summary -type f 2>/dev/null | wc -l)
if [ "$QC_FILES" -gt 0 ]; then
    echo ""
    echo "QC & Alignment Summary:"
    echo "  ✓ QC summary files: $QC_FILES"
    echo "  ✓ Results location: /data/results/QC_Summary/"
fi

echo ""
echo "Log Files:"
echo "  → All logs: /data/logs/"
echo "  → Latest run: /data/logs/*_${TIMESTAMP}.log"
echo ""
echo "Disk Space:"
df -h /data 2>/dev/null | tail -1 | awk '{print "  → Available: " $4 " / Used: " $3}' || echo "  → N/A"
echo ""
echo "──────────────────────────────────────────"
echo "✅ ALL PIPELINE STAGES COMPLETED"
echo "──────────────────────────────────────────"
echo ""
echo "Results are ready for downstream analysis!"
echo ""
