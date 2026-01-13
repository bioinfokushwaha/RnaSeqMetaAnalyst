#!/bin/bash
# ===================================================================
# FIXED: master.sh - DEG Analysis Fix
# ===================================================================
# KEY FIX: Ensure sampleinfo.txt is properly copied to DEG directory
# before running R scripts, even on second runs (--deg-only)
# ===================================================================

set -e
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

export LOG2FC_THRESHOLD="${LOG2FC_THRESHOLD:-0}"
export FDR_THRESHOLD="${FDR_THRESHOLD:-0.05}"
export PVALUE_THRESHOLD="${PVALUE_THRESHOLD:-0.05}"
export STAR_SJDB_OVERHANG="${STAR_SJDB_OVERHANG:-100}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-alignment)
            RUN_ALIGNMENT="false"
            shift
            ;;
        --skip-deg)
            RUN_DEG="false"
            RUN_COMPARISON="false"
            shift
            ;;
        --skip-comparison)
            RUN_COMPARISON="false"
            shift
            ;;
        --deg-only)
            RUN_ALIGNMENT="false"
            RUN_DEG="true"
            shift
            ;;
        --comparison-only)
            RUN_ALIGNMENT="false"
            RUN_DEG="false"
            RUN_COMPARISON="true"
            shift
            ;;
        --deg-and-comparison)
            RUN_ALIGNMENT="false"
            RUN_DEG="true"
            RUN_COMPARISON="true"
            shift
            ;;
        --log2fc)
            export LOG2FC_THRESHOLD="$2"
            shift 2
            ;;
        --fdr)
            export FDR_THRESHOLD="$2"
            shift 2
            ;;
        --pvalue)
            export PVALUE_THRESHOLD="$2"
            shift 2
            ;;
        --sjdb-overhang)
            export STAR_SJDB_OVERHANG="$2"
            shift 2
            ;;
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

DEG Thresholds:
  --log2fc <value>     log2 Fold Change threshold (default: 1)
  --fdr <value>        FDR threshold (default: 0.05)
  --pvalue <value>     p-value threshold (default: 0.05)

Alignment Options:
  --sjdb-overhang <value> STAR parameter (default: 100)

Note: sampleinfo.txt is automatically copied from /data/ to DEG/ directory
EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ===================================================================
# VALIDATE ENVIRONMENT VARIABLES (only if running alignment)
# ===================================================================
if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    echo ">>> Validating environment variables..."
    : "${THREADS:?THREADS is not set}"
    : "${MODE:?MODE is not set (should be SE or PE)}"
    : "${GENOME_DIR:?GENOME_DIR is not set}"
    : "${READ_DIR:?READ_DIR is not set}"
    : "${TRIM_DIR:?TRIM_DIR is not set}"
    : "${GTF:?GTF is not set}"
    : "${FASTA:?FASTA is not set}"
    : "${GFF:?GFF is not set}"

    MODE=$(echo "$MODE" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')
    if [[ "$MODE" != "SE" && "$MODE" != "PE" ]]; then
        echo "ERROR: MODE must be 'SE' or 'PE', got: $MODE" >&2
        exit 1
    fi

    # Validate critical files
    echo ">>> Checking for critical files..."
    for file in "$GTF" "$FASTA" "$GFF"; do
        if [ ! -f "$file" ]; then
            echo "ERROR: Required file not found: $file" >&2
            exit 1
        fi
        echo "  [OK] Found: $file"
    done

    if [ ! -d "$READ_DIR" ] || [ -z "$(ls -A $READ_DIR)" ]; then
        echo "ERROR: READ_DIR is empty or doesn't exist: $READ_DIR" >&2
        exit 1
    fi
    echo "  [OK] Found reads in: $READ_DIR"
fi

# ===================================================================
# KEY FIX: PREPARE SAMPLEINFO FOR DEG ANALYSIS
# ===================================================================
# This ensures sampleinfo.txt is in the DEG directory before R scripts run
# This works for BOTH first-time runs and --deg-only reruns
# ===================================================================

if [[ "$RUN_DEG" == "true" ]]; then
    echo ""
    echo ">>> Preparing sampleinfo.txt for DEG analysis..."
    
    # Try multiple locations for sampleinfo.txt
    SAMPLEINFO_FOUND=false
    SAMPLE_INFO_SOURCE=""
    
    # Location 1: Mounted as /data/sampleinfo.txt
    if [ -f "/data/sampleinfo.txt" ]; then
        SAMPLE_INFO_SOURCE="/data/sampleinfo.txt"
        SAMPLEINFO_FOUND=true
        echo "  [OK] Found at: /data/sampleinfo.txt"
    
    # Location 2: Mounted as /data/data/sampleinfo.txt (nested)
    elif [ -f "/data/data/sampleinfo.txt" ]; then
        SAMPLE_INFO_SOURCE="/data/data/sampleinfo.txt"
        SAMPLEINFO_FOUND=true
        echo "  [OK] Found at: /data/data/sampleinfo.txt"
    
    # Location 3: Already in DEG directory (from previous run)
    elif [ -f "./DEG/sampleinfo.txt" ]; then
        echo "  [OK] Found in DEG directory from previous run"
        SAMPLE_INFO_SOURCE="./DEG/sampleinfo.txt"
        SAMPLEINFO_FOUND=true
    fi
    
    # Copy sampleinfo.txt to DEG directory if found
    if [ "$SAMPLEINFO_FOUND" = true ] && [ -n "$SAMPLE_INFO_SOURCE" ]; then
        mkdir -p ./DEG
        cp "$SAMPLE_INFO_SOURCE" "./DEG/sampleinfo.txt"
        export SAMPLE_INFO_PATH="./DEG/sampleinfo.txt"
        echo "  [OK] Copied to: ./DEG/sampleinfo.txt"
        
        # Verify the copy
        if [ ! -f "./DEG/sampleinfo.txt" ]; then
            echo "ERROR: Failed to copy sampleinfo.txt" >&2
            exit 1
        fi
        
        # Show contents for verification
        echo ""
        echo "  Sample info content:"
        head -3 "./DEG/sampleinfo.txt" | sed 's/^/    /'
        
    else
        echo "WARNING: sampleinfo.txt not found in any location:"
        echo "  Checked: /data/sampleinfo.txt"
        echo "  Checked: /data/data/sampleinfo.txt"
        echo "  Checked: ./DEG/sampleinfo.txt"
        echo ""
        echo "  DEG analysis will not be available"
        RUN_DEG="false"
        RUN_COMPARISON="false"
    fi
fi

# ===================================================================
# EXPORT VARIABLES FOR SUBSHELLS
# ===================================================================
export THREADS MODE GENOME_DIR READ_DIR GTF FASTA GFF TRIM_DIR INDEX_DIR
export SAMPLE_INFO_PATH LOG2FC_THRESHOLD FDR_THRESHOLD PVALUE_THRESHOLD
export STAR_SJDB_OVERHANG

# ===================================================================
# SETUP OUTPUT DIRECTORIES
# ===================================================================
if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    echo ""
    echo ">>> Creating output directories..."
    mkdir -p ./DEG \
             ./Indices/Bowtie ./Mapping/Bowtie ./Quantification/Bowtie/{HT,FC,RSEM} \
             ./Indices/Hisat2 ./Mapping/Hisat2 ./Quantification/Hisat2/{HT,FC} \
             ./Indices/STAR ./Mapping/STAR ./Quantification/STAR/{HT,FC,RSEM} \
             ./qc/{fastqc_raw,fastqc_clean,fastp_reports,multiqc} \
             ./Trim ./logs ./pipeline_comparison
    
    chmod -R 777 ./DEG ./Indices ./Mapping ./Quantification ./qc ./Trim ./logs ./pipeline_comparison 2>/dev/null || true
    echo "  [OK] Directories created"
fi

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# ===================================================================
# STAGE 1: ALIGNMENT & QUANTIFICATION
# ===================================================================
if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    
    echo ""
    echo "=========================================="
    echo "STAGE 1: ALIGNMENT & QUANTIFICATION"
    echo "Environment: rnaseq-cli"
    echo "=========================================="
    
    # Quality Control
    echo ""
    echo ">>> Running Quality Control..."
    micromamba run -n rnaseq-cli \
        --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
        --env TRIM_DIR --env GTF --env FASTA --env GFF \
        /opt/project/scripts/quality_control.sh
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Quality control failed!" >&2
        exit 1
    fi
    echo "  [OK] Quality control completed"

    TRIM_COUNT=$(find "$TRIM_DIR" -name "*.fastq" | wc -l)
    if [ $TRIM_COUNT -eq 0 ]; then
        echo "ERROR: No trimmed FASTQ files found" >&2
        exit 1
    fi
    echo "  [OK] Found $TRIM_COUNT trimmed FASTQ files"

    # Alignment & Quantification (Parallel)
    echo ""
    echo ">>> Running Alignment & Quantification (Parallel)..."
    
    micromamba run -n rnaseq-cli parallel -j 3 --halt soon,fail=1 ::: \
        "/opt/project/scripts/bowtie.sh > logs/bowtie_${TIMESTAMP}.log 2>&1" \
        "/opt/project/scripts/hisat2.sh > logs/hisat2_${TIMESTAMP}.log 2>&1" \
        "/opt/project/scripts/star.sh > logs/star_${TIMESTAMP}.log 2>&1"

    if [ $? -ne 0 ]; then
        echo "ERROR: Alignment failed!" >&2
        exit 1
    fi
    echo "  [OK] All aligners completed"

    # Count Matrix Cleaning
    echo ""
    echo ">>> Cleaning Count Matrices..."
    
    cd DEG || exit 1
    
    micromamba run -n rnaseq-cli /opt/project/scripts/bowtie_clean.sh
    micromamba run -n rnaseq-cli /opt/project/scripts/hisat2_clean.sh
    micromamba run -n rnaseq-cli /opt/project/scripts/star_clean.sh
    
    cd ..
    
    echo "  [OK] Count matrices cleaned"
    
else
    echo ""
    echo "=========================================="
    echo "STAGE 1: ALIGNMENT & QUANTIFICATION SKIPPED"
    echo "=========================================="
fi

# ===================================================================
# STAGE 2: DEG ANALYSIS (WITH FIX)
# ===================================================================
if [[ "$RUN_DEG" == "true" ]]; then
    
    echo ""
    echo "=========================================="
    echo "STAGE 2: DEG ANALYSIS"
    echo "Environment: rnaseq-r"
    echo "=========================================="
    echo "Thresholds: log2FC=$LOG2FC_THRESHOLD, FDR=$FDR_THRESHOLD"
    echo ""
    
    # KEY FIX: Verify sampleinfo.txt exists in DEG directory
    if [ ! -f "./DEG/sampleinfo.txt" ]; then
        echo "ERROR: sampleinfo.txt not found in DEG directory!" >&2
        echo "  Expected: ./DEG/sampleinfo.txt" >&2
        exit 1
    fi
    echo "  [OK] sampleinfo.txt verified in DEG directory"
    
    # Count available matrices
    COUNT_FILES=$(find ./DEG -name '*_clean.txt' | wc -l)
    if [ "$COUNT_FILES" -eq 0 ]; then
        echo "ERROR: No count matrices found in ./DEG/" >&2
        exit 1
    fi
    echo "  [OK] Found $COUNT_FILES count matrices"
    
    # Run edgeR
    echo ""
    echo ">>> Running edgeR..."
    micromamba run -n rnaseq-r \
        --env SAMPLE_INFO_PATH \
        --env LOG2FC_THRESHOLD \
        --env FDR_THRESHOLD \
        --env PVALUE_THRESHOLD \
        /opt/project/scripts/run_all_edgeR.sh \
        > logs/edgeR_${TIMESTAMP}.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: edgeR failed!" >&2
        tail -20 logs/edgeR_${TIMESTAMP}.log >&2
        exit 1
    fi
    echo "  [OK] edgeR completed"
    
    # Run DESeq2
    echo ""
    echo ">>> Running DESeq2..."
    micromamba run -n rnaseq-r \
        --env SAMPLE_INFO_PATH \
        --env LOG2FC_THRESHOLD \
        --env FDR_THRESHOLD \
        --env PVALUE_THRESHOLD \
        /opt/project/scripts/run_all_deseq2.sh \
        > logs/deseq2_${TIMESTAMP}.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: DESeq2 failed!" >&2
        tail -20 logs/deseq2_${TIMESTAMP}.log >&2
        exit 1
    fi
    echo "  [OK] DESeq2 completed"
    
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
    
    DEG_COUNT=$(find ./DEG -name "*_UP.xlsx" -o -name "*_DOWN.xlsx" | wc -l)
    if [ "$DEG_COUNT" -eq 0 ]; then
        echo "WARNING: No DEG results found - skipping comparison" >&2
    else
        echo "  [OK] Found $DEG_COUNT DEG result files"
        
        echo ""
        echo ">>> Running comparison analysis..."
        micromamba run -n rnaseq-comparison \
            python /opt/project/scripts/run_pipeline_comparison.py \
            > logs/comparison_${TIMESTAMP}.log 2>&1
        
        if [ $? -ne 0 ]; then
            echo "ERROR: Comparison failed!" >&2
            exit 1
        fi
        echo "  [OK] Comparison completed"
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
echo "Pipeline Completed: $(date)"
echo "=========================================="
echo ""

if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    BAM_COUNT=$(find ./Mapping -name "*.bam" 2>/dev/null | wc -l)
    CLEAN_COUNT=$(find ./DEG -name "*_clean.txt" 2>/dev/null | wc -l)
    echo "Stage 1 Output:"
    echo "  - BAM files: $BAM_COUNT"
    echo "  - Clean count matrices: $CLEAN_COUNT"
fi

if [[ "$RUN_DEG" == "true" ]]; then
    DEG_UP=$(find ./DEG -name "*_UP.xlsx" 2>/dev/null | wc -l)
    DEG_DOWN=$(find ./DEG -name "*_DOWN.xlsx" 2>/dev/null | wc -l)
    echo "Stage 2 Output:"
    echo "  - DEG files: $((DEG_UP + DEG_DOWN))"
fi

if [[ "$RUN_COMPARISON" == "true" ]]; then
    COMP_FILES=$(find ./pipeline_comparison -type f 2>/dev/null | wc -l)
    echo "Stage 3 Output:"
    echo "  - Comparison files: $COMP_FILES"
fi

echo ""
echo "[OK] Pipeline execution completed successfully!"
echo ""