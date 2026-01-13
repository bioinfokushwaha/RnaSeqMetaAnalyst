#!/bin/bash
# ===================================================================
# RNA-Seq Pipeline Master Script - Complete Fixed Version
# ===================================================================
# This script orchestrates the complete RNA-seq analysis pipeline:
# 1. Quality Control (FastQC, fastp, MultiQC)
# 2. Alignment (STAR, HISAT2, Bowtie2)
# 3. Quantification (featureCounts, HTSeq, RSEM)
# 4. Count Matrix Cleaning
# 5. Differential Expression Analysis (DESeq2, edgeR)
# 6. Pipeline Comparison Analysis
# ===================================================================

set -e
set -o pipefail

echo "=========================================="
echo "RNA-Seq Pipeline Started: $(date)"
echo "=========================================="

# ===================================================================
# SET WORKING DIRECTORY
# ===================================================================
# CRITICAL: Change to /data directory and stay there
cd /data || { echo "ERROR: Cannot access /data directory"; exit 1; }
export WORK_DIR="/data"
echo "Working directory: $(pwd)"
echo ""

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
  --log2fc <value>     log2 Fold Change threshold (default: 0)
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
# EARLY SETUP: CREATE ALL DIRECTORIES INCLUDING DEG
# ===================================================================
echo ""
echo ">>> Step 1.1: Creating output directories..."
mkdir -p /data/DEG \
         /data/results/DEG \
         /data/Indices/Bowtie /data/Mapping/Bowtie /data/Quantification/Bowtie/{HT,FC,RSEM} \
         /data/Indices/Hisat2 /data/Mapping/Hisat2 /data/Quantification/Hisat2/{HT,FC} \
         /data/Indices/STAR /data/Mapping/STAR /data/Quantification/STAR/{HT,FC,RSEM} \
         /data/Indices/RSEM/{Bowtie,STAR} \
         /data/qc/{fastqc_raw,fastqc_clean,fastp_reports,multiqc} \
         /data/Trim /data/logs /data/pipeline_comparison

chmod -R 777 /data/DEG /data/results /data/Indices /data/Mapping /data/Quantification /data/qc /data/Trim /data/logs /data/pipeline_comparison 2>/dev/null || true
echo "  [OK] All directories created (including /data/DEG)"

# ===================================================================
# PREPARE SAMPLEINFO FOR DEG ANALYSIS
# ===================================================================
if [[ "$RUN_DEG" == "true" ]]; then
    echo ""
    echo ">>> Step 1.2: Preparing sampleinfo.txt for DEG analysis..."
    
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
    elif [ -f "/data/DEG/sampleinfo.txt" ]; then
        echo "  [OK] Found in DEG directory from previous run"
        SAMPLE_INFO_SOURCE="/data/DEG/sampleinfo.txt"
        SAMPLEINFO_FOUND=true
    fi
    
    # Copy sampleinfo.txt to DEG directory if found
    if [ "$SAMPLEINFO_FOUND" = true ] && [ -n "$SAMPLE_INFO_SOURCE" ]; then
        cp "$SAMPLE_INFO_SOURCE" "/data/DEG/sampleinfo.txt"
        export SAMPLE_INFO_PATH="/data/DEG/sampleinfo.txt"
        echo "  [OK] Copied to: /data/DEG/sampleinfo.txt"
        
        # Verify the copy
        if [ ! -f "/data/DEG/sampleinfo.txt" ]; then
            echo "ERROR: Failed to copy sampleinfo.txt" >&2
            exit 1
        fi
        
        # Show contents for verification
        echo ""
        echo "  Sample info content:"
        head -3 "/data/DEG/sampleinfo.txt" | sed 's/^/    /'
        
    else
        echo "WARNING: sampleinfo.txt not found in any location:"
        echo "  Checked: /data/sampleinfo.txt"
        echo "  Checked: /data/data/sampleinfo.txt"
        echo "  Checked: /data/DEG/sampleinfo.txt"
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
        bash -c "cd /data && /opt/project/scripts/quality_control.sh"
    
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
    
    micromamba run -n rnaseq-cli bash -c "
        cd /data
        export THREADS='$THREADS'
        export MODE='$MODE'
        export GENOME_DIR='$GENOME_DIR'
        export READ_DIR='$READ_DIR'
        export TRIM_DIR='$TRIM_DIR'
        export GTF='$GTF'
        export FASTA='$FASTA'
        export GFF='$GFF'
        export INDEX_DIR='$INDEX_DIR'
        
        parallel -j 3 --halt soon,fail=1 ::: \
            '/opt/project/scripts/bowtie.sh > /data/logs/bowtie_${TIMESTAMP}.log 2>&1' \
            '/opt/project/scripts/hisat2.sh > /data/logs/hisat2_${TIMESTAMP}.log 2>&1' \
            '/opt/project/scripts/star.sh > /data/logs/star_${TIMESTAMP}.log 2>&1'
    "

    if [ $? -ne 0 ]; then
        echo "ERROR: Alignment failed!" >&2
        echo "Check logs at: /data/logs/" >&2
        exit 1
    fi
    echo "  [OK] All aligners completed"

    # Copy Quantification Outputs to DEG Directory
    echo ""
    echo ">>> Step 1.3: Copying Quantification Outputs to DEG directory"
    
    # STAR outputs
    [ -f "/data/Quantification/STAR/FC/FC_Count.txt" ] && \
        cp /data/Quantification/STAR/FC/FC_Count.txt /data/DEG/S_FC_Count.txt && \
        echo "  ✓ Copied STAR FC"
    
    [ -f "/data/Quantification/STAR/HT/HTSeq_Count_union.txt" ] && \
        cp /data/Quantification/STAR/HT/HTSeq_Count_union.txt /data/DEG/S_HTSeq_Count_union.txt && \
        echo "  ✓ Copied STAR HTSeq"
    
    [ -f "/data/Quantification/STAR/RSEM/rsem_gene_counts_matrix.txt" ] && \
        cp /data/Quantification/STAR/RSEM/rsem_gene_counts_matrix.txt /data/DEG/S_rsem_gene_counts_matrix.txt && \
        echo "  ✓ Copied STAR RSEM"
    
    # HISAT2 outputs
    [ -f "/data/Quantification/Hisat2/FC/FC_Count.txt" ] && \
        cp /data/Quantification/Hisat2/FC/FC_Count.txt /data/DEG/H_FC_Count.txt && \
        echo "  ✓ Copied HISAT2 FC"
    
    [ -f "/data/Quantification/Hisat2/HT/HTSeq_Count_union.txt" ] && \
        cp /data/Quantification/Hisat2/HT/HTSeq_Count_union.txt /data/DEG/H_HTSeq_Count_union.txt && \
        echo "  ✓ Copied HISAT2 HTSeq"
    
    # Bowtie2 outputs
    [ -f "/data/Quantification/Bowtie/FC/FC_Count.txt" ] && \
        cp /data/Quantification/Bowtie/FC/FC_Count.txt /data/DEG/B_FC_Count.txt && \
        echo "  ✓ Copied Bowtie2 FC"
    
    [ -f "/data/Quantification/Bowtie/HT/HTSeq_Count_union.txt" ] && \
        cp /data/Quantification/Bowtie/HT/HTSeq_Count_union.txt /data/DEG/B_HTSeq_Count_union.txt && \
        echo "  ✓ Copied Bowtie2 HTSeq"
    
    [ -f "/data/Quantification/Bowtie/RSEM/rsem_gene_counts_matrix.txt" ] && \
        cp /data/Quantification/Bowtie/RSEM/rsem_gene_counts_matrix.txt /data/DEG/B_rsem_gene_counts_matrix.txt && \
        echo "  ✓ Copied Bowtie2 RSEM"
    
    # Verify at least some files were copied
    COPIED_COUNT=$(find /data/DEG -name "*.txt" -type f 2>/dev/null | wc -l)
    if [ "$COPIED_COUNT" -eq 0 ]; then
        echo "WARNING: No count files were copied to DEG directory!" >&2
    else
        echo "  ✓ Total files copied: $COPIED_COUNT"
    fi
    
    # Count Matrix Cleaning
    echo ""
    echo ">>> Step 1.4: Cleaning Count Matrices"
    
    # Run each cleaning script with explicit working directory
    echo "  - Cleaning Bowtie2 results..."
    micromamba run -n rnaseq-cli bash -c "cd /data && /opt/project/scripts/bowtie_clean.sh"
    
    echo "  - Cleaning HISAT2 results..."
    micromamba run -n rnaseq-cli bash -c "cd /data && /opt/project/scripts/hisat2_clean.sh"
    
    echo "  - Cleaning STAR results..."
    micromamba run -n rnaseq-cli bash -c "cd /data && /opt/project/scripts/star_clean.sh"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Cleaning failed!" >&2
        exit 1
    fi
    
    echo "  ✓ All count matrices cleaned"
    
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
    
    # Verify sampleinfo.txt exists in DEG directory
    if [ ! -f "/data/DEG/sampleinfo.txt" ]; then
        echo "ERROR: sampleinfo.txt not found in DEG directory!" >&2
        echo "  Expected: /data/DEG/sampleinfo.txt" >&2
        exit 1
    fi
    echo "  [OK] sampleinfo.txt verified in DEG directory"
    
    # Count available matrices
    COUNT_FILES=$(find /data/DEG -name '*_clean.txt' 2>/dev/null | wc -l)
    if [ "$COUNT_FILES" -eq 0 ]; then
        echo "ERROR: No count matrices found in /data/DEG/" >&2
        echo ""
        echo "Available files in /data/DEG/:"
        ls -lh /data/DEG/ 2>/dev/null || echo "  Directory is empty or inaccessible"
        exit 1
    fi
    echo "  [OK] Found $COUNT_FILES count matrices"
    
    # List available matrices
    echo ""
    echo "  Clean count matrices available:"
    find /data/DEG -name '*_clean.txt' -exec basename {} \; | sed 's/^/    - /'
    
    # Run edgeR
    echo ""
    echo ">>> Running edgeR Analysis..."
    micromamba run -n rnaseq-r bash -c "
        cd /data
        export SAMPLE_INFO_PATH='$SAMPLE_INFO_PATH'
        export LOG2FC_THRESHOLD='$LOG2FC_THRESHOLD'
        export FDR_THRESHOLD='$FDR_THRESHOLD'
        export PVALUE_THRESHOLD='$PVALUE_THRESHOLD'
        /opt/project/scripts/run_all_edgeR.sh
    " > /data/logs/edgeR_${TIMESTAMP}.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: edgeR failed!" >&2
        echo ""
        echo "Last 30 lines of edgeR log:"
        tail -30 /data/logs/edgeR_${TIMESTAMP}.log >&2
        exit 1
    fi
    echo "  ✓ edgeR completed successfully"
    
    # Run DESeq2
    echo ""
    echo ">>> Running DESeq2 Analysis..."
    micromamba run -n rnaseq-r bash -c "
        cd /data
        export SAMPLE_INFO_PATH='$SAMPLE_INFO_PATH'
        export LOG2FC_THRESHOLD='$LOG2FC_THRESHOLD'
        export FDR_THRESHOLD='$FDR_THRESHOLD'
        export PVALUE_THRESHOLD='$PVALUE_THRESHOLD'
        /opt/project/scripts/run_all_deseq2.sh
    " > /data/logs/deseq2_${TIMESTAMP}.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: DESeq2 failed!" >&2
        echo ""
        echo "Last 30 lines of DESeq2 log:"
        tail -30 /data/logs/deseq2_${TIMESTAMP}.log >&2
        exit 1
    fi
    echo "  ✓ DESeq2 completed successfully"
    
    # Verify DEG results were created
    DEG_FILES=$(find /data/results/DEG -name "*.xlsx" 2>/dev/null | wc -l)
    echo ""
    echo "  Generated $DEG_FILES DEG result files"
    
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
    
    # Check for DEG results
    DEG_COUNT=$(find /data/results/DEG -name "*_UP.xlsx" -o -name "*_DOWN.xlsx" 2>/dev/null | wc -l)
    
    if [ "$DEG_COUNT" -eq 0 ]; then
        echo "WARNING: No DEG results found - skipping comparison" >&2
        echo ""
        echo "Expected location: /data/results/DEG/"
        echo "Files looked for: *_UP.xlsx, *_DOWN.xlsx"
    else
        echo "  [OK] Found $DEG_COUNT DEG result files"
        
        echo ""
        echo ">>> Running Pipeline Comparison Analysis..."
        micromamba run -n rnaseq-comparison bash -c "
            cd /data
            python /opt/project/scripts/run_pipeline_comparison.py
        " > /data/logs/comparison_${TIMESTAMP}.log 2>&1
        
        if [ $? -ne 0 ]; then
            echo "ERROR: Comparison analysis failed!" >&2
            echo ""
            echo "Last 30 lines of comparison log:"
            tail -30 /data/logs/comparison_${TIMESTAMP}.log >&2
            exit 1
        fi
        echo "  ✓ Comparison analysis completed successfully"
        
        # Verify comparison outputs
        COMP_FILES=$(find /data/results/pipeline_comparison -type f 2>/dev/null | wc -l)
        echo "  Generated $COMP_FILES comparison output files"
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
echo "=========================================="
echo "Completion Time: $(date)"
echo ""

# Generate comprehensive summary
echo "Output Summary:"
echo "──────────────────────────────────────────"

if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    BAM_COUNT=$(find /data/Mapping -name "*.bam" 2>/dev/null | wc -l)
    CLEAN_COUNT=$(find /data/DEG -name "*_clean.txt" 2>/dev/null | wc -l)
    echo ""
    echo "Stage 1 - Alignment & Quantification:"
    echo "  ✓ BAM files created: $BAM_COUNT"
    echo "  ✓ Clean count matrices: $CLEAN_COUNT"
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

# Show log file locations
echo ""
echo "Log Files:"
echo "  → All logs: /data/logs/"
echo "  → Latest run: /data/logs/*_${TIMESTAMP}.log"

# Final status
echo ""
echo "──────────────────────────────────────────"
echo "✅ ALL PIPELINE STAGES COMPLETED"
echo "──────────────────────────────────────────"
echo ""
echo "Results are ready for downstream analysis!"
echo ""
