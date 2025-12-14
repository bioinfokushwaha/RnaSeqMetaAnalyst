#!/bin/bash

set -e
set -o pipefail
set -x

echo "=========================================="
echo "RNA-Seq Pipeline Started: $(date)"
echo "=========================================="

# === Validate Required Environment Variables ===
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

# === Validate Critical Files ===
echo ">>> Checking for critical files..."
for file in "$GTF" "$FASTA" "$GFF"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file" >&2
        exit 1
    fi
    echo "  ✓ Found: $file"
done

if [ ! -d "$READ_DIR" ] || [ -z "$(ls -A $READ_DIR)" ]; then
    echo "ERROR: READ_DIR is empty or doesn't exist: $READ_DIR" >&2
    exit 1
fi
echo "  ✓ Found reads in: $READ_DIR"

SAMPLE_INFO_PATH="/data/sampleinfo.txt"
if [ ! -f "$SAMPLE_INFO_PATH" ]; then
    echo "ERROR: Sample info file not found: $SAMPLE_INFO_PATH" >&2
    exit 1
fi
echo "  ✓ Found sample info: $SAMPLE_INFO_PATH"

# === Export Variables ===
export THREADS MODE GENOME_DIR READ_DIR GTF FASTA GFF TRIM_DIR INDEX_DIR SAMPLE_INFO_PATH

# === Setup Output Directories ===
echo ">>> Creating output directories..."
mkdir -p ./DEG \
         ./Indices/Bowtie \
         ./Mapping/Bowtie \
         ./Quantification/Bowtie/HT \
         ./Quantification/Bowtie/FC \
         ./Quantification/Bowtie/RSEM \
         ./Indices/Hisat2 \
         ./Mapping/Hisat2 \
         ./Quantification/Hisat2/HT \
         ./Quantification/Hisat2/FC \
         ./Indices/STAR \
         ./Mapping/STAR \
         ./Quantification/STAR/HT \
         ./Quantification/STAR/FC \
         ./Quantification/STAR/RSEM \
         ./qc/fastqc_raw \
         ./qc/fastqc_clean \
         ./qc/fastp_reports \
         ./qc/multiqc \
         ./Trim

chmod -R 777 ./DEG ./Indices ./Mapping ./Quantification ./qc ./Trim
echo "  ✓ Directories created and permissions set"

# === Step 1: Quality Control ===
echo ""
echo "=========================================="
echo "STEP 1: Quality Control"
echo "=========================================="
micromamba run -n rnaseq-cli \
    --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
    --env TRIM_DIR --env GTF --env FASTA --env GFF \
    /opt/project/scripts/quality_control.sh

if [ $? -ne 0 ]; then
    echo "ERROR: Quality control failed!" >&2
    exit 1
fi
echo "  ✓ Quality control completed successfully"

TRIM_COUNT=$(find "$TRIM_DIR" -name "*.fastq" | wc -l)
if [ $TRIM_COUNT -eq 0 ]; then
    echo "ERROR: No trimmed FASTQ files found in $TRIM_DIR" >&2
    exit 1
fi
echo "  ✓ Found $TRIM_COUNT trimmed FASTQ files"

# === Step 2: Alignment & Quantification (Parallel) ===
echo ""
echo "=========================================="
echo "STEP 2: Alignment & Quantification (Parallel)"
echo "=========================================="

mkdir -p ./logs
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo ">>> Starting HISAT2 pipeline..."
micromamba run -n rnaseq-cli \
    --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
    --env TRIM_DIR --env GTF --env FASTA --env GFF --env INDEX_DIR \
    /opt/project/scripts/hisat2.sh > ./logs/hisat2_${TIMESTAMP}.log 2>&1 &
HISAT2_PID=$!

echo ">>> Starting Bowtie2 pipeline..."
micromamba run -n rnaseq-cli \
    --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
    --env TRIM_DIR --env GTF --env FASTA --env GFF --env INDEX_DIR \
    /opt/project/scripts/bowtie.sh > ./logs/bowtie_${TIMESTAMP}.log 2>&1 &
BOWTIE_PID=$!

echo ">>> Starting STAR pipeline..."
micromamba run -n rnaseq-cli \
    --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
    --env TRIM_DIR --env GTF --env FASTA --env GFF --env INDEX_DIR \
    /opt/project/scripts/star.sh > ./logs/star_${TIMESTAMP}.log 2>&1 &
STAR_PID=$!

echo ""
echo ">>> Waiting for alignment jobs to complete..."
echo "  PIDs: HISAT2=$HISAT2_PID, Bowtie2=$BOWTIE_PID, STAR=$STAR_PID"

FAILED_JOBS=""

wait $HISAT2_PID
HISAT2_EXIT=$?
if [ $HISAT2_EXIT -eq 0 ]; then
    echo "  ✓ HISAT2 completed successfully"
else
    echo "  ✗ HISAT2 failed with exit code $HISAT2_EXIT" >&2
    FAILED_JOBS="$FAILED_JOBS HISAT2"
fi

wait $BOWTIE_PID
BOWTIE_EXIT=$?
if [ $BOWTIE_EXIT -eq 0 ]; then
    echo "  ✓ Bowtie2 completed successfully"
else
    echo "  ✗ Bowtie2 failed with exit code $BOWTIE_EXIT" >&2
    FAILED_JOBS="$FAILED_JOBS Bowtie2"
fi

wait $STAR_PID
STAR_EXIT=$?
if [ $STAR_EXIT -eq 0 ]; then
    echo "  ✓ STAR completed successfully"
else
    echo "  ✗ STAR failed with exit code $STAR_EXIT" >&2
    FAILED_JOBS="$FAILED_JOBS STAR"
fi

if [ -n "$FAILED_JOBS" ]; then
    echo ""
    echo "ERROR: The following alignment jobs failed:$FAILED_JOBS" >&2
    echo "Check log files in ./logs/ for details" >&2
    exit 1
fi

echo ""
echo "  ✓ All alignment and quantification jobs completed successfully"

# ===================================================================
# === Step 2.5: Copy Quantification Outputs to DEG Directory ===
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 2.5: Copying Quantification Outputs"
echo "=========================================="

echo ">>> Copying count matrices to DEG directory..."

# STAR outputs (FC + HTSeq + RSEM)
if [ -f "Quantification/STAR/FC/FC_Count.txt" ]; then
    cp Quantification/STAR/FC/FC_Count.txt DEG/S_FC_Count.txt
    echo "  ✓ Copied STAR featureCounts"
else
    echo "  ⚠ STAR FC output not found"
fi

if [ -f "Quantification/STAR/HT/HTSeq_Count_union.txt" ]; then
    cp Quantification/STAR/HT/HTSeq_Count_union.txt DEG/S_HTSeq_Count_union.txt
    echo "  ✓ Copied STAR HTSeq"
else
    echo "  ⚠ STAR HTSeq output not found"
fi

if [ -f "Quantification/STAR/RSEM/rsem_gene_counts_matrix.txt" ]; then
    cp Quantification/STAR/RSEM/rsem_gene_counts_matrix.txt DEG/S_rsem_gene_counts_matrix.txt
    echo "  ✓ Copied STAR RSEM"
else
    echo "  ⚠ STAR RSEM output not found"
fi

# HISAT2 outputs (FC + HTSeq + RSEM)
if [ -f "Quantification/Hisat2/FC/FC_Count.txt" ]; then
    cp Quantification/Hisat2/FC/FC_Count.txt DEG/H_FC_Count.txt
    echo "  ✓ Copied HISAT2 featureCounts"
else
    echo "  ⚠ HISAT2 FC output not found"
fi

if [ -f "Quantification/Hisat2/HT/HTSeq_Count_union.txt" ]; then
    cp Quantification/Hisat2/HT/HTSeq_Count_union.txt DEG/H_HTSeq_Count_union.txt
    echo "  ✓ Copied HISAT2 HTSeq"
else
    echo "  ⚠ HISAT2 HTSeq output not found"
fi

if [ -f "Quantification/Hisat2/RSEM/rsem_gene_counts_matrix.txt" ]; then
    cp Quantification/Hisat2/RSEM/rsem_gene_counts_matrix.txt DEG/H_rsem_gene_counts_matrix.txt
    echo "  ✓ Copied HISAT2 RSEM"
else
    echo "  ⚠ HISAT2 RSEM output not found"
fi

# Bowtie2 outputs (FC + HTSeq + RSEM)
if [ -f "Quantification/Bowtie/FC/FC_Count.txt" ]; then
    cp Quantification/Bowtie/FC/FC_Count.txt DEG/B_FC_Count.txt
    echo "  ✓ Copied Bowtie2 featureCounts"
else
    echo "  ⚠ Bowtie2 FC output not found"
fi

if [ -f "Quantification/Bowtie/HT/HTSeq_Count_union.txt" ]; then
    cp Quantification/Bowtie/HT/HTSeq_Count_union.txt DEG/B_HTSeq_Count_union.txt
    echo "  ✓ Copied Bowtie2 HTSeq"
else
    echo "  ⚠ Bowtie2 HTSeq output not found"
fi

if [ -f "Quantification/Bowtie/RSEM/rsem_gene_counts_matrix.txt" ]; then
    cp Quantification/Bowtie/RSEM/rsem_gene_counts_matrix.txt DEG/B_rsem_gene_counts_matrix.txt
    echo "  ✓ Copied Bowtie2 RSEM"
else
    echo "  ⚠ Bowtie2 RSEM output not found"
fi

echo ""
echo "  ✓ All available quantification outputs copied"

# ===================================================================
# === Step 3: Clean Count Matrices (Unified Script) ===
# ===================================================================
echo ""
echo "=========================================="
echo "STEP 3: Cleaning Count Matrices"
echo "=========================================="

echo ">>> Running unified cleaning script..."
micromamba run -n rnaseq-cli \
    --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
    --env TRIM_DIR --env GTF --env FASTA --env GFF \
    /opt/project/scripts/bowtie_clean.sh
micromamba run -n rnaseq-cli \
    --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
    --env TRIM_DIR --env GTF --env FASTA --env GFF \
    /opt/project/scripts/hisat2_clean.sh
micromamba run -n rnaseq-cli \
    --env THREADS --env MODE --env GENOME_DIR --env READ_DIR \
    --env TRIM_DIR --env GTF --env FASTA --env GFF \
    /opt/project/scripts/star_clean.sh

if [ $? -ne 0 ]; then
    echo "ERROR: Cleaning failed!" >&2
    exit 1
fi

echo ""
echo "  ✓ All count matrices cleaned successfully"

# === Step 4: DEG Analysis ===
echo ""
echo "=========================================="
echo "STEP 4: Differential Expression Analysis"
echo "=========================================="

echo ">>> Running edgeR DEG Analysis..."
micromamba run -n rnaseq-r \
    --env SAMPLE_INFO_PATH \
    /opt/project/scripts/run_all_edgeR.sh > ./logs/edgeR_${TIMESTAMP}.log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: edgeR analysis failed!" >&2
    echo "Check ./logs/edgeR_${TIMESTAMP}.log for details" >&2
    exit 1
fi
echo "  ✓ edgeR analysis completed"

echo ""
echo ">>> Running DESeq2 DEG Analysis..."
micromamba run -n rnaseq-r \
    --env SAMPLE_INFO_PATH \
    /opt/project/scripts/run_all_deseq2.sh > ./logs/deseq2_${TIMESTAMP}.log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: DESeq2 analysis failed!" >&2
    echo "Check ./logs/deseq2_${TIMESTAMP}.log for details" >&2
    exit 1
fi
echo "  ✓ DESeq2 analysis completed"

# === Step 5: Final Verification ===
echo ""
echo "=========================================="
echo "STEP 5: Final Verification"
echo "=========================================="

DEG_RESULTS=$(find ./DEG -name "*_UP.xlsx" -o -name "*_DOWN.xlsx" | wc -l)
SUMMARY_FILES=$(find ./DEG -name "*_Summary.xlsx" | wc -l)
PCA_FILES=$(find ./DEG -name "*_PCA*.xlsx" -o -name "*_PCA*.png" | wc -l)
VOLCANO_FILES=$(find ./DEG -name "*_Volcano.png" | wc -l)

echo ">>> DEG Analysis Results Summary:"
echo "  • DEG files (UP/DOWN): $DEG_RESULTS"
echo "  • Summary files: $SUMMARY_FILES"
echo "  • PCA files: $PCA_FILES"
echo "  • Volcano plots: $VOLCANO_FILES"

if [ $DEG_RESULTS -eq 0 ]; then
    echo ""
    echo "WARNING: No DEG result files were generated!" >&2
    echo "Check logs in ./logs/ for R script errors" >&2
else
    echo ""
    echo "  ✓ DEG results successfully generated"
fi

# === Finalize ===
echo ""
echo ">>> Finalizing permissions..."
chmod -R 777 ./DEG ./Indices ./Mapping ./Quantification ./qc ./Trim ./logs 2>/dev/null || true

echo ""
echo "=========================================="
echo "RNA-Seq Pipeline Completed: $(date)"
echo "=========================================="
echo ""
echo "✅ Summary:"
echo "  ✓ Quality Control: Complete"
echo "  ✓ Alignment (HISAT2, Bowtie2, STAR): Complete"
echo "  ✓ Quantification (featureCounts + HTSeq + RSEM): Complete"
echo "  ✓ Count Matrix Cleaning: Complete"
echo "  ✓ DEG Analysis (edgeR, DESeq2): Complete"
echo ""
echo "📊 Pipeline Combinations Generated:"
echo "  • Bowtie2-FC-edgeR, Bowtie2-FC-DESeq2"
echo "  • Bowtie2-HT-edgeR, Bowtie2-HT-DESeq2"
echo "  • Bowtie2-RSEM-edgeR, Bowtie2-RSEM-DESeq2"
echo "  • HISAT2-FC-edgeR, HISAT2-FC-DESeq2"
echo "  • STAR-FC-edgeR, STAR-FC-DESeq2"
echo "  • STAR-HT-edgeR, STAR-HT-DESeq2"
echo "  • STAR-RSEM-edgeR, STAR-RSEM-DESeq2"
echo "  TOTAL: 16 combinations"
echo ""
echo "  • Total DEG result files: $DEG_RESULTS"
echo "  • Total Volcano plots: $VOLCANO_FILES"
echo ""
echo "📁 Results Location:"
echo "  • Results: /data/DEG/"
echo "  • Logs: /data/logs/"
echo ""
echo "✅ Pipeline execution completed successfully!"
echo ""

set +x
