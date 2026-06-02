#!/bin/bash
# ===================================================================
# RNA-Seq Meta-Analyst Pipeline Launcher - OPTIMIZED VERSION
# ===================================================================
# Modes:
#   (default)          Full RNA-seq pipeline  → master.sh
#   --busco            BUSCO assessment        → master_busco.sh
#
# BUSCO can be run standalone via run_busco.sh or through this script
# with --busco plus its required arguments.  Both paths produce
# identical output; use whichever is more convenient.
# ===================================================================
set -e

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================
function show_help {
    cat << 'EOF'
=== RNA-Seq Meta-Analyst Pipeline Launcher - OPTIMIZED ===

Usage: ./run_pipeline1.sh [OPTIONS]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 RNA-SEQ PIPELINE MODE  (default)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Required Arguments:
  --project-dir    <path>  Absolute path to the project directory.
  --mode           <SE|PE> Sequencing mode: Single-End (SE) or Paired-End (PE).

Path Arguments (auto-detected from --project-dir if standard structure):
  --reads-dir      <path>  Path to raw FASTQ files (default: <project-dir>/data/raw)
  --genome-dir     <path>  Path to genome files (default: <project-dir>/data/genome)
  --sample-info    <path>  Path to sampleinfo.txt (default: <project-dir>/data/sampleinfo.txt)
  --output-dir     <path>  Results directory (default: <project-dir>/results)

Pre-built Index Support:
  --index-dir      <path>  Path to pre-built indices directory
                           Expected structure:
                             <path>/Bowtie/genome.*.bt2
                             <path>/Hisat2/genome.*.ht2
                             <path>/STAR/SA, Genome, etc.
                             <path>/RSEM/Bowtie/rsem_ref.*
                             <path>/RSEM/STAR/rsem_ref.*

Pipeline Control:
  --skip-alignment      Skip alignment and quantification
  --skip-deg            Skip DEG analysis
  --skip-comparison     Skip comparison analysis
  --deg-only            Only run DEG (requires existing count matrices)
  --comparison-only     Only run comparison (requires existing DEG results)

DEG Thresholds:
  --log2fc <value>      log2 Fold Change threshold (default: 0)
  --fdr <value>         FDR threshold (default: 0.05)
  --pvalue <value>      p-value threshold (default: 0.05)

Alignment Options:
  --sjdb-overhang <int> STAR SJDB overhang (default: 100)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 BUSCO MODE  (--busco)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  --busco               Run BUSCO assessment + GFF/GTF filtering
                        instead of the RNA-seq pipeline.

Required with --busco:
  --protein      <path> Path to protein FASTA file (protein.faa)
                        Relative to --project-dir or absolute
  --lineage      <name> BUSCO lineage dataset
                        (e.g. mammalia_odb12, bacteria_odb10)

Optional with --busco:
  --genome-dir   <path> Directory containing GFF/GTF files
                        (default: <project-dir>/data/genome)
  --gff          <path> Explicit path to GFF3 file
  --gtf          <path> Explicit path to GTF file
  --output-dir   <path> Results directory (default: <project-dir>/results)

  Outputs written to: <output-dir>/busco_output/
  housekeeping_genes.xlsx auto-copied to <output-dir>/ for
  downstream DEG analysis.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 SHARED OPTIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  --threads      <int>  Number of threads (default: 8)
  --image-name   <name> Docker image name (default: rnaseqmetaanalyst:latest)
  --detached            Run container in detached/background mode
  --help                Show this help message

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Full RNA-seq pipeline with pre-built indices
  ./run_pipeline1.sh --project-dir $(pwd) --mode PE \
    --index-dir /path/to/shared/indices

  # DEG only (after alignment is done)
  ./run_pipeline1.sh --project-dir $(pwd) --mode PE --deg-only

  # BUSCO — auto-detects GFF/GTF from genome-dir
  ./run_pipeline1.sh --project-dir $(pwd) \
    --busco \
    --protein  data/genome/protein.faa \
    --lineage  mammalia_odb12

  # BUSCO — explicit paths, background, 16 threads
  ./run_pipeline1.sh --project-dir $(pwd) \
    --busco \
    --protein  data/genome/protein.faa \
    --lineage  mammalia_odb12 \
    --gff      data/genome/genomic.gff \
    --gtf      data/genome/genomic.gtf \
    --threads  16 \
    --detached

  # Recommended full workflow:
  #   Step 1 — BUSCO (generates housekeeping_genes.xlsx)
  ./run_pipeline1.sh --project-dir $(pwd) \
    --busco --protein data/genome/protein.faa --lineage mammalia_odb12
  #   Step 2 — RNA-seq pipeline (picks up housekeeping_genes.xlsx automatically)
  ./run_pipeline1.sh --project-dir $(pwd) --mode PE

Available BUSCO lineages (common):
  mammalia_odb12     vertebrata_odb10    eukaryota_odb10
  bacteria_odb10     fungi_odb10         insecta_odb10
  embryophyta_odb10  actinobacteria_odb10  firmicutes_odb10
  Full list: https://busco.ezlab.org/list_of_lineages.html
EOF
}

# ===================================================================
# DEFAULT VALUES
# ===================================================================
THREADS=8
IMAGE_NAME="rnaseqmetaanalyst:latest"
DETACHED_MODE=false
LOG2FC_THRESHOLD=0
FDR_THRESHOLD=0.05
PVALUE_THRESHOLD=0.05
STAR_SJDB_OVERHANG=100
PIPELINE_FLAGS=""
INDEX_DIR=""

# BUSCO-specific defaults
RUN_BUSCO=false
PROTEIN_FILE=""
BUSCO_LINEAGE=""
HOST_GFF=""
HOST_GTF=""

# ===================================================================
# PARSE COMMAND-LINE ARGUMENTS
# ===================================================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        # ── Shared ────────────────────────────────────────────────
        --project-dir)   PROJECT_DIR="$2";      shift ;;
        --genome-dir)    HOST_GENOME_DIR="$2";  shift ;;
        --output-dir)    HOST_OUTPUT_DIR="$2";  shift ;;
        --threads)       THREADS="$2";          shift ;;
        --image-name)    IMAGE_NAME="$2";       shift ;;
        --detached)      DETACHED_MODE=true ;;
        --help)          show_help; exit 0 ;;

        # ── RNA-seq pipeline ──────────────────────────────────────
        --reads-dir)      HOST_READS_DIR="$2";    shift ;;
        --sample-info)    HOST_SAMPLE_INFO="$2";  shift ;;
        --index-dir)      HOST_INDEX_DIR="$2";    shift ;;
        --mode)           MODE="$2";              shift ;;
        --skip-alignment) PIPELINE_FLAGS+=" --skip-alignment" ;;
        --skip-deg)       PIPELINE_FLAGS+=" --skip-deg" ;;
        --skip-comparison)PIPELINE_FLAGS+=" --skip-comparison" ;;
        --deg-only)       PIPELINE_FLAGS+=" --deg-only" ;;
        --comparison-only)PIPELINE_FLAGS+=" --comparison-only" ;;
        --log2fc)
            LOG2FC_THRESHOLD="$2"
            PIPELINE_FLAGS+=" --log2fc $2"
            shift ;;
        --fdr)
            FDR_THRESHOLD="$2"
            PIPELINE_FLAGS+=" --fdr $2"
            shift ;;
        --pvalue)
            PVALUE_THRESHOLD="$2"
            PIPELINE_FLAGS+=" --pvalue $2"
            shift ;;
        --sjdb-overhang)
            STAR_SJDB_OVERHANG="$2"
            PIPELINE_FLAGS+=" --sjdb-overhang $2"
            shift ;;

        # ── BUSCO mode ────────────────────────────────────────────
        --busco)    RUN_BUSCO=true ;;
        --protein)  PROTEIN_FILE="$2";  shift ;;
        --lineage)  BUSCO_LINEAGE="$2"; shift ;;
        --gff)      HOST_GFF="$2";      shift ;;
        --gtf)      HOST_GTF="$2";      shift ;;

        *) echo "ERROR: Unknown parameter: $1"; echo ""; show_help; exit 1 ;;
    esac
    shift
done

# ===================================================================
# VALIDATE SHARED REQUIRED ARGS
# ===================================================================
if [ -z "$PROJECT_DIR" ]; then
    echo "ERROR: --project-dir is required."
    echo ""
    show_help
    exit 1
fi

# Make PROJECT_DIR absolute
PROJECT_DIR=$(realpath "$PROJECT_DIR")

# ===================================================================
# ██████╗ ██╗   ██╗███████╗ ██████╗ ██████╗     ███╗   ███╗ ██████╗ ██████╗ ███████╗
# ██╔══██╗██║   ██║██╔════╝██╔════╝██╔═══██╗    ████╗ ████║██╔═══██╗██╔══██╗██╔════╝
# ██████╔╝██║   ██║███████╗██║     ██║   ██║    ██╔████╔██║██║   ██║██║  ██║█████╗
# ██╔══██╗██║   ██║╚════██║██║     ██║   ██║    ██║╚██╔╝██║██║   ██║██║  ██║██╔══╝
# ██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝    ██║ ╚═╝ ██║╚██████╔╝██████╔╝███████╗
# ╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝     ╚═╝     ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
# ===================================================================
if [ "$RUN_BUSCO" = true ]; then

    # ------------------------------------------------------------------
    # If run_busco.sh lives next to this script, delegate to it — this
    # avoids duplicating logic and keeps BUSCO launch behaviour in sync.
    # ------------------------------------------------------------------
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    BUSCO_LAUNCHER="${SCRIPT_DIR}/run_busco.sh"

    if [ -f "$BUSCO_LAUNCHER" ] && [ -x "$BUSCO_LAUNCHER" ]; then
        echo ""
        echo "=========================================="
        echo "BUSCO MODE  →  delegating to run_busco.sh"
        echo "=========================================="
        echo ""

        # Build the argument list for run_busco.sh
        BUSCO_ARGS=(
            --project-dir "$PROJECT_DIR"
            --threads     "$THREADS"
            --image-name  "$IMAGE_NAME"
        )

        [ -n "$PROTEIN_FILE"    ] && BUSCO_ARGS+=(--protein    "$PROTEIN_FILE")
        [ -n "$BUSCO_LINEAGE"   ] && BUSCO_ARGS+=(--lineage    "$BUSCO_LINEAGE")
        [ -n "$HOST_GENOME_DIR" ] && BUSCO_ARGS+=(--genome-dir "$HOST_GENOME_DIR")
        [ -n "$HOST_GFF"        ] && BUSCO_ARGS+=(--gff        "$HOST_GFF")
        [ -n "$HOST_GTF"        ] && BUSCO_ARGS+=(--gtf        "$HOST_GTF")
        [ -n "$HOST_OUTPUT_DIR" ] && BUSCO_ARGS+=(--output-dir "$HOST_OUTPUT_DIR")
        [ "$DETACHED_MODE" = true ] && BUSCO_ARGS+=(--detached)

        exec "$BUSCO_LAUNCHER" "${BUSCO_ARGS[@]}"
        # exec replaces this process — nothing below runs
    fi

    # ------------------------------------------------------------------
    # run_busco.sh not found — run BUSCO logic inline
    # (identical behaviour, no external dependency)
    # ------------------------------------------------------------------
    echo ""
    echo "=========================================="
    echo "BUSCO PIPELINE LAUNCHER  (inline mode)"
    echo "=========================================="
    echo ""

    # ── Validate BUSCO-required args ──────────────────────────────────
    BUSCO_MISSING=0
    if [ -z "$PROTEIN_FILE" ]; then
        echo "ERROR: --protein is required in --busco mode" >&2
        BUSCO_MISSING=1
    fi
    if [ -z "$BUSCO_LINEAGE" ]; then
        echo "ERROR: --lineage is required in --busco mode" >&2
        BUSCO_MISSING=1
    fi
    [ $BUSCO_MISSING -ne 0 ] && { echo ""; show_help; exit 1; }

    # ── Resolve paths ─────────────────────────────────────────────────
    # Protein
    if [[ "$PROTEIN_FILE" != /* ]]; then
        PROTEIN_FILE="${PROJECT_DIR}/${PROTEIN_FILE}"
    fi

    # Genome dir
    if [ -z "$HOST_GENOME_DIR" ]; then
        if [ -d "${PROJECT_DIR}/data/genome" ]; then
            HOST_GENOME_DIR="${PROJECT_DIR}/data/genome"
        else
            HOST_GENOME_DIR="${PROJECT_DIR}/data"
        fi
    fi

    # GFF
    if [ -z "$HOST_GFF" ]; then
        HOST_GFF="${HOST_GENOME_DIR}/genomic.gff"
        [ ! -f "$HOST_GFF" ] && HOST_GFF="${HOST_GENOME_DIR}/genomic.gff3"
    fi

    # GTF
    if [ -z "$HOST_GTF" ]; then
        HOST_GTF="${HOST_GENOME_DIR}/genomic.gtf"
    fi

    # Output dir
    HOST_OUTPUT_DIR="${HOST_OUTPUT_DIR:-${PROJECT_DIR}/results}"

    # ── Pre-run validation ────────────────────────────────────────────
    echo ">>> Performing pre-run checks..."

    if ! command -v docker &> /dev/null; then
        echo "ERROR: Docker not found." >&2; exit 1
    fi
    echo "  [OK] Docker available"

    if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
        echo "ERROR: Docker image '$IMAGE_NAME' not found." >&2
        echo "       Build it: docker build -t $IMAGE_NAME ." >&2
        exit 1
    fi
    echo "  [OK] Docker image: $IMAGE_NAME"

    [ ! -f "$PROTEIN_FILE" ] && { echo "ERROR: Protein file not found: $PROTEIN_FILE" >&2; exit 1; }
    echo "  [OK] Protein file: $PROTEIN_FILE"

    [ ! -f "$HOST_GFF" ] && {
        echo "ERROR: GFF file not found: $HOST_GFF" >&2
        echo "       Use --gff to specify path explicitly." >&2
        exit 1
    }
    echo "  [OK] GFF file: $HOST_GFF"

    [ ! -f "$HOST_GTF" ] && {
        echo "ERROR: GTF file not found: $HOST_GTF" >&2
        echo "       Use --gtf to specify path explicitly." >&2
        exit 1
    }
    echo "  [OK] GTF file: $HOST_GTF"
    echo "  [OK] Lineage:  $BUSCO_LINEAGE"
    echo "  [OK] Threads:  $THREADS"

    # ── Prepare output dirs ───────────────────────────────────────────
    echo ""
    echo ">>> Preparing output directories..."
    mkdir -p "${HOST_OUTPUT_DIR}/busco_output"
    mkdir -p "${HOST_OUTPUT_DIR}/logs"
    chmod -R 777 "$HOST_OUTPUT_DIR" 2>/dev/null || true
    echo "  [OK] Output directory: $HOST_OUTPUT_DIR"

    # ── Resolve container-internal paths ─────────────────────────────
    PROTEIN_FILENAME=$(basename "$PROTEIN_FILE")
    PROTEIN_HOST_DIR=$(dirname "$PROTEIN_FILE")
    GFF_FILENAME=$(basename "$HOST_GFF")
    GTF_FILENAME=$(basename "$HOST_GTF")

    CONTAINER_GFF="/data/genome/${GFF_FILENAME}"
    CONTAINER_GTF="/data/genome/${GTF_FILENAME}"

    if [ "$(realpath "$PROTEIN_HOST_DIR")" = "$(realpath "$HOST_GENOME_DIR")" ]; then
        CONTAINER_PROTEIN="/data/genome/${PROTEIN_FILENAME}"
        EXTRA_PROTEIN_MOUNT=""
    else
        CONTAINER_PROTEIN="/data/protein/${PROTEIN_FILENAME}"
        EXTRA_PROTEIN_MOUNT="-v ${PROTEIN_HOST_DIR}:/data/protein:ro"
    fi

    # ── Build Docker args ─────────────────────────────────────────────
    echo ""
    echo ">>> Preparing Docker configuration..."
    HOST_UID=$(id -u)
    HOST_GID=$(id -g)
    echo "  [INFO] Running as UID:GID = ${HOST_UID}:${HOST_GID}"

    DOCKER_ARGS="--rm"
    DOCKER_ARGS+=" --user ${HOST_UID}:${HOST_GID}"
    [ "$DETACHED_MODE" = true ] && DOCKER_ARGS+=" -d"
    DOCKER_ARGS+=" -v ${HOST_GENOME_DIR}:/data/genome:ro"
    DOCKER_ARGS+=" -v ${HOST_OUTPUT_DIR}:/data"
    [ -n "$EXTRA_PROTEIN_MOUNT" ] && DOCKER_ARGS+=" $EXTRA_PROTEIN_MOUNT"
    DOCKER_ARGS+=" -e PROTEIN_FILE=${CONTAINER_PROTEIN}"
    DOCKER_ARGS+=" -e BUSCO_LINEAGE=${BUSCO_LINEAGE}"
    DOCKER_ARGS+=" -e GFF_FILE=${CONTAINER_GFF}"
    DOCKER_ARGS+=" -e GTF_FILE=${CONTAINER_GTF}"
    DOCKER_ARGS+=" -e THREADS=${THREADS}"

    # ── Config summary ────────────────────────────────────────────────
    echo ""
    echo "Configuration Summary:"
    echo "──────────────────────────────────────────"
    echo "  Image        : $IMAGE_NAME"
    echo "  Lineage      : $BUSCO_LINEAGE"
    echo "  Threads      : $THREADS"
    echo "  Protein file : $PROTEIN_FILE"
    echo "  GFF file     : $HOST_GFF"
    echo "  GTF file     : $HOST_GTF"
    echo "  Output dir   : $HOST_OUTPUT_DIR"
    echo "  Running as   : UID=$HOST_UID, GID=$HOST_GID"
    echo "  Mode         : $([ "$DETACHED_MODE" = true ] && echo BACKGROUND || echo FOREGROUND)"
    echo "──────────────────────────────────────────"
    echo ""

    # ── Launch ───────────────────────────────────────────────────────
    echo ">>> Launching BUSCO container..."
    LOG_FILE="${HOST_OUTPUT_DIR}/busco_run.log"

    if [ "$DETACHED_MODE" = true ]; then
        CONTAINER_ID=$(docker run $DOCKER_ARGS "$IMAGE_NAME" \
            /bin/bash -c "/opt/project/scripts/master_busco.sh")

        echo "=========================================="
        echo "✅ BUSCO CONTAINER STARTED IN BACKGROUND"
        echo "=========================================="
        echo "  Container ID : $CONTAINER_ID"
        echo "  Monitor      : docker logs -f $CONTAINER_ID"
        echo "  Log file     : tail -f $LOG_FILE"
        echo "  Stop         : docker stop $CONTAINER_ID"
        echo ""
        echo "$CONTAINER_ID" > "${HOST_OUTPUT_DIR}/busco_container_id.txt"
        docker logs -f "$CONTAINER_ID" > "$LOG_FILE" 2>&1 &
    else
        docker run $DOCKER_ARGS "$IMAGE_NAME" \
            /bin/bash -c "/opt/project/scripts/master_busco.sh" 2>&1 | tee "$LOG_FILE"
        EXIT_CODE=${PIPESTATUS[0]}

        echo ""
        if [ $EXIT_CODE -eq 0 ]; then
            echo "=========================================="
            echo "✅ BUSCO PIPELINE COMPLETED SUCCESSFULLY"
            echo "=========================================="
            echo ""
            echo "Results:"
            echo "  BUSCO run dir  : ${HOST_OUTPUT_DIR}/busco_output/busco_run/"
            echo "  Core gene list : ${HOST_OUTPUT_DIR}/busco_output/busco_core_genes1.txt"
            echo "  Filtered GFF3  : ${HOST_OUTPUT_DIR}/busco_output/filtered/busco_filtered.gff"
            echo "  Filtered GTF   : ${HOST_OUTPUT_DIR}/busco_output/filtered/busco_filtered.gtf"
            echo "  Match table    : ${HOST_OUTPUT_DIR}/busco_output/filtered/busco_matched.tsv"
            echo "  HKG Excel      : ${HOST_OUTPUT_DIR}/housekeeping_genes.xlsx"
            echo "  Full log       : $LOG_FILE"
            echo ""
            echo "  ✓ housekeeping_genes.xlsx is ready for the RNA-seq pipeline."
            echo "    Run the pipeline next:"
            echo "      ./run_pipeline1.sh --project-dir $PROJECT_DIR --mode PE"
            echo ""
            # Quick BUSCO summary if available
            SUMMARY=$(find "${HOST_OUTPUT_DIR}/busco_output" -name "short_summary*.txt" 2>/dev/null | head -1)
            if [ -n "$SUMMARY" ]; then
                echo "BUSCO Summary:"
                echo "──────────────────────────────────────────"
                grep -A 10 "Results:" "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
                echo "──────────────────────────────────────────"
            fi
        else
            echo "=========================================="
            echo "❌ BUSCO PIPELINE FAILED"
            echo "=========================================="
            echo "  Exit code : $EXIT_CODE"
            echo "  Check log : $LOG_FILE"
            echo "  or logs/  : ${HOST_OUTPUT_DIR}/logs/"
        fi
        exit $EXIT_CODE
    fi

    # BUSCO mode ends here — nothing below should run
    exit 0
fi

# ===================================================================
# ██████╗ ███╗   ██╗ █████╗     ███████╗███████╗ ██████╗
# ██╔══██╗████╗  ██║██╔══██╗    ██╔════╝██╔════╝██╔═══██╗
# ██████╔╝██╔██╗ ██║███████║    ███████╗█████╗  ██║   ██║
# ██╔══██╗██║╚██╗██║██╔══██║    ╚════██║██╔══╝  ██║▄▄ ██║
# ██║  ██║██║ ╚████║██║  ██║    ███████║███████╗╚██████╔╝
# ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝   ╚══════╝╚══════╝ ╚══▀▀═╝
# ===================================================================

# ── Validate RNA-seq required args ────────────────────────────────
if [ -z "$MODE" ]; then
    echo "ERROR: --mode <SE|PE> is required for the RNA-seq pipeline."
    echo ""
    show_help
    exit 1
fi

MODE=$(echo "$MODE" | tr '[:lower:]' '[:upper:]')
if [[ "$MODE" != "SE" && "$MODE" != "PE" ]]; then
    echo "ERROR: --mode must be 'SE' or 'PE', got: $MODE"
    exit 1
fi

# ── Set default paths ─────────────────────────────────────────────
HOST_READS_DIR=${HOST_READS_DIR:-"${PROJECT_DIR}/data/raw"}
if [ -d "${PROJECT_DIR}/data/genome" ]; then
    HOST_GENOME_DIR=${HOST_GENOME_DIR:-"${PROJECT_DIR}/data/genome"}
else
    HOST_GENOME_DIR=${HOST_GENOME_DIR:-"${PROJECT_DIR}/data"}
fi
HOST_SAMPLE_INFO=${HOST_SAMPLE_INFO:-"${PROJECT_DIR}/data/sampleinfo.txt"}
HOST_OUTPUT_DIR=${HOST_OUTPUT_DIR:-"${PROJECT_DIR}/results"}

# ===================================================================
# PRE-RUN VALIDATION
# ===================================================================
echo ""
echo "=========================================="
echo "RNA-SEQ PIPELINE LAUNCHER - OPTIMIZED"
echo "=========================================="
echo ""
echo ">>> Performing pre-run checks..."

if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker not found. Is Docker installed and running?" >&2
    exit 1
fi

if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo "ERROR: Docker image '$IMAGE_NAME' not found" >&2
    echo "       Please build the image first:" >&2
    echo "       docker build -t $IMAGE_NAME ." >&2
    exit 1
fi

[ ! -d "$HOST_READS_DIR"  ] && { echo "ERROR: Reads directory not found: $HOST_READS_DIR" >&2;   exit 1; }
[ ! -d "$HOST_GENOME_DIR" ] && { echo "ERROR: Genome directory not found: $HOST_GENOME_DIR" >&2; exit 1; }
[ ! -f "$HOST_SAMPLE_INFO"] && { echo "ERROR: Sample info not found: $HOST_SAMPLE_INFO" >&2;     exit 1; }

echo "[OK] Docker is available"
echo "[OK] Docker image found: $IMAGE_NAME"
echo "[OK] Reads directory:    $HOST_READS_DIR"
echo "[OK] Genome directory:   $HOST_GENOME_DIR"
echo "[OK] Sample info:        $HOST_SAMPLE_INFO"

# ===================================================================
# VALIDATE INDEX DIRECTORY (if provided)
# ===================================================================
if [ -n "$HOST_INDEX_DIR" ]; then
    echo ""
    echo ">>> Validating pre-built index directory..."
    echo "    Path: $HOST_INDEX_DIR"

    [ ! -d "$HOST_INDEX_DIR" ] && { echo "ERROR: Index directory not found: $HOST_INDEX_DIR" >&2; exit 1; }

    HAS_INDICES=false
    [ -f "$HOST_INDEX_DIR/Bowtie/genome.1.bt2"          ] && { echo "[OK] Bowtie2 index found";     HAS_INDICES=true; }
    [ -f "$HOST_INDEX_DIR/Hisat2/genome.1.ht2"          ] && { echo "[OK] HISAT2 index found";      HAS_INDICES=true; }
    [ -f "$HOST_INDEX_DIR/STAR/SA"                       ] && { echo "[OK] STAR index found";        HAS_INDICES=true; }
    [ -f "$HOST_INDEX_DIR/RSEM/Bowtie/rsem_ref.idx.fa"  ] && { echo "[OK] RSEM-Bowtie index found"; HAS_INDICES=true; }
    [ -f "$HOST_INDEX_DIR/RSEM/STAR/rsem_ref.idx.fa"    ] && { echo "[OK] RSEM-STAR index found";   HAS_INDICES=true; }

    if [ "$HAS_INDICES" = false ]; then
        echo "WARNING: No indices found in $HOST_INDEX_DIR" >&2
        echo "         Continue anyway? Indices will be built locally. (y/N): "
        read -r CONTINUE
        [[ ! "$CONTINUE" =~ ^[Yy]$ ]] && exit 1
    else
        echo "[OK] Pre-built indices validated (missing ones will be built automatically)"
    fi
fi

# ===================================================================
# AUTO-DETECT GENOME FILES
# ===================================================================
echo ""
echo ">>> Auto-detecting genome files in: $HOST_GENOME_DIR"

GENOME_FASTA=""
for pattern in "*.fna" "*.fa" "*.fasta" "*.fna.gz" "*.fa.gz" "*.fasta.gz"; do
    GENOME_FASTA=$(find "$HOST_GENOME_DIR" -maxdepth 1 -type f -name "$pattern" 2>/dev/null | head -1)
    [ -n "$GENOME_FASTA" ] && { GENOME_FASTA=$(basename "$GENOME_FASTA"); break; }
done
[ -z "$GENOME_FASTA" ] && { echo "ERROR: Genome FASTA not found in $HOST_GENOME_DIR" >&2; exit 1; }
echo "[OK] Genome FASTA: $GENOME_FASTA"

GENOME_GTF=""
for pattern in "*.gtf" "*.gtf.gz" "genomic.gtf"; do
    GENOME_GTF=$(find "$HOST_GENOME_DIR" -maxdepth 1 -type f -name "$pattern" 2>/dev/null | head -1)
    [ -n "$GENOME_GTF" ] && { GENOME_GTF=$(basename "$GENOME_GTF"); break; }
done
[ -z "$GENOME_GTF" ] && { echo "ERROR: GTF file not found in $HOST_GENOME_DIR" >&2; exit 1; }
echo "[OK] GTF: $GENOME_GTF"

GENOME_GFF=""
for pattern in "*.gff" "*.gff3" "*.gff.gz" "*.gff3.gz" "genomic.gff"; do
    GENOME_GFF=$(find "$HOST_GENOME_DIR" -maxdepth 1 -type f -name "$pattern" 2>/dev/null | head -1)
    [ -n "$GENOME_GFF" ] && { GENOME_GFF=$(basename "$GENOME_GFF"); break; }
done
[ -z "$GENOME_GFF" ] && { echo "ERROR: GFF file not found in $HOST_GENOME_DIR" >&2; exit 1; }
echo "[OK] GFF: $GENOME_GFF"

FASTQ_COUNT=$(find "$HOST_READS_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null | wc -l)
[ "$FASTQ_COUNT" -eq 0 ] && { echo "ERROR: No FASTQ files found in $HOST_READS_DIR" >&2; exit 1; }
echo "[OK] Found $FASTQ_COUNT FASTQ files"

# ===================================================================
# CREATE OUTPUT DIRECTORIES
# ===================================================================
echo ""
echo ">>> Preparing output directories..."
mkdir -p "$HOST_OUTPUT_DIR"
mkdir -p "${HOST_OUTPUT_DIR}/results/DEG"
mkdir -p "${HOST_OUTPUT_DIR}/results/pipeline_comparison"
mkdir -p "${HOST_OUTPUT_DIR}/logs"
mkdir -p "${HOST_OUTPUT_DIR}/DEG"
chmod -R 777 "$HOST_OUTPUT_DIR" 2>/dev/null || true
echo "[OK] Results directory: $HOST_OUTPUT_DIR"

# Warn if housekeeping_genes.xlsx is missing — BUSCO should be run first
HKG_CANDIDATES=(
    "${HOST_OUTPUT_DIR}/housekeeping_genes.xlsx"
    "${HOST_OUTPUT_DIR}/busco_output/filtered/housekeeping_genes.xlsx"
    "${PROJECT_DIR}/data/housekeeping_genes.xlsx"
)
HKG_HOST_FOUND=""
for candidate in "${HKG_CANDIDATES[@]}"; do
    if [ -f "$candidate" ]; then
        HKG_HOST_FOUND="$candidate"
        break
    fi
done

if [ -n "$HKG_HOST_FOUND" ]; then
    echo "[OK] housekeeping_genes.xlsx found: $HKG_HOST_FOUND"
else
    echo ""
    echo "  ⚠  housekeeping_genes.xlsx not found."
    echo "     HKG analysis inside DESeq2/edgeR will be skipped."
    echo "     To generate it, run BUSCO first:"
    echo "       ./run_pipeline1.sh --project-dir $PROJECT_DIR \\"
    echo "         --busco --protein <protein.faa> --lineage <lineage>"
    echo "     (or use the standalone: ./run_busco.sh)"
    echo ""
fi

# ===================================================================
# PREPARE DOCKER VOLUMES AND ENVIRONMENT
# ===================================================================
echo ""
echo ">>> Preparing Docker configuration..."

HOST_UID=$(id -u)
HOST_GID=$(id -g)
echo "[INFO] Container will run as UID:GID = ${HOST_UID}:${HOST_GID}"

DOCKER_ARGS="--rm"
DOCKER_ARGS+=" --user ${HOST_UID}:${HOST_GID}"
[ "$DETACHED_MODE" = true ] && DOCKER_ARGS+=" -d"

DOCKER_ARGS+=" -v ${HOST_READS_DIR}:/data/reads:ro"
DOCKER_ARGS+=" -v ${HOST_GENOME_DIR}:/data/genome:ro"
DOCKER_ARGS+=" -v ${HOST_SAMPLE_INFO}:/data/sampleinfo.txt:ro"
DOCKER_ARGS+=" -v ${HOST_OUTPUT_DIR}:/data"

if [ -n "$HOST_INDEX_DIR" ]; then
    echo "[OK] Mounting pre-built indices: $HOST_INDEX_DIR"
    DOCKER_ARGS+=" -v ${HOST_INDEX_DIR}:/indices:ro"
    DOCKER_ARGS+=" -e INDEX_DIR=/indices"
    PIPELINE_FLAGS+=" --index-dir /indices"
fi

DOCKER_ARGS+=" -e THREADS=$THREADS"
DOCKER_ARGS+=" -e MODE=$MODE"
DOCKER_ARGS+=" -e READ_DIR=/data/reads"
DOCKER_ARGS+=" -e GENOME_DIR=/data/genome"
DOCKER_ARGS+=" -e TRIM_DIR=/data/Trim"
DOCKER_ARGS+=" -e GTF=/data/genome/$GENOME_GTF"
DOCKER_ARGS+=" -e FASTA=/data/genome/$GENOME_FASTA"
DOCKER_ARGS+=" -e GFF=/data/genome/$GENOME_GFF"
DOCKER_ARGS+=" -e LOG2FC_THRESHOLD=$LOG2FC_THRESHOLD"
DOCKER_ARGS+=" -e FDR_THRESHOLD=$FDR_THRESHOLD"
DOCKER_ARGS+=" -e PVALUE_THRESHOLD=$PVALUE_THRESHOLD"
DOCKER_ARGS+=" -e STAR_SJDB_OVERHANG=$STAR_SJDB_OVERHANG"

# ===================================================================
# CONFIGURATION SUMMARY
# ===================================================================
echo ""
echo "Configuration Summary:"
echo "──────────────────────────────────────────"
echo "  Image         : $IMAGE_NAME"
echo "  Mode          : $MODE"
echo "  Threads       : $THREADS"
echo "  Log2FC        : $LOG2FC_THRESHOLD"
echo "  FDR           : $FDR_THRESHOLD"
echo "  SJDB overhang : $STAR_SJDB_OVERHANG"
echo "  Genome FASTA  : $GENOME_FASTA"
echo "  GTF           : $GENOME_GTF"
echo "  GFF           : $GENOME_GFF"
echo "  Index source  : $([ -n "$HOST_INDEX_DIR" ] && echo "PRE-BUILT ($HOST_INDEX_DIR)" || echo "BUILD LOCALLY")"
echo "  HKG file      : $([ -n "$HKG_HOST_FOUND" ] && echo "$HKG_HOST_FOUND" || echo "not found (HKG analysis skipped)")"
echo "  Pipeline flags: $PIPELINE_FLAGS"
echo "  Running as    : UID=$HOST_UID, GID=$HOST_GID"
echo "──────────────────────────────────────────"
echo ""

# ===================================================================
# RUN THE CONTAINER
# ===================================================================
echo ">>> Launching container..."
echo ""

LOG_FILE="${HOST_OUTPUT_DIR}/pipeline_run.log"

if [ "$DETACHED_MODE" = true ]; then
    CONTAINER_ID=$(docker run $DOCKER_ARGS "$IMAGE_NAME" \
        /bin/bash -c "/opt/project/scripts/master.sh $PIPELINE_FLAGS")

    echo "=========================================="
    echo "✅ CONTAINER STARTED IN BACKGROUND"
    echo "=========================================="
    echo "  Container ID : $CONTAINER_ID"
    echo "  Monitor      : docker logs -f $CONTAINER_ID"
    echo "  Log file     : tail -f $LOG_FILE"
    echo "  Stop         : docker stop $CONTAINER_ID"
    echo ""
    echo "$CONTAINER_ID" > "${HOST_OUTPUT_DIR}/container_id.txt"
    docker logs -f "$CONTAINER_ID" > "$LOG_FILE" 2>&1 &
else
    docker run $DOCKER_ARGS "$IMAGE_NAME" \
        /bin/bash -c "/opt/project/scripts/master.sh $PIPELINE_FLAGS" 2>&1 | tee "$LOG_FILE"

    EXIT_CODE=${PIPESTATUS[0]}

    echo ""
    if [ $EXIT_CODE -eq 0 ]; then
        echo "=========================================="
        echo "✅ PIPELINE COMPLETED SUCCESSFULLY"
        echo "=========================================="
        echo ""
        echo "Results location: $HOST_OUTPUT_DIR"
        echo ""
        echo "DEG Results:"
        DEG_COUNT=$(find "${HOST_OUTPUT_DIR}/results/DEG" -name "*.xlsx" 2>/dev/null | wc -l)
        echo "  Found $DEG_COUNT Excel files"
        [ $DEG_COUNT -gt 0 ] && find "${HOST_OUTPUT_DIR}/results/DEG" -name "*.xlsx" 2>/dev/null | head -5

        echo ""
        echo "Performance:"
        [ -n "$HOST_INDEX_DIR" ] && echo "  ✓ Used pre-built indices" || echo "  ✓ Built indices locally"
        echo "  ✓ Parallel alignment execution"
        echo ""
        echo "Full log: $LOG_FILE"
    else
        echo "=========================================="
        echo "❌ PIPELINE FAILED"
        echo "=========================================="
        echo "  Exit code : $EXIT_CODE"
        echo "  Check log : $LOG_FILE"
        echo "  or logs/  : ${HOST_OUTPUT_DIR}/logs/"
    fi
    exit $EXIT_CODE
fi
