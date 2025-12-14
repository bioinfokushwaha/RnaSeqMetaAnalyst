#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

function show_help {
    echo "=== RNA-Seq Meta-Analyst Pipeline Launcher ==="
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Required Arguments:"
    echo "  --project-dir    <path>  Absolute path to the project directory."
    echo "  --mode           <SE|PE> Sequencing mode: Single-End (SE) or Paired-End (PE)."
    echo ""
    echo "Path Arguments (can be inferred from --project-dir if standard structure is used):"
    echo "  --reads-dir      <path>  Path to the directory with raw FASTQ files (default: <project-dir>/data/raw)."
    echo "  --genome-dir     <path>  Path to the directory with genome FASTA/GTF/GFF files (default: <project-dir>/data/genome)."
    echo "  --sample-info    <path>  Path to the sampleinfo.txt file (default: <project-dir>/data/sampleinfo.txt)."
    echo "  --output-dir     <path>  Path where the results will be saved (default: <project-dir>/results)."
    echo "  --index-dir      <path>  Optional: Path to a directory with pre-built global indices."
    echo ""
    echo "Other Options:"
    echo "  --threads        <int>   Number of threads to use (default: 8)."
    echo "  --image-name     <name>  Name of the Docker image to use (default: rnaseqmetaanalyst:latest)."
    echo "  --help                   Show this help message."
    echo ""
    echo "Example:"
    echo "  $0 --project-dir /home/user/my_rnaseq_project --mode SE --threads 16"
}

# --- Default values ---
THREADS=8
IMAGE_NAME="rnaseqmetaanalyst:latest"

# --- Parse Command-Line Arguments ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --project-dir) PROJECT_DIR="$2"; shift ;;
        --reads-dir) HOST_READS_DIR="$2"; shift ;;
        --genome-dir) HOST_GENOME_DIR="$2"; shift ;;
        --sample-info) HOST_SAMPLE_INFO="$2"; shift ;;
        --output-dir) HOST_OUTPUT_DIR="$2"; shift ;;
        --index-dir) HOST_INDEX_DIR="$2"; shift ;;
        --threads) THREADS="$2"; shift ;;
        --mode) MODE="$2"; shift ;;
        --image-name) IMAGE_NAME="$2"; shift ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done

# --- Validate Required Arguments & Set Defaults ---
if [ -z "$PROJECT_DIR" ] || [ -z "$MODE" ]; then
    echo "ERROR: --project-dir and --mode are required."
    show_help
    exit 1
fi

# Set default paths based on project directory if they were not provided
HOST_READS_DIR=${HOST_READS_DIR:-"${PROJECT_DIR}/data/raw"}
HOST_GENOME_DIR=${HOST_GENOME_DIR:-"${PROJECT_DIR}/data/genome"}
HOST_SAMPLE_INFO=${HOST_SAMPLE_INFO:-"${PROJECT_DIR}/data/sampleinfo.txt"}
HOST_OUTPUT_DIR=${HOST_OUTPUT_DIR:-"${PROJECT_DIR}/results"}

# Define filenames for inside the container (these are assumed to be consistent)
GTF="/data/genome/genomic.gtf"
FASTA="/data/genome/GCF_002263795.2_ARS-UCD1.3_genomic.fna"
GFF="/data/genome/genomic.gff"

# --- Pre-run Checks ---
echo ">>> Performing pre-run checks..."
if ! command -v docker &> /dev/null; then
    echo "ERROR: docker command could not be found. Is Docker installed and running?"
    exit 1
fi
if [ ! -d "$HOST_READS_DIR" ]; then
    echo "ERROR: Reads directory not found at $HOST_READS_DIR"
    exit 1
fi
if [ ! -d "$HOST_GENOME_DIR" ]; then
    echo "ERROR: Genome directory not found at $HOST_GENOME_DIR"
    exit 1
fi
if [ ! -f "$HOST_SAMPLE_INFO" ]; then
    echo "ERROR: Sample info file not found at $HOST_SAMPLE_INFO"
    exit 1
fi
echo ">>> All checks passed."

# Create the output directory if it doesn't exist
mkdir -p "$HOST_OUTPUT_DIR"
echo ">>> Results will be saved in: $HOST_OUTPUT_DIR"

# --- Prepare Docker Arguments ---
DOCKER_ARGS="--rm -d \
  -v ${HOST_READS_DIR}:/data/reads:ro \
  -v ${HOST_GENOME_DIR}:/data/genome:ro \
  -v ${HOST_SAMPLE_INFO}:/data/sampleinfo.txt:ro \
  -v ${HOST_OUTPUT_DIR}:/data"

# If a shared index directory is provided and exists, mount it and set the env var
if [ -n "$HOST_INDEX_DIR" ] && [ -d "$HOST_INDEX_DIR" ]; then
    echo ">>> Mounting pre-built indices from: $HOST_INDEX_DIR"
    DOCKER_ARGS+=" -v ${HOST_INDEX_DIR}:/indices:ro -e INDEX_DIR=/indices"
else
    echo ">>> Using local project indices (will be built if not found)."
fi

# --- Run the container in detached mode ---
echo ">>> Launching Docker container '$IMAGE_NAME' in the background..."

CONTAINER_ID=$(docker run $DOCKER_ARGS \
  -e THREADS="$THREADS" \
  -e MODE="$MODE" \
  -e READ_DIR=/data/reads \
  -e GENOME_DIR=/data/genome \
  -e TRIM_DIR=/data/Trim \
  -e GTF="$GTF" \
  -e FASTA="$FASTA" \
  -e GFF="$GFF" \
  "$IMAGE_NAME")

LOG_FILE="${HOST_OUTPUT_DIR}/pipeline_run.log"

echo "=== Docker container started with ID: $CONTAINER_ID ==="
echo "The pipeline will continue to run even if you close this terminal."
echo "Logs are being saved to: $LOG_FILE"
echo
echo "To monitor its progress, run:"
echo "  docker logs -f $CONTAINER_ID"
echo
echo "To stop the container, run:"
echo "  docker stop $CONTAINER_ID"
echo

# Start saving logs in the background
docker logs -f "$CONTAINER_ID" > "$LOG_FILE" 2>&1 &
