#!/bin/bash
# ===================================================================
# WRAPPER: Run DEG Analysis Only
# Uses rnaseq-r conda environment
# ===================================================================

set -e

echo "=========================================="
echo "DEG ANALYSIS ONLY"
echo "Environment: rnaseq-r"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Check if we're inside container
if [ ! -f "/opt/project/scripts/master.sh" ]; then
    echo "ERROR: This script should be run inside the container" >&2
    echo ""
    echo "Usage from host:" >&2
    echo "  docker run --rm \\" >&2
    echo "    -v \$(pwd)/data:/data/data:ro \\" >&2
    echo "    -v \$(pwd)/results:/data \\" >&2
    echo "    rnaseqmetaanalyst:latest \\" >&2
    echo "    /opt/project/scripts/run_deg_only.sh" >&2
    echo ""
    exit 1
fi

# Call master.sh with --deg-only flag
exec /opt/project/scripts/master.sh --deg-only "$@"
