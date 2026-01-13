#!/bin/bash
# ===================================================================
# WRAPPER: Run Comparison Analysis Only
# Uses rnaseq-comparison conda environment (Python)
# ===================================================================

set -e

echo "=========================================="
echo "COMPARISON ANALYSIS ONLY"
echo "Environment: rnaseq-comparison (Python)"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Check if we're inside container
if [ ! -f "/opt/project/scripts/master.sh" ]; then
    echo "ERROR: This script should be run inside the container" >&2
    echo ""
    echo "Usage from host:" >&2
    echo "  docker run --rm \\" >&2
    echo "    -v \$(pwd)/results:/data \\" >&2
    echo "    rnaseqmetaanalyst:latest \\" >&2
    echo "    /opt/project/scripts/run_comparison_only.sh" >&2
    echo ""
    exit 1
fi

# Call master.sh with --comparison-only flag
exec /opt/project/scripts/master.sh --comparison-only "$@"
