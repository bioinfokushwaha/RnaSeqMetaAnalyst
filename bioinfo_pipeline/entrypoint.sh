#!/bin/bash

# Exit on any error
set -e

# Name of your conda environment
ENV_NAME="rnaseq"

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

# Echo environment info for debug
echo "🧬 Conda environment '$ENV_NAME' activated."
echo "📁 Working directory: $(pwd)"
echo "🧾 Running command: $@"

# If no command is given, start a shell
if [ $# -eq 0 ]; then
    exec bash
else
    exec "$@"
fi

