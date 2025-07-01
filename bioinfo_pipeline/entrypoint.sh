#!/bin/bash
set -e

ENV_NAME="16_ppl_bioinfo-env"
source /opt/conda/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

echo "🧬 Conda environment '$ENV_NAME' activated."
echo "📁 Working directory: $(pwd)"
echo "🧾 Running command: $@"

if [ $# -eq 0 ]; then
    exec bash
else
    exec "$@"
fi
