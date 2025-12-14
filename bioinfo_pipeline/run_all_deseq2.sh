#!/bin/bash
set -e

# Path to your R script
RSCRIPT="/opt/project/scripts/run_deseq2_analysis.R"
SAMPLEINFO="DEG/sampleinfo.txt"
# Validate that the sample info path is set
if [ -z "$SAMPLE_INFO_PATH" ] || [ ! -f "$SAMPLE_INFO_PATH" ]; then
    echo "ERROR: SAMPLE_INFO_PATH is not set or sampleinfo.txt not found at '$SAMPLE_INFO_PATH'" >&2
    exit 1
fi

# Loop through all cleaned CSV count files in the DEG directory
for count_file in DEG/*_clean.txt; do
    echo "Running DESeq2 for $count_file"
    Rscript "$RSCRIPT" "$count_file" "$SAMPLE_INFO_PATH"
    echo "Finished $count_file"
    echo "-----------------------------"
done
