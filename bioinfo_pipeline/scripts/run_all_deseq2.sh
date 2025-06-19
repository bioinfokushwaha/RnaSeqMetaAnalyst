#!/bin/bash

# Path to your R script
RSCRIPT="run_deseq2_analysis.R"

# Sample info file
SAMPLEINFO="DEG/sampleinfo.txt"

# Loop over all count files in DEG/
for count_file in DEG/*.txt; do
    # Skip sampleinfo.txt
    if [[ "$count_file" == *"sampleinfo.txt" ]]; then
        continue
    fi

    echo "Running DESeq2 for $count_file"
    Rscript "$RSCRIPT" "$count_file" "$SAMPLEINFO"
    echo "Finished $count_file"
    echo "-----------------------------"
done
