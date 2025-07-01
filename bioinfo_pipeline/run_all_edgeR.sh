#!/bin/bash

RSCRIPT="run_edgeR_analysis.R"
SAMPLEINFO="DEG/sampleinfo.txt"

for count_file in DEG/*.txt; do
    if [[ "$count_file" == *"sampleinfo.txt" ]]; then
        continue
    fi

    echo "Running edgeR for $count_file"
    Rscript "$RSCRIPT" "$count_file" "$SAMPLEINFO"
    echo "Finished $count_file"
    echo "-----------------------------"
done
