#!/bin/bash

# Improved script to clean RNA-seq count matrices
echo "Cleaning RNA-seq count matricex"
# Clean featureCounts file
awk -F'\t' '
BEGIN { OFS = "," }
NR == 2 {
    printf "Geneid"
    for (i = 7; i <= NF; i++) {
        gsub(/^.*\//, "", $i)    # remove path
        gsub(/\.bam$/, "", $i)   # remove .bam extension
        printf ",%s", $i
    }
    printf "\n"
    next
}
NR > 2 {
    printf "%s", $1
    for (i = 7; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
}
' H_FC_Count.txt > H_FC_counts_clean.txt
rm H_FC_Count.txt

echo "featureCounts file cleaned successfully: H_FC_counts_clean.csv"
#######
# Clean HTSeq union file
# Clean HTSeq union file
echo "Processing HTSeq file..."
awk -F'\t' '
BEGIN { OFS = "," }
NR == 1 {
    printf "Geneid"
    for (i = 2; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
    next
}
NR > 1 && $1 !~ /^__/ {
    printf "%s", $1
    for (i = 2; i <= NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
}
' H_HTSeq_Count_union.txt > H_HTSeq_counts_clean.txt
rm H_HTSeq_Count_union.txt

#mv  B_rsem_gene_counts_matrix.txt B_RSEM_counts_clean.txt
echo "Cleaning complete!"
echo "Output files created:"
