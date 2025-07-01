#!/bin/bash

# Improved script to clean RNA-seq count matrices
echo "Cleaning RNA-seq count matrices..."

# Improved script to clean RNA-seq count matrices
echo "Cleaning RNA-seq count matrices..."

# Clean featureCounts file
echo "Processing featureCounts file..."
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
' B_FC_Count.txt > B_FC_counts_clean.txt
rm B_FC_Count.txt
echo "featureCounts file cleaned successfully: B_FC_counts_clean.txt"
######
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
' B_HTSeq_Count_union.txt > B_HTSeq_counts_clean.txt
rm B_HTSeq_Count_union.txt
# Clean RSEM gene counts matrix
sed 's/\t/,/g' B_rsem_gene_counts_matrix.txt > B_RSEM_counts_clean.txt
rm B_rsem_gene_counts_matrix.txt
echo "Cleaning complete!"
echo "Output files created:"
echo " - B_FC_counts_clean.txt ($(wc -l < B_FC_counts_clean.txt) lines)"
echo " - B_HTSeq_counts_clean.txt ($(wc -l < B_HTSeq_counts_clean.txt) lines)"
echo " - B_RSEM_counts_clean.txt ($(wc -l < B_RSEM_counts_clean.txt) lines)"

# Quick validation
echo -e "\nValidation check:"
echo "featureCounts header:"
head -n 1 B_FC_counts_clean.txt
echo "HTSeq header:"
head -n 1 B_HTSeq_counts_clean.txt
echo "RSEM header:"
head -n 1 B_RSEM_counts_clean.txt
