#!/bin/bash
echo "Cleaning RNA-seq count matrices..."
# Clean HTSeq union file
#!/bin/bash

echo "🧼 Cleaning HTSeq count matrix and converting to comma-separated format..."

# Input and output paths
INPUT="S_HTSeq_Count_union.txt"
OUTPUT="S_HTSeq_Count_union_clean.txt"

# Create clean CSV-style header
HEADER=$(head -n1 "$INPUT" | awk '{
    printf "Geneid"
    for (i=2; i<=NF; i++) {
        gsub(/_Aligned.*$/, "", $i)
        printf ",%s", $i
    }
    printf "\n"
}')

# Write cleaned header to output
echo "$HEADER" > "$OUTPUT"

# Process and append cleaned data (excluding __rows)
awk '!/^__/ {
    printf $1
    for (i=2; i<=NF; i++) {
        printf ",%s", $i
    }
    printf "\n"
}' "$INPUT" | tail -n +2 >> "$OUTPUT"

echo "✅ Final cleaned CSV-style matrix saved to: $OUTPUT"
rm S_HTSeq_Count_union.txt
#######
awk -F'\t' '
BEGIN { OFS = "," }
NR == 2 {
    printf "Geneid"
    for (i = 7; i <= NF; i++) {
        name = $i
        gsub(/^.*\//, "", name)  # Remove path
        gsub(/_Aligned.*/, "", name)
        printf ",%s", name
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
' S_FC_Count.txt > S_FC_counts_clean.txt
rm S_FC_Count.txt
# Clean RSEM gene counts matrix
sed 's/\t/,/g' S_rsem_gene_counts_matrix.txt > S_RSEM_counts_clean.txt
rm S_rsem_gene_counts_matrix.txt 
echo "Cleaning complete!"
echo "Output files created:"

# Quick validation
echo -e "\nValidation check:"
