#!/bin/bash

THREADS=10
MODE="SE"
GENOME_DIR="data/genome"
READ_DIR="data/reads"
GTF="data/genome"
FASTA="data/genome"

echo "DEBUG: MODE is [$MODE]"
echo "GENOME_DIR is [$GENOME_DIR]"

mkdir -p Indices/Hisat2 Mapping/Hisat2 Quantification/Hisat2/{HT,FC,RSEM}

# Index genome (uncomment if not already done)
hisat2-build -p $THREADS $FASTA/*.fna Indices/Hisat2/Hisat2

echo "Genome indexing completed."

if [[ $MODE == "PE" ]]; then
    echo "Running in Paired-End mode..."

    for i in "$READ_DIR"/*_R1.fastq; do
        SAMPLE=$(basename "$i" _R1.fastq)
        R1="$READ_DIR/${SAMPLE}_R1.fastq"
        R2="$READ_DIR/${SAMPLE}_R2.fastq"

        echo "Processing sample: $SAMPLE"
        hisat2 -p $THREADS --dta -x Indices/Hisat2/Hisat2 \
            -1 "$R1" -2 "$R2" \
            -S Mapping/Hisat2/${SAMPLE}.sam \
            --summary-file Mapping/Hisat2/${SAMPLE}.summary

        echo "Converting SAM to BAM for $SAMPLE..."
        samtools view -@ $THREADS -b Mapping/Hisat2/${SAMPLE}.sam > Mapping/Hisat2/${SAMPLE}.bam
        rm Mapping/Hisat2/${SAMPLE}.sam

        echo "Sorting and indexing BAM for $SAMPLE..."
        samtools sort -@ $THREADS -o Mapping/Hisat2/${SAMPLE}_sorted.bam Mapping/Hisat2/${SAMPLE}.bam
        mv Mapping/Hisat2/${SAMPLE}_sorted.bam Mapping/Hisat2/${SAMPLE}.bam
        samtools index Mapping/Hisat2/${SAMPLE}.bam
    done

else
    echo "Running in Single-End mode..."

   for i in "$READ_DIR"/*.fastq; do
        SAMPLE=$(basename "$i" .fastq)
        R="$READ_DIR/${SAMPLE}.fastq"
        echo "Processing sample: $SAMPLE"
        hisat2 -p $THREADS --dta -x Indices/Hisat2/Hisat2 \
          -U "$R" \
            -S Mapping/Hisat2/${SAMPLE}.sam \
            --summary-file Mapping/Hisat2/${SAMPLE}.summary

        echo "Converting SAM to BAM for $SAMPLE..."
        samtools view -@ $THREADS -b Mapping/Hisat2/${SAMPLE}.sam > Mapping/Hisat2/${SAMPLE}.bam
        rm Mapping/Hisat2/${SAMPLE}.sam

        echo "Sorting and indexing BAM for $SAMPLE..."
        samtools sort -@ $THREADS -o Mapping/Hisat2/${SAMPLE}_sorted.bam Mapping/Hisat2/${SAMPLE}.bam
        mv Mapping/Hisat2/${SAMPLE}_sorted.bam Mapping/Hisat2/${SAMPLE}.bam
        samtools index Mapping/Hisat2/${SAMPLE}.bam
    done
fi
#############
##########################################
# Index BAM files before HTSeq
##########################################
echo "Indexing BAM files..."
for BAM in Mapping/Hisat2/*.bam; do
    if [ -f "$BAM" ] && [ ! -f "${BAM}.bai" ]; then
        echo "Indexing $BAM"
        samtools index "$BAM"
    fi
done
#########################
echo "Running htseq-count..."
echo "Running htseq-count sample-wise..."
HT_DIR="Quantification/Hisat2/"
mkdir -p "$HT_DIR"

for BAM in Mapping/Hisat2/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "Counting for sample: $SAMPLE"

    # Name sort the BAM
    samtools sort -n -@ $THREADS "$BAM" -o "$HT_DIR/${SAMPLE}_name_sorted.bam"

    # Convert to SAM
    samtools view -h "$HT_DIR/${SAMPLE}_name_sorted.bam" > "$HT_DIR/${SAMPLE}.sam"

    # Count
    htseq-count -f sam -r name -m union "$HT_DIR/${SAMPLE}.sam" "$GTF"/*.gtf > "$HT_DIR/${SAMPLE}.htseq.txt" 2> "$HT_DIR/${SAMPLE}.err"

    # Clean up
    rm "$HT_DIR/${SAMPLE}_name_sorted.bam" "$HT_DIR/${SAMPLE}.sam"
done

# Merge all HTSeq-count outputs
echo "Merging htseq-count files into a count matrix..."

cut -f1 "$HT_DIR/$(ls $HT_DIR/*.htseq.txt | head -n1)" > "$HT_DIR/gene_column.txt"

for FILE in "$HT_DIR"/*.htseq.txt; do
    SAMPLE=$(basename "$FILE" .htseq.txt)
    cut -f2 "$FILE" > "$HT_DIR/${SAMPLE}_counts.txt"
done

paste "$HT_DIR/gene_column.txt" "$HT_DIR/"*_counts.txt > "$HT_DIR/HTSeq_Count_union.txt"

HEADER="Gene"
for FILE in "$HT_DIR/"*_counts.txt; do
    SAMPLE=$(basename "$FILE" _counts.txt)
    HEADER+="\t$SAMPLE"
done
sed -i "1s/^/${HEADER}\n/" "$HT_DIR/HTSeq_Count_union.txt"

rm "$HT_DIR/"*_counts.txt "$HT_DIR/gene_column.txt"


# Merge all HTSeq-count outputs
echo "Merging htseq-count files into a count matrix..."

# Extract gene names (assuming they are in column 1) from the first file
cut -f1 "$HT_DIR/$(ls $HT_DIR | head -n1)" > "$HT_DIR/gene_column.txt"

# Extract count columns and paste them together
for FILE in "$HT_DIR"/*.htseq.txt; do
    SAMPLE=$(basename "$FILE" .htseq.txt)
    cut -f2 "$FILE" > "$HT_DIR/${SAMPLE}_counts.txt"
done

# Combine all into one file
paste "$HT_DIR/gene_column.txt" "$HT_DIR/"*_counts.txt > "$HT_DIR/HTSeq_Count_union.txt"

# Add header
HEADER="Gene"
for FILE in "$HT_DIR/"*_counts.txt; do
    SAMPLE=$(basename "$FILE" _counts.txt)
    HEADER+="\t$SAMPLE"
done
sed -i "1s/^/${HEADER}\n/" "$HT_DIR/HTSeq_Count_union.txt"

# Clean intermediate count columns
rm "$HT_DIR/"*_counts.txt "$HT_DIR/gene_column.txt"

###########
echo "Running featureCounts..."
if [[ $MODE == "PE" ]]; then
    echo "Running in Paired-End mode..."
    featureCounts -p -g gene_id -T $THREADS -O -M -a "$GTF"/*.gtf -o Quantification/Hisat2/FC/FC_Count.txt Mapping/Hisat2/*.bam
else
    echo "Running in Single-End mode..."
    featureCounts -g gene_id -T $THREADS -O -M -a "$GTF"/*.gtf -o Quantification/Hisat2/FC/FC_Count.txt Mapping/Hisat2/*.bam
fi

# Optional DE analysis
# Rscript DE_edgr.r Sampleinfo.txt FC_Count.txt logFC.txt Pval.txt
# Rscript DE_DESe.r Sampleinfo.txt FC_Count.txt logFC.txt Pval.txt

echo "All processing completed."

#!/bin/bash
echo "Running RSEM quantification..."
# Create RSEM reference if not present
# RSEM reference preparation (run once)
mkdir -p Indices/RSEM/Bowtie
#rsem-prepare-reference --hisat2 --gtf $GTF/*.gtf $FASTA/*.fna Indices/RSEM/Hisat2/rsem_ref 
echo "RSEM reference preparation completed."
####

