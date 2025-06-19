
#!/bin/bash

THREADS=8
MODE="SE"
GENOME_DIR="data/genome"
READ_DIR="data/reads"
GTF="data/genome"
FASTA="data/genome"

echo "DEBUG: MODE is [$MODE]"
echo "GENOME_DIR is [$GENOME_DIR]"

mkdir -p Indices/Bowtie Mapping/Bowtie Quantification/Bowtie/{HT,FC,RSEM}

# Index genome (uncomment if needed)
#bowtie2-build "$FASTA"/*fna Indices/Bowtie/Bowtie

#echo "Genome indexing completed."

#if [[ $MODE == "PE" ]]; then
#    echo "Running in Paired-End mode..."
#
#    for i in "$READ_DIR"/*_R1.fastq; do
#        SAMPLE=$(basename "$i" _R1.fastq)
#        R1="$READ_DIR/${SAMPLE}_R1.fastq"
#        R2="$READ_DIR/${SAMPLE}_R2.fastq"
#
#        echo "Processing sample: $SAMPLE"
#        bowtie2 -p $THREADS -x Indices/Bowtie/Bowtie \
#            -1 "$R1" -2 "$R2" \
#            -S Mapping/Bowtie/${SAMPLE}.sam \
#            2> Mapping/Bowtie/${SAMPLE}.log
#
#        echo "Converting SAM to BAM for $SAMPLE..."
#        samtools view -@ $THREADS -b Mapping/Bowtie/${SAMPLE}.sam > Mapping/Bowtie/${SAMPLE}.bam
#        rm Mapping/Bowtie/${SAMPLE}.sam
#
#        echo "Sorting and indexing BAM for $SAMPLE..."
#        samtools sort -@ $THREADS -o Mapping/Bowtie/${SAMPLE}_sorted.bam Mapping/Bowtie/${SAMPLE}.bam
#        mv Mapping/Bowtie/${SAMPLE}_sorted.bam Mapping/Bowtie/${SAMPLE}.bam
#        samtools index Mapping/Bowtie/${SAMPLE}.bam
#    done
#
#else
#    echo "Running in Single-End mode..."
#
#    for i in "$READ_DIR"/*.fastq; do
#        SAMPLE=$(basename "$i" .fastq)
#        R1="$READ_DIR/${SAMPLE}.fastq"
#
#        echo "Processing sample: $SAMPLE"
#        bowtie2 -p $THREADS -x Indices/Bowtie/Bowtie \
#            -U "$R1" \
#            -S Mapping/Bowtie/${SAMPLE}.sam \
#            2> Mapping/Bowtie/${SAMPLE}.log
#
#        echo "Converting SAM to BAM for $SAMPLE..."
#        samtools view -@ $THREADS -b Mapping/Bowtie/${SAMPLE}.sam > Mapping/Bowtie/${SAMPLE}.bam
#        rm Mapping/Bowtie/${SAMPLE}.sam
#
#        echo "Sorting and indexing BAM for $SAMPLE..."
#        samtools sort -@ $THREADS -o Mapping/Bowtie/${SAMPLE}_sorted.bam Mapping/Bowtie/${SAMPLE}.bam
#        mv Mapping/Bowtie/${SAMPLE}_sorted.bam Mapping/Bowtie/${SAMPLE}.bam
#        samtools index Mapping/Bowtie/${SAMPLE}.bam
#    done
#fi

########## HTSeq-count ##########
#echo "Running htseq-count sample-wise..."
#
#BOW_DIR="Quantification/Bowtie/HT"
#mkdir -p "$BOW_DIR"
#
## Generate counts
#for BAM in Mapping/Bowtie/*.bam; do
#    SAMPLE=$(basename "$BAM" .bam)
#    echo "Counting for sample: $SAMPLE"
#    htseq-count -f bam -m union "$BAM" "$GTF"/*gtf > "$BOW_DIR/${SAMPLE}.htseq.txt"
#done

# Merge htseq-count results
#echo "Merging htseq-count files into a count matrix..."
#
#cut -f1 "$(ls $BOW_DIR/*.htseq.txt | head -n1)" > "$BOW_DIR/gene_column.txt"
#
#for FILE in "$BOW_DIR"/*.htseq.txt; do
#    SAMPLE=$(basename "$FILE" .htseq.txt)
#    cut -f2 "$FILE" > "$BOW_DIR/${SAMPLE}_counts.txt"
#done
#
#paste "$BOW_DIR/gene_column.txt" "$BOW_DIR/"*_counts.txt > "$BOW_DIR/HTSeq_Count_union.txt"

# Add header
#HEADER="Gene"
#for FILE in "$BOW_DIR"/*_counts.txt; do
#    SAMPLE=$(basename "$FILE" _counts.txt)
#    HEADER+="\t$SAMPLE"
#done
#sed -i "1s/^/${HEADER}\n/" "$BOW_DIR/HTSeq_Count_union.txt"
#
# Clean up
#rm "$BOW_DIR/"*_counts.txt "$BOW_DIR/gene_column.txt"
#
##### featureCounts ##########
echo "Running featureCounts..."
if [[ $MODE == "PE" ]]; then
    echo "Running in Paired-End mode..."
    featureCounts -p -g gene_id -T $THREADS -O -M -a "$GTF"/*gtf -o Quantification/Bowtie/FC/FC_Count.txt Mapping/Bowtie/*.bam
else
    echo "Running in Single-End mode..."
    featureCounts -g gene_id -T $THREADS -O -M -a "$GTF"/*gtf -o Quantification/Bowtie/FC/FC_Count.txt Mapping/Bowtie/*.bam

fi

echo "All processing completed."

########### RSEM Quantification ##########
echo "Running RSEM quantification..."
echo "Running RSEM quantification..."
# Create RSEM reference if not present
# RSEM reference preparation (run once)
mkdir -p Indices/RSEM/Bowtie
rsem-prepare-reference --gtf $GTF/*.gtf $FASTA/*.fna Indices/RSEM/Bowtie/rsem_ref --bowtie2
echo "RSEM reference preparation completed."
####
########## RSEM ##########
echo "Running RSEM quantification..."

mkdir -p Quantification/Bowtie/RSEM

if [[ $MODE == "PE" ]]; then
    echo "Running RSEM in Paired-End mode..."
    for R1 in "$READ_DIR"/*_R1.fastq; do
        SAMPLE=$(basename "$R1" _R1.fastq)
        R2="$READ_DIR/${SAMPLE}_R2.fastq"

        echo "Quantifying sample: $SAMPLE"
        rsem-calculate-expression \
            --paired-end \
            --bowtie2 \
            --bowtie2-path /usr/bin/ \
            --estimate-rspd \
            --append-names \
            -p $THREADS \
            "$R1" "$R2" \
            Indices/RSEM/Bowtie/rsem_ref \
            Quantification/Bowtie/RSEM/"$SAMPLE"
    done

else
    echo "Running RSEM in Single-End mode..."
    for R in "$READ_DIR"/*.fastq; do
        SAMPLE=$(basename "$R" .fastq)

        echo "Quantifying sample: $SAMPLE"
        rsem-calculate-expression \
            --bowtie2 \
            --bowtie2-path /usr/bin/ \
            --estimate-rspd \
            --append-names \
            -p $THREADS \
            "$R" \
            Indices/RSEM/Bowtie/rsem_ref \
            Quantification/Bowtie/RSEM/"$SAMPLE"
    done
fi


########
#!/bin/bash

# Directory where RSEM results are 
#!/bin/bash

# Directory where RSEM results are stored
RSEM_DIR="Quantification/Bowtie/RSEM"
OUTPUT="$RSEM_DIR/rsem_gene_counts_matrix.txt"

# Get gene_id column from the first file (excluding header)
FIRST_FILE=$(ls "$RSEM_DIR"/*.genes.results | head -n1)
cut -f1 "$FIRST_FILE" | tail -n +2 > "$RSEM_DIR/gene_ids.txt"

# Create a temporary directory for counts
TEMP_DIR="$RSEM_DIR/temp_counts"
mkdir -p "$TEMP_DIR"

# Extract raw counts (5th column) from each file
for FILE in "$RSEM_DIR"/*.genes.results; do
    SAMPLE=$(basename "$FILE" .genes.results)
    echo "Processing $SAMPLE..."
    cut -f5 "$FILE" | tail -n +2 > "$TEMP_DIR/${SAMPLE}.txt"
done

# Merge gene_ids with all sample count files
paste "$RSEM_DIR/gene_ids.txt" "$TEMP_DIR"/*.txt > "$OUTPUT"

# Create header row
HEADER="Geneid"
for FILE in "$TEMP_DIR"/*.txt; do
    SAMPLE=$(basename "$FILE" .txt)
    HEADER="$HEADER,$SAMPLE"
done

# Prepend header to the final matrix
(echo "$HEADER"; cat "$OUTPUT") > "$RSEM_DIR/tmpfile" && mv "$RSEM_DIR/tmpfile" "$OUTPUT"

# Clean up temporary files
rm -r "$TEMP_DIR" "$RSEM_DIR/gene_ids.txt"

echo "✅ Cleaned RSEM gene count matrix saved to: $OUTPUT"

