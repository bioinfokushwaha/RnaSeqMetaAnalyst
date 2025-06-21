
# Index genome (uncomment if needed)
#STAR --runMode genomeGenerate --runThreadN $THREADS  --genomeDir Indices/STAR/  --genomeFastaFiles $FASTA/*.fna --sjdbGTFfile $GTF/*.gff    --sjdbOverhang 100

echo "Genome indexing completed."

if [[ $MODE == "PE" ]]; then
    echo "Running in Paired-End mode..."

    for i in "$READ_DIR"/*_R1.fastq; do
        SAMPLE=$(basename "$i" _R1.fastq)
        R1="$READ_DIR/${SAMPLE}_R1.fastq"
        R2="$READ_DIR/${SAMPLE}_R2.fastq"

        echo "Processing sample: $SAMPLE"
        STAR --genomeDir Indices/STAR \
            --sjdbGTFfile $GTF/*.gff \
            --readFilesIn "$R1" "$R2" \
            --outSAMtype BAM SortedByCoordinate \
           --outFilterMultimapNmax 10 \
            --outSAMunmapped Within \
            --quantMode TranscriptomeSAM \
            --twopassMode Basic \
            --runThreadN $THREADS \
            --outFileNamePrefix Mapping/STAR/"$SAMPLE"_ \
            --sjdbOverhang 100 \
            --limitIObufferSize 150000000 150000000 \
            --outSAMattrRGline ID:$SAMPLE PL:illumina PU:${SAMPLE}-201710 LIB:KAPA SM:$SAMPLE LB:${SAMPLE}-3RNAseq

        echo "Indexing BAM for $SAMPLE..."
        samtools index Mapping/STAR/"${SAMPLE}_Aligned.sortedByCoord.out.bam"
    done

else
    echo "Running in Single-End mode..."

   for i in "$READ_DIR"/*.fastq; do
        SAMPLE=$(basename "$i" .fastq)
        R="$READ_DIR/${SAMPLE}.fastq"

        echo "Processing sample: $SAMPLE"
        STAR --genomeDir Indices/STAR \
            --sjdbGTFfile $GTF/*.gff \
            --readFilesIn "$R" \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 10 \
            --quantMode TranscriptomeSAM \
            --outSAMunmapped Within \
             --twopassMode Basic \
            --runThreadN $THREADS \
            --outFileNamePrefix Mapping/STAR/"$SAMPLE"_ \
            --sjdbOverhang 100 \
            --limitIObufferSize 150000000 150000000 \
            --outSAMattrRGline ID:$SAMPLE PL:illumina PU:${SAMPLE}-201710 LIB:KAPA SM:$SAMPLE LB:${SAMPLE}-3RNAseq

        echo "Indexing BAM for $SAMPLE..."
        samtools index Mapping/STAR/"${SAMPLE}_Aligned.sortedByCoord.out.bam"
    done
fi
#
######@####
##########################################
# Index BAM files before HTSeq
##########################################
echo "Indexing BAM files..."
for BAM in Mapping/STAR/*_Aligned.sortedByCoord.out.bam; do
    if [ -f "$BAM" ] && [ ! -f "${BAM}.bai" ]; then
        echo "Indexing $BAM"
        samtools index "$BAM"
    fi
done
######### HTSeq-count ##########
echo "Running htseq-count sample-wise..."

STAR_DIR="Quantification/STAR/HT"
mkdir -p "$STAR_DIR"

# Generate counts
for BAM in Mapping/STAR/*_Aligned.sortedByCoord.out.bam; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "Counting for sample: $SAMPLE"
    htseq-count -f bam -m union "$BAM" "$GTF"/*gtf > "$STAR_DIR/${SAMPLE}.htseq.txt"
done

# Merge htseq-count results
echo "Merging htseq-count files into a count matrix..."

cut -f1 "$(ls $STAR_DIR/*.htseq.txt | head -n1)" > "$STAR_DIR/gene_column.txt"

for FILE in "$STAR_DIR"/*.htseq.txt; do
    SAMPLE=$(basename "$FILE" .htseq.txt)
  cut -f2 "$FILE" > "$STAR_DIR/${SAMPLE}_counts.txt"
done

paste "$STAR_DIR/gene_column.txt" "$STAR_DIR/"*_counts.txt > "$STAR_DIR/HTSeq_Count_union.txt"

#### Add header
HEADER="Gene"
for FILE in "$STAR_DIR"/*_counts.txt; do
    SAMPLE=$(basename "$FILE" _counts.txt)
    HEADER+="\t$SAMPLE"
done
sed -i "1s/^/${HEADER}\n/" "$STAR_DIR/HTSeq_Count_union.txt"
# Clean up
rm "$STAR_DIR/"*_counts.txt "$STAR_DIR/gene_column.txt"

###### featureCounts ##########
echo "Running featureCounts..."
if [[ $MODE == "PE" ]]; then
    echo "Running in Paired-End mode..."
    featureCounts -p -g gene_id -T $THREADS -O -M -a "$GTF"/*gtf -o Quantification/STAR/FC/FC_Count.txt Mapping/STAR/*_Aligned.sortedByCoord.out.bam
else
    echo "Running in Single-End mode..."
    featureCounts -g gene_id -T $THREADS -O -M -a "$GTF"/*gtf -o Quantification/STAR/FC/FC_Count.txt Mapping/STAR/*_Aligned.sortedByCoord.out.bam

fi
############ rsem reference building 
# STAR-based RSEM reference preparation
mkdir -p Indices/RSEM/STAR
echo "Running RSEM quantification..."
# Create RSEM reference if not present
# RSEM reference preparation (run once)
mkdir -p Indices/RSEM/STAR
rsem-prepare-reference --gtf $GTF/*.gtf $FASTA/*.fna Indices/RSEM/STAR/rsem_ref --star
echo "RSEM reference preparation completed."
####


#mkdir -p Quantification/STAR/RSEM
echo "Running RSEM quantification..."

if [[ $MODE == "PE" ]]; then
    echo "Running RSEM in Paired-End mode..."
    for R1 in "$READ_DIR"/*_R1.fastq; do
        SAMPLE=$(basename "$R1" _R1.fastq)
        R2="$READ_DIR/${SAMPLE}_R2.fastq"

        echo "Quantifying sample: $SAMPLE"
        rsem-calculate-expression \
            --paired-end \
            --star \
            --star-path /usr/bin/ \
            --estimate-rspd \
            --append-names \
            -p $THREADS \
            "$R1" "$R2" \
            Indices/RSEM/STAR/rsem_ref \
            Quantification/STAR/RSEM/"$SAMPLE"
    done

else
    echo "Running RSEM in Single-End mode..."
    for R in "$READ_DIR"/*.fastq; do
        SAMPLE=$(basename "$R" .fastq)

        echo "Quantifying sample: $SAMPLE"
        rsem-calculate-expression \
            --star \
            --star-path /usr/bin/ \
            --estimate-rspd \
            --append-names \
            -p $THREADS \
            "$R" \
            Indices/RSEM/STAR/rsem_ref \
            Quantification/STAR/RSEM/"$SAMPLE"
    done
fi
##
#!/bin/bash
# Directory where RSEM results are stored
RSEM_DIR="Quantification/STAR/RSEM"
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
######
# Define destination folder
DEST="DEG"
mkdir -p "$DEST"

echo "Copying STAR merged files..."
cp Quantification/STAR/FC/FC_Count.txt "$DEST/S_FC_Count.txt"
cp Quantification/STAR/HT/HTSeq_Count_union.txt "$DEST/S_HTSeq_Count_union.txt"
cp Quantification/STAR/RSEM/rsem_gene_counts_matrix.txt "$DEST/S_rsem_gene_counts_matrix.txt"

