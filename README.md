# RnaSeqMetaAnalyst
![image](https://github.com/user-attachments/assets/3f5cde0a-61b9-4ed2-af77-9ffa789507de)

A bioinformatics pipeline to perfrom Meta-analysis of RNA-Seq data using uniform processing.
This repository contains a complete, reproducible RNA-Seq analysis pipeline packaged in a Docker container. It automates alignment (HISAT2, STAR, Bowtie2), Quantificatin (featureCount, htseq-count, rsem) and differential gene expression analysis using DESeq2 and edgeR.

📦 Features
🔁 End-to-end automation: fromAlignment to DEG analysis

🐳 Dockerized: no dependency issues

🧬 Supports HISAT2, STAR, Bowtie2 for alignment

📊 Generates Excel summaries for DEG results (edgeR & DESeq2)

📁 Output organized in user-friendly structure


# Requirements
## Files

## Tools Used
The following tools and packages were used in this RNA-seq analysis pipeline:

    FastQC – for raw sequence quality control

    Fastp – for trimming and filtering of reads

    MultiQC – for aggregating QC reports

    STAR – for spliced alignment to the genome

    Hisat2 – for efficient alignment of RNA-seq reads

    Bowtie2 – for general sequence alignment

    HTSeq – for counting reads mapped to genes

    FeatureCounts – for efficient read summarization

    RSEM – for transcript quantification

    R Packages:

        DESeq2 – differential expression analysis

        edgeR – differential expression analysis

        dplyr – data manipulation

        openxlsx – Excel output handling
## 📁 Project Directory Structure

```
project_root/
│
├── reference/                     # Reference genome and annotation
│   ├── genome.fa                  # Genome FASTA file
│   ├── genome.gtf
|   |               # Annotation GTF file
│   ├── hisat_index/               # HISAT2 index
│   ├── bowtie2_index/             # Bowtie2 index
│   └── star_index/                # STAR index
│
├── data/                          # Raw sequencing files
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   └── ...
│
├── sampleinfo/
│   └── sampleinfo.txt             # Sample condition metadata
│
├── DEG/                           # DEG results (edgeR/DESeq2)
│   └── (generated output folders per comparison)
│
├── bioinfo_pipeline/
│   ├── Dockerfile                 # Docker build file
│   └── scripts/                   # Pipeline scripts
│       ├── master.py              # Master controller script
│       ├── run_all_edgeR.sh       # edgeR batch runner
│       ├── run_all_deseq2.sh      # DESeq2 batch runner
│       ├── run_edgeR_analysis.R   # edgeR DEG logic
│       ├── run_deseq_222.R        # DESeq2 DEG logic
│       ├── hisat_ok.sh            # HISAT2 aligner
│       ├── bowtie.sh              # Bowtie2 aligner
│       ├── star.sh                # STAR aligner
│       └── (any other helper scripts)
```

   ## 1) Clone the repository
```
git clone https://github.com/bioinfokushwaha/RnaSeqMetaAnalyst.git
```

## 2) Move into the pipeline directory
```
cd RnaSeqMetaAnalyst/bioinfo_pipeline
```

## 3) Build the Docker image
```
docker build -t rnaseq-metaanalyst .
```

## 4) Run the pipeline (Adjust paths if needed as per your directory. Also the MODE as single or paired end)
```
docker run --rm -it -e THREADS=8 -e MODE=SE -e GENOME_DIR=/genome/ -e READ_DIR=/raw/ -e TRIM_DIR=/Trim -e GTF=/genome/GCF_003369695.1_UOA_Brahman_1_genomic.gtf -e FASTA=/genome/GCF_003369695.1_UOA_Brahman_1_genomic.fna -e GFF=/genome/GCF_003369695.1_UOA_Brahman_1_genomic.gff -v "$(pwd)/raw/":/raw/ -v "$(pwd)/Trim/":/Trim/ -v "$(pwd)/genome/":/genome/ -v "$(pwd)/":/data rnaseq-metaanalyst bash /opt/project/scripts/master.sh


```
Note: MODE="SE" if single end files; MODE="PE" if paired end files; The command given above is an example to understand. Change the input file names as required
