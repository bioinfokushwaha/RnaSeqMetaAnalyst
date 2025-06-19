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

## Tools 
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
<pre> 
   # Clone the repository
git clone https://github.com/bioinfokushwaha/RnaSeqMetaAnalyst.git

# Move into the pipeline directory
cd RnaSeqMetaAnalyst/bioinfo_pipeline

# Build the Docker image
docker build -t rnaseq-metaanalyst .

# Run the pipeline (adjust paths if needed)
docker run -it -v "$PWD/data":/data -v "$PWD/counts":/counts  -v "$PWD/output":/output  rnaseq-metaanalyst python /opt/bioinfo/scripts/master.py
 </pre>
