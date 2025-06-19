# RnaSeqMetaAnalyst
![image](https://github.com/user-attachments/assets/3f5cde0a-61b9-4ed2-af77-9ffa789507de)

A bioinformatics pipeline to perfrom Meta-analysis of RNA-Seq data using uniform processing.
This repository contains a complete, reproducible RNA-Seq analysis pipeline packaged in a Docker container. It automates alignment (HISAT2, STAR, Bowtie2) and differential gene expression analysis using DESeq2 and edgeR.
# Requirements
## Files
project_root/

├── reference/

│   ├── genome.fa

│   ├── genome.gtf

│   ├── hisat_index/

│   ├── bowtie2_index/

│   └── star_index/

│
|
├── data/

│   ├── sample1_R1.fastq.gz

│   ├── sample1_R2.fastq.gz

│   └── ...

│

├── sampleinfo/

│   └── sampleinfo.txt

│

├── DEG/
   
│   └── (output will be written here)

│

├── bioinfo_pipeline/
|

│   ├── Dockerfile

│   └── scripts/

│       └── (all your .sh, .R, and .py files)

## Tools 
