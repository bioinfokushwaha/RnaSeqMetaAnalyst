# RnaSeqMetaAnalyst
![image](https://github.com/user-attachments/assets/3f5cde0a-61b9-4ed2-af77-9ffa789507de)

A bioinformatics pipeline to perfrom Meta-analysis of RNA-Seq data using uniform processing.
This repository contains a complete, reproducible RNA-Seq analysis pipeline packaged in a Docker container. It automates alignment (HISAT2, STAR, Bowtie2) and differential gene expression analysis using DESeq2 and edgeR.
# Requirements
## Files
bioinfo_pipeline/
│
├── Dockerfile             # Builds the environment
├── environment.yml        # Conda dependencies
├── scripts/
│   ├── hisat_ok.sh
│   ├── bowtie2222_ok.sh
│   ├── run_all_deseq2.sh
│   ├── run_all_edgeR.sh
│   ├── run_deseq2_analysis.R
│   ├── run_edgeR_analysis.R
│   └── controller.py      # Or master.py – controls the full pipeline

## Tools 
