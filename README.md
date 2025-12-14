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
# RNA-Seq Meta-Analyst Pipeline - Complete Documentation

A comprehensive bioinformatics pipeline for RNA-sequencing analysis with support for multiple alignment tools and quantification methods. This guide is designed for beginners and advanced users alike.

---

## Table of Contents

1. [Overview](#overview)
2. [System Requirements](#system-requirements)
3. [Installation](#installation)
4. [Pipeline Structure](#pipeline-structure)
5. [Detailed Workflow](#detailed-workflow)
6. [Running the Pipeline](#running-the-pipeline)
7. [Output Files](#output-files)
8. [Troubleshooting](#troubleshooting)
9. [Advanced Usage](#advanced-usage)

---

## Overview

### What is RNA-Seq Meta-Analyst?

RNA-Seq Meta-Analyst is an automated containerized pipeline that processes raw RNA-sequencing data from raw FASTQ files to differential gene expression results. It runs multiple alignment and quantification tools in parallel, allowing you to compare different analysis strategies.

### Key Features

- **Multiple Aligners**: HISAT2, Bowtie2, and STAR
- **Multiple Quantifiers**: featureCounts, HTSeq, and RSEM
- **Automated Processing**: From quality control to differential expression analysis
- **Parallel Execution**: Aligners run simultaneously to save time
- **Comprehensive Output**: DEG lists, volcano plots, PCA analysis, and more
- **Docker Containerized**: No dependency conflicts, works on any system

### Pipeline Combinations Generated

The pipeline produces results from multiple analysis approaches:

```
Bowtie2 + featureCounts + edgeR
Bowtie2 + featureCounts + DESeq2
Bowtie2 + HTSeq + edgeR
Bowtie2 + HTSeq + DESeq2
Bowtie2 + RSEM + edgeR
Bowtie2 + RSEM + DESeq2

HISAT2 + featureCounts + edgeR
HISAT2 + featureCounts + DESeq2

STAR + featureCounts + edgeR
STAR + featureCounts + DESeq2
STAR + HTSeq + edgeR
STAR + HTSeq + DESeq2
STAR + RSEM + edgeR
STAR + RSEM + DESeq2

Total: 16 analysis combinations
```

---

## System Requirements

### Minimum Hardware Specifications

| Requirement | Minimum | Recommended |
|---|---|---|
| CPU Cores | 4 | 16+ |
| RAM | 16 GB | 64 GB |
| Disk Space | 500 GB | 1-2 TB |
| GPU | Not required | Not required |

### Software Prerequisites

- **Docker**: Version 20.10 or newer
  - [Install Docker](https://docs.docker.com/get-docker/)
  - Verify installation: `docker --version`
  
- **Linux/macOS/Windows (with WSL2)**:
  - Linux: Native support
  - macOS: Docker Desktop
  - Windows: Docker Desktop with WSL2

### Check Docker Installation

```bash
# Verify Docker is running
docker ps

# You should see: CONTAINER ID   IMAGE   COMMAND   STATUS

# If Docker daemon is not running, start it:
# macOS: Applications > Docker > Open
# Linux: sudo systemctl start docker
```

---

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/bioinfokushwaha/RnaSeqMetaAnalyst/tree/main/bioinfo_pipeline
cd rnaseq-meta-analyst
```

### Step 2: Build the Docker Image

```bash
# Build the container (this takes 10-30 minutes)
docker build -t rnaseqmetaanalyst:latest .

# Check if the build was successful
docker images | grep rnaseqmetaanalyst
```

### Step 3: Prepare Your Data Structure

Create a project directory with the following structure:

```
my_rnaseq_project/
├── data/
│   ├── raw/                          # Your raw FASTQ files
│   │   ├── sample1_1.fastq.gz
│   │   ├── sample1_2.fastq.gz
│   │   ├── sample2_1.fastq.gz
│   │   ├── sample2_2.fastq.gz
│   │   └── ...
│   ├── genome/                       # Reference genome files
│   │   ├── genomic.fasta
│   │   ├── genomic.gtf
│   │   └── genomic.gff
│   └── sampleinfo.txt                # Sample metadata (see format below)
└── results/                          # Output directory (created automatically)
```

### Step 4: Create Sample Information File

Create `data/sampleinfo.txt` with the following format:

```
SampleName	CONDITION
sample1_1	Control
sample1_2	Control
sample2_1	Treatment
sample2_2	Treatment
```

**Column Definitions:**
- **SampleName**: Must match your FASTQ filename (without extensions)
- **CONDITION**: The experimental condition/group (e.g., Control, Treatment, Disease)

**Example for paired-end reads:**
```
sample_1	Control
sample_2	Control
sample_1	Treatment
sample_2	Treatment
```

---

## Pipeline Structure

### Directory Organization

```
Pipeline Components:
├── Quality Control (QC)
│   ├── FastQC (raw reads)
│   ├── FastQC (trimmed reads)
│   ├── fastp (trimming)
│   └── MultiQC (summary)
│
├── Alignment (Parallel)
│   ├── HISAT2
│   │   ├── Index building
│   │   └── Alignment → BAM files
│   ├── Bowtie2
│   │   ├── Index building
│   │   └── Alignment → BAM files
│   └── STAR
│       ├── Index building
│       └── Alignment → BAM files
│
├── Quantification
│   ├── featureCounts
│   ├── HTSeq
│   └── RSEM (Bowtie2 & STAR only)
│
└── Differential Expression (DEG)
    ├── edgeR analysis
    ├── DESeq2 analysis
    ├── Volcano plots
    ├── PCA analysis
    └── Gene lists (UP/DOWN regulated)
```

### Tools Explained

#### Alignment Tools

**HISAT2** (Hierarchical Indexing for Spliced Alignment of Transcripts)
- Fast, memory-efficient genome aligner
- Best for: Mammalian genomes
- Does NOT support RSEM quantification
- Output: SAM/BAM files

**Bowtie2**
- General-purpose short-read aligner
- Balanced accuracy and speed
- Supports RSEM quantification
- Output: SAM/BAM files

**STAR** (Splicing Transcription Alignment from RNA-seq)
- Specialized for RNA-seq with accurate splice junction detection
- Best for: Accurate splice site detection
- Supports RSEM quantification
- Output: BAM files

#### Quantification Tools

**featureCounts**
- Fast read counting at exon/gene level
- Counts reads overlapping gene features
- Output: Raw gene counts

**HTSeq**
- Rigorous gene quantification with strict mode
- Uses "union" mode: handles ambiguous reads carefully
- Output: Gene counts

**RSEM** (RNA-Seq by Expectation-Maximization)
- Estimates isoform-level and gene-level expression
- Handles ambiguous reads probabilistically
- Produces: Expected counts, normalized values
- Only compatible with: Bowtie2, STAR

#### Differential Expression Tools

**edgeR**
- Uses negative binomial distribution
- Good for: Smaller sample sizes
- Produces: Fold change, p-values, FDR

**DESeq2**
- Uses Bayesian approach for normalization
- Good for: General RNA-seq analysis
- Produces: Fold change, p-values, FDR
- Note: Requires integer counts (RSEM counts are auto-converted)

---

## Detailed Workflow

### Step 1: Quality Control (QC)

**What happens:**
- FastQC analyzes raw FASTQ files for quality metrics
- fastp trims low-quality reads and adapters
- FastQC re-analyzes trimmed reads
- MultiQC creates a comprehensive summary report

**Key metrics to check:**
- Per-base quality scores (should be >30)
- GC content distribution
- Adapter contamination (should be <5%)
- Read length distribution

**Output files:**
```
qc/
├── fastqc_raw/          # Initial quality reports
├── fastqc_clean/        # Quality after trimming
├── fastp_reports/       # Trimming statistics
└── multiqc/             # Summary report
    └── multiqc_report.html  # View in browser
```

### Step 2: Alignment (Parallel Processing)

**What happens:**

1. **Index Building** (runs once, cached for reuse)
   - Creates searchable genome index for each aligner
   - Time: 1-30 minutes (depends on genome size)
   
2. **Read Alignment** (parallel across tools)
   - Maps each trimmed read to the genome
   - Generates alignment files (SAM format)
   
3. **SAM to BAM Conversion**
   - Compresses alignments to binary format
   - Sorts by genomic coordinate
   - Creates index files for fast access

**Alignment statistics:**
- % of reads successfully aligned
- Uniquely mapped reads
- Multi-mapped reads
- Unmapped reads

**Output files:**
```
Mapping/
├── Bowtie/
│   └── *_sorted.bam          # Coordinate-sorted alignments
├── Hisat2/
│   └── *_sorted.bam
└── STAR/
    └── *_sorted.bam
```

### Step 3: Quantification

**What happens:**

1. **featureCounts** (gene-level counting)
   - Counts reads overlapping known genes
   - Handles paired-end reads and strand information
   
2. **HTSeq** (alternative gene-level counting)
   - Uses union mode for ambiguous reads
   - Provides comparison to featureCounts
   
3. **RSEM** (probabilistic quantification)
   - Estimates isoform-level expression
   - Aggregates to gene-level
   - Handles multi-mapped reads probabilistically

**Output files:**
```
Quantification/
├── Bowtie/
│   ├── FC/          # featureCounts output
│   ├── HT/          # HTSeq output
│   └── RSEM/        # RSEM output
├── Hisat2/
│   ├── FC/
│   └── HT/          # Note: No RSEM for HISAT2
└── STAR/
    ├── FC/
    ├── HT/
    └── RSEM/
```

### Step 4: Data Cleaning

**What happens:**

The pipeline cleans count matrices:
- Removes batch metadata (sample names standardized)
- Removes file extensions (.bam, .fastq, etc.)
- Converts to CSV format
- RSEM: Maintains decimal values (DESeq2 automatically handles conversion)
- Validates: Ensures all expected samples present

**Output files:**
```
DEG/
├── *_FC_counts_clean.txt       # featureCounts cleaned
├── *_HTSeq_counts_clean.txt    # HTSeq cleaned
├── *_RSEM_counts_clean.txt     # RSEM cleaned
└── sampleinfo.txt              # Sample metadata
```

### Step 5: Differential Expression Analysis

**edgeR Analysis:**
1. Filters low-count genes
2. Calculates normalization factors
3. Estimates gene-wise dispersion
4. Performs statistical testing (GLM-based)
5. Generates volcano plots

**DESeq2 Analysis:**
1. Filters low-count genes
2. Estimates size factors (normalization)
3. Estimates gene-wise dispersion
4. Performs statistical testing (Wald test)
5. Generates volcano plots

**For each tool:**
- PCA plot (shows sample clustering)
- Normalized count matrix
- Housekeeping genes analysis (if available)
- Volcano plots (high resolution: 3400×3400 @ 600 DPI)
- DEG lists filtered for non-zero expression

**Output files:**
```
DEG/
└── [tool_name]_[quantifier]/
    ├── [name]_Normalized_Counts.xlsx
    ├── [name]_PCA_Scores.xlsx
    ├── [name]_PCA_Plot.png
    ├── [name]_HKG_NonZero_Normalized.xlsx
    ├── [name]_[Condition1]_vs_[Condition2]_UP.xlsx
    ├── [name]_[Condition1]_vs_[Condition2]_DOWN.xlsx
    ├── [name]_[Condition1]_vs_[Condition2]_Volcano.png
    └── [name]_DEG_Summary.xlsx
```

---

## Running the Pipeline

### Quick Start (5 minutes)

```bash
# Navigate to the pipeline directory
cd /path/to/my_rnaseq_project

# Run the pipeline
./run_pipeline.sh \
  --project-dir /path/to/my_rnaseq_project \
  --mode PE \
  --threads 16
```

### Detailed Command Reference

```bash
./run_pipeline.sh [OPTIONS]
```

#### Required Arguments

```bash
--project-dir <path>    # Absolute path to your project directory
--mode <SE|PE>          # SE = Single-End, PE = Paired-End reads
```

#### Optional Arguments

```bash
--reads-dir <path>      # Raw FASTQ directory
                        # Default: <project-dir>/data/raw

--genome-dir <path>     # Genome files (FASTA, GTF, GFF)
                        # Default: <project-dir>/data/genome

--sample-info <path>    # Sample metadata file
                        # Default: <project-dir>/data/sampleinfo.txt

--output-dir <path>     # Results directory
                        # Default: <project-dir>/results

--threads <int>         # CPU threads to use (default: 8)
                        # Recommended: Total cores on your system

--index-dir <path>      # Pre-built indices (optional)

--image-name <name>     # Docker image name
                        # Default: rnaseqmetaanalyst:latest

--help                  # Show help message
```

#### Example Commands

**Single-End, 8 threads:**
```bash
./run_pipeline.sh \
  --project-dir /home/user/my_rnaseq \
  --mode SE \
  --threads 8
```

**Paired-End, 32 threads, custom paths:**
```bash
./run_pipeline.sh \
  --project-dir /home/user/my_rnaseq \
  --mode PE \
  --threads 32 \
  --reads-dir /mnt/data/fastq \
  --genome-dir /mnt/reference/genome \
  --output-dir /mnt/results/my_analysis
```

**With pre-built indices:**
```bash
./run_pipeline.sh \
  --project-dir /home/user/my_rnaseq \
  --mode PE \
  --threads 16 \
  --index-dir /home/shared_indices
```

### Monitoring Progress

```bash
# Watch the container in real-time (shows logs)
docker logs -f <CONTAINER_ID>

# Get the container ID if you forgot it
docker ps -a | grep rnaseqmetaanalyst

# Check file creation progress
ls -lh /path/to/results/DEG/
```

### Estimated Runtime

| Phase | Time | Notes |
|---|---|---|
| QC + Trimming | 30 min - 2 hrs | Depends on read count |
| Index Building | 30 min - 2 hrs | First run only, then cached |
| Alignment | 1-4 hrs | Parallel (3 aligners simultaneously) |
| Quantification | 30 min - 2 hrs | Depends on alignment size |
| DEG Analysis | 5-15 min | Per combination, usually fast |
| **Total** | **3-12 hours** | Single sample, typical hardware |

---

## Output Files

### Project Results Structure

```
my_rnaseq_project/results/
├── qc/
│   ├── fastqc_raw/
│   │   ├── sample1_1_fastqc.html       # View in browser
│   │   ├── sample1_2_fastqc.html
│   │   └── ...
│   ├── fastqc_clean/
│   │   ├── sample1_R1_fastqc.html
│   │   ├── sample1_R2_fastqc.html
│   │   └── ...
│   ├── fastp_reports/
│   │   ├── sample1_fastp.html
│   │   ├── sample1_fastp.json
│   │   └── ...
│   └── multiqc/
│       └── multiqc_report.html         # Master QC report
│
├── Trim/                               # Cleaned FASTQ files
│   └── *.fastq (uncompressed)
│
├── Mapping/
│   ├── Bowtie/
│   │   └── *_sorted.bam               # BAM files
│   ├── Hisat2/
│   │   └── *_sorted.bam
│   └── STAR/
│       └── *_sorted.bam
│
├── Quantification/
│   ├── Bowtie/
│   │   ├── FC/FC_Count.txt            # Gene counts
│   │   ├── HT/HTSeq_Count_union.txt
│   │   └── RSEM/*.genes.results       # RSEM outputs
│   ├── Hisat2/
│   │   └── ... (no RSEM)
│   └── STAR/
│       └── ... (includes RSEM)
│
├── DEG/                                # MAIN RESULTS FOLDER
│   ├── B_FC_counts_clean.txt           # Cleaned count matrices
│   ├── B_HTSeq_counts_clean.txt
│   ├── B_RSEM_counts_clean.txt
│   ├── H_FC_counts_clean.txt
│   ├── S_FC_counts_clean.txt
│   ├── S_HTSeq_counts_clean.txt
│   ├── S_RSEM_counts_clean.txt
│   │
│   ├── B_FC_counts_clean/              # DEG results by method
│   │   ├── B_FC_counts_clean_DESeq2_Normalized_Counts.xlsx
│   │   ├── B_FC_counts_clean_DESeq2_Control_vs_Treatment_UP.xlsx
│   │   ├── B_FC_counts_clean_DESeq2_Control_vs_Treatment_DOWN.xlsx
│   │   ├── B_FC_counts_clean_DESeq2_Control_vs_Treatment_Volcano.png
│   │   └── ...
│   ├── B_HTSeq_counts_clean/
│   ├── B_RSEM_counts_clean/
│   └── ... (directories for each analysis combination)
│
├── Indices/                            # Built genome indices (cached)
│   ├── Bowtie/genome/*.bt2
│   ├── Hisat2/genome/*.ht2
│   └── STAR/genome/
│
└── logs/
    ├── hisat2_TIMESTAMP.log            # Detailed logs per tool
    ├── bowtie_TIMESTAMP.log
    ├── star_TIMESTAMP.log
    ├── edgeR_TIMESTAMP.log
    ├── deseq2_TIMESTAMP.log
    └── pipeline_run.log                # Master log
```

### Key Output Files Explained

#### Count Matrices (in DEG/)
These are your primary quantification results:

```
S_FC_counts_clean.txt        # STAR + featureCounts
```

Format:
```
Geneid,Sample1,Sample2,Sample3,Sample4
ENSG00000000003,100,120,50,60
ENSG00000000005,0,0,5,10
ENSG00000000419,5000,4800,3000,3200
...
```

#### Differential Expression Results

**UP genes** (upregulated in first condition):
```
B_FC_counts_clean_edgeR_Control_vs_Treatment_UP.xlsx
```

Columns:
- `Trans`: Gene ID
- `logFC`: Log2 fold change
- `logCPM`: Log2 counts per million
- `LR`: Likelihood ratio
- `PValue`: P-value
- `FDR`: Adjusted p-value (false discovery rate)
- `NormalizedCount_Mean`: Average normalized expression

#### Visualization Files

**Volcano Plot** (PNG image):
```
B_FC_counts_clean_edgeR_Control_vs_Treatment_Volcano.png
```

Shows:
- X-axis: Log2 fold change (left=downregulated, right=upregulated)
- Y-axis: -log10(FDR) (higher=more significant)
- Red points: Significant DEGs (FDR < 0.05, |FC| > 1)
- Blue points: Non-significant genes

**PCA Plot** (PNG image):
```
B_FC_counts_clean_edgeR_PCA_Plot.png
```

Shows:
- Sample clustering by principal components
- Color-coded by experimental condition
- Helps identify batch effects

---

## Troubleshooting

### Common Issues

#### 1. "ERROR: docker command could not be found"

**Problem:** Docker is not installed or not in PATH.

**Solution:**
```bash
# Check if Docker is installed
which docker

# If not found, install Docker:
# macOS: https://docs.docker.com/desktop/install/mac-install/
# Linux: sudo apt-get install docker.io
# Windows: https://docs.docker.com/desktop/install/windows-install/

# Verify installation
docker --version  # Should show version
docker ps        # Should connect to daemon
```

#### 2. "ERROR: Reads directory not found"

**Problem:** FASTQ files are in wrong location.

**Solution:**
```bash
# Check your directory structure
ls -la /path/to/my_rnaseq_project/data/raw/

# You should see files like:
# sample1_1.fastq.gz
# sample1_2.fastq.gz

# If not, correct the --reads-dir path:
./run_pipeline.sh \
  --project-dir /path/to/my_rnaseq_project \
  --reads-dir /correct/path/to/fastq \
  --mode PE
```

#### 3. "ERROR: No matching samples found"

**Problem:** Sample names in FASTQ don't match sampleinfo.txt.

**Solution:**
```bash
# Check FASTQ filenames
ls data/raw/ | head -5

# Check sampleinfo.txt
cat data/sampleinfo.txt

# Example - if your files are:
# sample1_1.fastq.gz
# sample1_2.fastq.gz

# Your sampleinfo.txt should have:
# SampleName     CONDITION
# sample1        Control

# NOT "sample1_1" or "sample1_1.fastq.gz"
```

#### 4. "ERROR: gtf file not found"

**Problem:** Genome files missing or in wrong location.

**Solution:**
```bash
# Check genome directory
ls -la data/genome/

# You need:
# - genomic.fasta (or .fa)
# - genomic.gtf
# - genomic.gff

# If files have different names, edit run_pipeline.sh:
# Find these lines:
GTF="/data/genome/genomic.gtf"
FASTA="/data/genome/GCF_002263795.2_ARS-UCD1.3_genomic.fna"
GFF="/data/genome/genomic.gff"

# Update with your filenames
```

#### 5. Container runs but stops with no output

**Problem:** Container ran but produced no results.

**Solution:**
```bash
# Check container logs
docker logs <CONTAINER_ID>

# Look for specific error messages
# Common causes:
# - Insufficient disk space: df -h
# - Out of memory: top
# - Permission issues: chmod 777 data/

# Re-run with more verbose output:
./run_pipeline.sh ... 2>&1 | tee pipeline.log
```

#### 6. "No DEG result files were generated"

**Problem:** edgeR/DESeq2 analysis failed silently.

**Solution:**
```bash
# Check sample info format
cat data/sampleinfo.txt

# Ensure:
# - No extra spaces
# - Condition names match across samples
# - At least 2 samples per condition

# Example CORRECT format:
# SampleName     CONDITION
# sample1        ControlGroup
# sample2        ControlGroup
# sample3        TreatmentGroup

# Example INCORRECT format (SPACES):
# SampleName     CONDITION
# sample1 _R1    ControlGroup    # Extra underscore
# sample2_R2     ControlGroup    # Inconsistent naming
```

### Getting Help

1. **Check the logs:**
   ```bash
   tail -f /path/to/results/logs/pipeline_run.log
   ```

2. **Check per-tool logs:**
   ```bash
   cat /path/to/results/logs/hisat2_*.log
   cat /path/to/results/logs/bowtie_*.log
   cat /path/to/results/logs/deseq2_*.log
   ```

3. **Report issues** (include):
   - Full command you ran
   - Last 50 lines of relevant log file
   - Output from `docker ps -a`

---

## Advanced Usage

### Running with Pre-built Indices

If you're analyzing multiple projects with the same genome, build indices once:

```bash
# Create shared indices directory
mkdir -p /mnt/shared_indices

# Run once with index saving
./run_pipeline.sh \
  --project-dir /home/user/project1 \
  --mode PE \
  --index-dir /mnt/shared_indices

# Next projects reuse these indices (10x faster)
./run_pipeline.sh \
  --project-dir /home/user/project2 \
  --mode PE \
  --index-dir /mnt/shared_indices

# Indices are cached and reused automatically
```

### Comparing Analysis Methods

The pipeline generates 16 combinations. To compare them:

```bash
# After analysis completes, compare DEG overlap:
# All UP genes across methods
ls DEG/*/edgeR_*_UP.xlsx | wc -l

# Find genes consistent across methods
# (Use your favorite spreadsheet tool or R/Python)
```

### Memory and Performance Optimization

**For limited resources:**
```bash
# Use fewer threads
./run_pipeline.sh \
  --project-dir ... \
  --mode PE \
  --threads 4  # Use only 4 cores instead of 8
```

**For maximum speed:**
```bash
# Use all available cores
./run_pipeline.sh \
  --project-dir ... \
  --mode PE \
  --threads $(nproc)  # Auto-detect core count
```

### Using Different Genome Files

Edit `run_pipeline.sh` to update genome file paths:

```bash
# Inside run_pipeline.sh, find:
GTF="/data/genome/genomic.gtf"
FASTA="/data/genome/genomic.fasta"
GFF="/data/genome/genomic.gff"

# Change to your filenames:
GTF="/data/genome/your_genome.gtf"
FASTA="/data/genome/your_genome.fa"
GFF="/data/genome/your_genome.gff"
```

### Custom Quality Control Settings

Edit `quality_control.sh` to adjust:

```bash
# Trimming parameters (in fastp section)
fastp -w $THREADS \
    -i "$R1" -I "$R2" \
    -o ... \
    --cut_front 20 \        # Remove 20bp from 5' end
    --cut_tail 20 \         # Remove 20bp from 3' end
    --length_required 50 \  # Keep reads >50bp
    --qualified_quality_phred 30  # Quality threshold
```

---

## Understanding Statistical Results

### Fold Change (FC)

**What it means:** How much a gene is over/under-expressed

```
logFC = +2   →  2^2 = 4x upregulated (4 times higher)
logFC = +1   →  2^1 = 2x upregulated (2 times higher)
logFC = 0    →  No change
logFC = -1   →  2^-1 = 0.5x (down to half)
logFC = -2   →  2^-2 = 0.25x (down to quarter)
```

### P-value vs FDR

**P-value:** Probability that this result occurred by chance
- Raw p-value: Not corrected for multiple testing
- Can be misleading when testing thousands of genes

**FDR (False Discovery Rate):** Adjusted p-value
- Corrects for testing multiple genes
- More reliable
- **Use FDR < 0.05** (standard threshold)

### Volcano Plot Interpretation

**Upper Right Quadrant:** Significantly upregulated genes
- High fold change (right of center)
- Low FDR (high on y-axis)
- Target genes for investigation

**Upper Left Quadrant:** Significantly downregulated genes
- Negative fold change (left of center)
- Low FDR (high on y-axis)

**Bottom Center:** Non-significant genes
- Either small fold changes or high p-values
- Not of interest

### PCA Plot Interpretation

**Good PCA plot:**
- Samples from same condition cluster together
- Different conditions are well-separated
- No obvious outliers

**Problem PCA plot:**
- Conditions overlap significantly → weak signal
- One sample far from others → possible batch effect or contamination

---


