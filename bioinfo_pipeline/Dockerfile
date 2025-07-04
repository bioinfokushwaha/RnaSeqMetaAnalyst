FROM bioconductor/bioconductor_docker:RELEASE_3_16

ENV DEBIAN_FRONTEND=noninteractive

# Install CLI tools from Bioconda
RUN apt-get update && apt-get install -y curl git vim && apt-get clean

# Install bioinformatics tools and R packages with conda
RUN conda install -y -c bioconda -c conda-forge \
    fastqc=0.11.9 \
    fastp=0.23.2 \
    multiqc=1.14 \
    sortmerna=2.1b \
    hisat2=2.2.1 \
    star=2.7.10b \
    bowtie2=2.4.5 \
    subread=2.0.6 \
    htseq=0.13.5 \
    rsem=1.3.3 \
    samtools=1.15.1 \
    r-dplyr \
    r-openxlsx \
    r-tidyverse \
    && conda clean -afy

# Install Bioconductor packages using R
RUN Rscript -e "BiocManager::install(c( \
    'DESeq2', 'edgeR', 'GenomeInfoDb', 'Biobase', \
    'S4Vectors', 'SummarizedExperiment', 'BiocGenerics', 'MatrixGenerics' \
))"

# Set default command
WORKDIR /data
CMD ["R"]
