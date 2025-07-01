FROM continuumio/miniconda3

LABEL maintainer="Sandeep Kumar Kushwaha"
LABEL description="RNA-Seq Meta Analysis from alignment to DEG (HISAT2/STAR + DESeq2/edgeR)"

# Set working directory
WORKDIR /opt/bioinfo

# Install conda environment
COPY environment.yml .
RUN conda update -n base -c defaults conda -y && \
    conda env create -f environment.yml && \
    conda clean --all --yes

# Set shell to default bash (we'll handle env manually via entrypoint)
SHELL ["/bin/bash", "-c"]

# Copy pipeline scripts
COPY . .

# Make scripts executable
RUN chmod +x *.sh *.R

# Add entrypoint to run everything inside the conda environment
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

ENTRYPOINT ["./entrypoint.sh"]
CMD ["bash"]

