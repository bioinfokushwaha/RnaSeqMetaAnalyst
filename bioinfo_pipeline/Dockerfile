# Base image with Miniconda
FROM continuumio/miniconda3:latest

# Set environment variables
ENV CONDA_ENV_NAME=16_ppl_bioinfo-env
ENV DEBIAN_FRONTEND=noninteractive

# Copy environment.yml into the container
COPY environment.yml /tmp/environment.yml

# Install system packages required by bioinformatics tools
RUN apt-get update && apt-get install -y \
    build-essential \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    wget \
    curl \
    git \
    unzip \
    libncurses5-dev \
    libncursesw5-dev \
    locales \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Create conda environment from environment.yml
RUN SSL_NO_VERIFY=1 conda env create -f /tmp/environment.yml && conda clean -afy

# Activate environment by default
SHELL ["conda", "run", "--no-capture-output", "-n", "16_ppl_bioinfo-env", "/bin/bash", "-c"]

# Set working directory
WORKDIR /pipeline

# Entrypoint into conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "16_ppl_bioinfo-env", "/bin/bash"]

