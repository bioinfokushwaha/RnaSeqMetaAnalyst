# Base image with Miniconda
FROM continuumio/miniconda3:23.3.1-0

# Set environment variables
ENV CONDA_ENV_NAME=16_ppl_bioinfo-env
ENV DEBIAN_FRONTEND=noninteractive

# Set locale
RUN apt-get update && apt-get install -y \
    locales && \
    locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN SSL_NO_VERIFY=1 conda env create -f /tmp/environment.yml && \
    conda clean -afy && \
    echo "conda activate $CONDA_ENV_NAME" >> /root/.bashrc

# Install system packages
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
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Optional: reduce Conda footprint
RUN find /opt/conda/ -name '*.a' -delete \
 && find /opt/conda/ -name '__pycache__' -type d -exec rm -r {} + \
 && rm -rf /opt/conda/pkgs/*

# Set working directory
WORKDIR /pipeline

# Set shell to use the conda environment
SHELL ["conda", "run", "--no-capture-output", "-n", "16_ppl_bioinfo-env", "/bin/bash", "-c"]

# Default entry
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "16_ppl_bioinfo-env", "/bin/bash"]
CMD ["bash"]
