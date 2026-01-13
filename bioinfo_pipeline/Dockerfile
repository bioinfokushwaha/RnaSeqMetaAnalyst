FROM ubuntu:22.04

LABEL maintainer="RNA-Seq Meta-Analyst"
LABEL description="Unified container with three separate conda environments: CLI, R, and Python"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    bzip2 \
    bash \
    ca-certificates \
    locales \
    wget \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Set locale to prevent encoding issues
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8

# Install micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | \
    tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# Set micromamba environment variables
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
ENV SHELL=/bin/bash

# Set proper shell behavior for micromamba
SHELL ["/bin/bash", "-c"]

# Copy all three environment files
COPY environment-cli.yml /tmp/environment-cli.yml
COPY environment-r.yml /tmp/environment-r.yml
COPY environment-comparison.yml /tmp/environment-comparison.yml

# Create all three conda environments in one layer
RUN micromamba create -y -n rnaseq-cli -f /tmp/environment-cli.yml && \
    micromamba create -y -n rnaseq-r -f /tmp/environment-r.yml && \
    micromamba create -y -n rnaseq-comparison -f /tmp/environment-comparison.yml && \
    micromamba clean -a -y && \
    rm -rf /tmp/*.yml

# Force reinstall stringi for R environment (if needed)
RUN micromamba run -n rnaseq-r Rscript -e \
    "if (!require('stringi', quietly=TRUE)) install.packages('stringi', type='source', repos='https://cloud.r-project.org')" && \
    micromamba clean -a -y

# Create project structure
RUN mkdir -p /opt/project/scripts /data

# Copy ALL scripts to the container
COPY *.sh /opt/project/scripts/
COPY *.R /opt/project/scripts/
COPY *.py /opt/project/scripts/

# Make all shell scripts executable
RUN chmod +x /opt/project/scripts/*.sh

# Set working directory to /data (where volumes will be mounted)
WORKDIR /data

# Default command runs master.sh
CMD ["/bin/bash", "-c", "/opt/project/scripts/master.sh"]
