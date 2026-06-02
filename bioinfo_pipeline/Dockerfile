FROM ubuntu:22.04

LABEL maintainer="RNA-Seq Meta-Analyst"
LABEL description="Unified container: rnaseq-cli, rnaseq-r, rnaseq-comparison, rnaseq-busco, rnaseq-busco-py"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
# Added libssl-dev and zlib1g-dev needed by BUSCO / Augustus at runtime
RUN apt-get update && apt-get install -y \
    curl \
    bzip2 \
    bash \
    ca-certificates \
    locales \
    wget \
    gzip \
    libssl-dev \
    zlib1g-dev \
    libncurses5-dev \
    libboost-all-dev \
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

# ===================================================================
# LAYER 1-3: Original three environments (CLI, R, Python comparison)
# ===================================================================
COPY environment-cli.yml        /tmp/environment-cli.yml
COPY environment-r.yml          /tmp/environment-r.yml
COPY environment-comparison.yml /tmp/environment-comparison.yml

RUN micromamba create -y -n rnaseq-cli        -f /tmp/environment-cli.yml        && \
    micromamba create -y -n rnaseq-r          -f /tmp/environment-r.yml          && \
    micromamba create -y -n rnaseq-comparison -f /tmp/environment-comparison.yml && \
    micromamba clean -a -y && \
    rm -rf /tmp/environment-cli.yml /tmp/environment-r.yml /tmp/environment-comparison.yml

# Force reinstall stringi for R environment (robustness fix)
RUN micromamba run -n rnaseq-r Rscript -e \
    "if (!require('stringi', quietly=TRUE)) install.packages('stringi', type='source', repos='https://cloud.r-project.org')" && \
    micromamba clean -a -y

# ===================================================================
# LAYER 4: rnaseq-busco
# Purpose : Run BUSCO 6.0.0 + SEPP 4.5.5 on protein/genome/RNA input
# Python  : 3.9 (pinned — BUSCO 6 incompatible with Python 3.11)
# pandas  : 1.5.3 (pinned — BUSCO 6 incompatible with pandas 2.x)
# Isolated from rnaseq-busco-py to avoid solver conflicts
# Usage   :
#   micromamba run -n rnaseq-busco \
#     busco -i protein.faa -l mammalia_odb12 -o busco_output \
#           -m proteins -f -c 16
#   micromamba run -n rnaseq-busco \
#     awk '$2=="Complete"{printf $3"\t"; for(i=7;i<=NF;i++){printf $i" "} print ""}' \
#       full_table.tsv > busco_core_genes1.txt
# ===================================================================
COPY environment-busco.yml /tmp/environment-busco.yml

RUN micromamba create -y -n rnaseq-busco -f /tmp/environment-busco.yml && \
    micromamba clean -a -y && \
    rm -f /tmp/environment-busco.yml

# ===================================================================
# LAYER 5: rnaseq-busco-py
# Purpose : Run the GFF/GTF BUSCO-filtering Python script
#           (uses gffpandas, pandas 2.x, Python 3.11)
# Kept separate from rnaseq-busco because:
#   - BUSCO 6 requires pandas <=1.5.x; gffpandas works best with pandas 2.x
#   - Python 3.11 needed for gffpandas 1.2.1 full feature support
# Usage   :
#   micromamba run -n rnaseq-busco-py \
#     python /opt/project/scripts/busco_filter.py
# ===================================================================
COPY environment-busco-py.yml /tmp/environment-busco-py.yml

RUN micromamba create -y -n rnaseq-busco-py -f /tmp/environment-busco-py.yml && \
    micromamba clean -a -y && \
    rm -f /tmp/environment-busco-py.yml

# ===================================================================
# PROJECT STRUCTURE & SCRIPTS
# ===================================================================
# Includes: master.sh, master_busco.sh, busco_filter.py,
#           bowtie.sh, hisat2.sh, star.sh, *_clean.sh,
#           run_deseq2_analysis.R, run_edgeR_analysis.R,
#           run_pipeline_comparison.py, and all wrappers
# NOTE: run_busco_pipeline.sh is NOT used — master_busco.sh
#       is the correct runtime entrypoint for BUSCO analysis.
RUN mkdir -p /opt/project/scripts /data

COPY *.sh  /opt/project/scripts/
COPY *.R   /opt/project/scripts/
COPY *.py  /opt/project/scripts/

RUN chmod +x /opt/project/scripts/*.sh

# ===================================================================
# PERMISSIONS: allow non-root execution
# ===================================================================
RUN mkdir -p /.cache/mamba /.cache/pip /opt/conda/pkgs && \
    chmod -R 777 /.cache       && \
    chmod -R 777 /opt/conda    && \
    chmod -R 777 /data         && \
    chmod -R 777 /opt/project

WORKDIR /data

# ===================================================================
# DEFAULT COMMAND — runs main RNA-seq pipeline.
# Override at runtime for BUSCO:
#
#   docker run --rm \
#     -v $(pwd)/results:/data \
#     -v $(pwd)/data/genome:/data/genome:ro \
#     -e PROTEIN_FILE=/data/genome/protein.faa \
#     -e BUSCO_LINEAGE=mammalia_odb12 \
#     -e GFF_FILE=/data/genome/genomic.gff \
#     -e GTF_FILE=/data/genome/genomic.gtf \
#     -e THREADS=16 \
#     rnaseqmetaanalyst:latest \
#     /bin/bash -c "/opt/project/scripts/master_busco.sh"
# ===================================================================
CMD ["/bin/bash", "-c", "/opt/project/scripts/master.sh"]
