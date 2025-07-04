FROM continuumio/miniconda3

ENV ENV_NAME=16_ppl_bioinfo-env
ENV DEBIAN_FRONTEND=noninteractive

# Install mamba for faster environment resolution
RUN conda install -n base -c conda-forge mamba

# Copy environment definition
COPY environment.yml /tmp/environment.yml

# Create conda environment using mamba
RUN mamba env create -f /tmp/environment.yml

# Activate environment by default on shell startup
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Ensure path is set correctly
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Copy entrypoint
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
