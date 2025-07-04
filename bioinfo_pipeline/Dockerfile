FROM continuumio/miniconda3

ENV ENV_NAME=16_ppl_bioinfo-env
ENV DEBIAN_FRONTEND=noninteractive

COPY environment.yml /tmp/environment.yml

# Create the conda environment
RUN conda env create -f /tmp/environment.yml

# Activate environment and install DESeq2 + GenomeInfoDbData from R
SHELL ["conda", "run", "-n", "16_ppl_bioinfo-env", "/bin/bash", "-c"]

RUN Rscript -e "\
  install.packages('BiocManager', repos='https://cloud.r-project.org'); \
  BiocManager::install(c('DESeq2', 'GenomeInfoDbData'), ask = FALSE, update = TRUE)"

# Activate environment by default on shell startup
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Ensure path is set correctly
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Copy entrypoint
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
