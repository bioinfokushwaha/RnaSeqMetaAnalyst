FROM continuumio/miniconda3

ENV ENV_NAME=16_ppl_bioinfo-env
ENV DEBIAN_FRONTEND=noninteractive

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN conda env create -f /tmp/environment.yml

# Install DESeq2 and GenomeInfoDbData manually for latest version
SHELL ["conda", "run", "-n", "16_ppl_bioinfo-env", "/bin/bash", "-c"]
RUN Rscript -e "\
  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
  BiocManager::install(c('DESeq2', 'GenomeInfoDbData'), ask = FALSE, update = TRUE)"

# Set path and activate env on container start
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Copy entrypoint script
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
