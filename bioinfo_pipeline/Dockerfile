FROM continuumio/miniconda3

ENV ENV_NAME=16_ppl_bioinfo-env
ENV DEBIAN_FRONTEND=noninteractive

COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml

SHELL ["conda", "run", "-n", "16_ppl_bioinfo-env", "/bin/bash", "-c"]

# Install latest BiocManager and DESeq2 + fix GenomeInfoDbData
RUN Rscript -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('DESeq2', 'GenomeInfoDbData'), ask = FALSE, update = TRUE)"

# Set environment as default on startup
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH
