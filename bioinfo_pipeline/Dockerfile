FROM continuumio/miniconda3

LABEL maintainer="Sandeep Kumar Kushwaha"
LABEL description="RNA-Seq Meta Analysis from alignment to DEG (HISAT2/STAR + DESeq2/edgeR)"

WORKDIR /opt/bioinfo

COPY environment.yml .
RUN conda update -n base -c defaults conda -y && \
    conda env create -f environment.yml && \
    conda clean --all --yes

SHELL ["conda", "run", "--no-capture-output", "-n", "16_ppl_bioinfo-env", "/bin/bash", "-c"]

COPY . .

RUN chmod +x *.sh && chmod +x *.R

WORKDIR /opt/bioinfo/

CMD ["bash"]

