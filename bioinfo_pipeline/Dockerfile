FROM ubuntu:22.04

LABEL maintainer="you@example.com"

# Install micromamba
RUN apt-get update && apt-get install -y curl bzip2 bash ca-certificates && \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
ENV SHELL=/bin/bash

# Set proper shell behavior for micromamba
SHELL ["/bin/bash", "-c"]

# Copy environment files
COPY environment-cli.yml /tmp/environment-cli.yml
COPY environment-r.yml /tmp/environment-r.yml

# Create both environments
RUN micromamba create -y -n rnaseq-cli -f /tmp/environment-cli.yml && \
    micromamba create -y -n rnaseq-r -f /tmp/environment-r.yml && \
    micromamba clean -a -y


# Copy all .sh and .R files from current directory into the image
COPY *.sh /opt/project/scripts/
COPY *.R /opt/project/scripts/

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"]

