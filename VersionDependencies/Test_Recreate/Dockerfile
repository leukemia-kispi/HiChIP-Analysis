# --------------------------------------------------------------
# Ubuntu 18.04 + Miniconda + Legacy Bioinformatics Environments
# --------------------------------------------------------------
FROM ubuntu:18.04
LABEL maintainer="Valdemar Priebe <valdemarpriebe@gmail.com>"

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    curl \
    ca-certificates \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda (Python 3.7)
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# Copy all Legacy environment YAMLs
COPY MACS2_Legacy.yml /envs/MACS2_Legacy.yml
COPY IDR_Legacy.yml /envs/IDR_Legacy.yml
COPY DovetailHiChIP_Legacy.yml /envs/DovetailHiChIP_Legacy.yml
COPY coolbox_Legacy.yml /envs/coolbox_Legacy.yml
COPY trimgalore_Legacy.yml /envs/trimgalore_Legacy.yml

# --------------------------------------------------------------
# Update conda base environment
# --------------------------------------------------------------
RUN conda install -y conda=4.8.5 && conda clean -afy

# --------------------------------------------------------------
# Create environments one by one (MACS2 last)
# --------------------------------------------------------------

# IDR_Legacy environment
RUN conda env create -f /envs/IDR_Legacy.yml && conda clean -afy

# DovetailHiChIP_Legacy environment
RUN conda env create -f /envs/DovetailHiChIP_Legacy.yml && conda clean -afy

# coolbox_Legacy environment
RUN conda env create -f /envs/coolbox_Legacy.yml && conda clean -afy

# --------------------------------------------------------------
# Install Java (needed for Picard)
# --------------------------------------------------------------
RUN apt-get update && apt-get install -y openjdk-8-jre-headless && rm -rf /var/lib/apt/lists/*
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64

# --------------------------------------------------------------
# Create Picard environment manually (stable version)
# --------------------------------------------------------------
RUN conda create -y -n Picard -c bioconda picard=2.25.7 && conda clean -afy

# trimgalore_Legacy environment
RUN conda env create -f /envs/trimgalore_Legacy.yml && conda clean -afy

# --------------------------------------------------------------
# Create MACS2 environment manually (bioconda macs2=2.2.6)
# --------------------------------------------------------------
RUN conda create -y -n MACS2 -c bioconda \
    macs2=2.2.6 \
    bowtie2=2.3.5.1 \
    samtools=1.7 \
    bedtools=2.30.0 \
    && conda clean -afy

# --------------------------------------------------------------
# Set bash as default shell
# --------------------------------------------------------------
SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash"]

