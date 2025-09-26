# Base image
FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------
# Install system packages including Java 17
# -----------------------------
RUN apt-get update && apt-get install -y \
    wget curl unzip bzip2 ca-certificates git \
    build-essential zlib1g-dev libbz2-dev liblzma-dev \
    fastqc fastp multiqc \
    kraken2 bowtie2 \
    openjdk-17-jdk \
    python3 python3-pip \
    r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------
# Install Miniforge3 (Conda)
# -----------------------------
RUN wget https://github.com/conda-forge/miniforge/releases/download/24.11.0-0/Miniforge3-24.11.0-0-Linux-x86_64.sh -O /tmp/miniforge3.sh && \
    bash /tmp/miniforge3.sh -b -p /opt/conda && \
    rm /tmp/miniforge3.sh

# -----------------------------
# Set JAVA_HOME and PATH so Java 17 is detected first
# -----------------------------
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:/opt/conda/bin:$PATH

# Verify Java
RUN java -version

# -----------------------------
# Configure Conda channels
# -----------------------------
RUN conda config --system --add channels conda-forge && \
    conda config --system --add channels bioconda && \
    conda config --system --set channel_priority flexible

# -----------------------------
# Install bioinformatics tools with pinned versions
# -----------------------------
RUN conda install -c bioconda -c conda-forge -y \
    humann=3.9 \
    kneaddata=0.12.3 \
    diamond=2.1.10 \
   # metaphlan=3.0.14 \
    && conda clean -afy

# -----------------------------
# Install Python packages
# -----------------------------
RUN pip3 install --no-cache-dir \
    pandas numpy scipy biopython matplotlib seaborn

# -----------------------------
# Install R packages
# -----------------------------
RUN Rscript -e "install.packages(c('tidyverse','data.table','ggplot2'), repos='https://cloud.r-project.org')"

# -----------------------------
# Install Nextflow
# -----------------------------
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && chmod +x /usr/local/bin/nextflow

# Verify Nextflow
RUN nextflow -v


# -----------------------------
# Set working directory
# -----------------------------
WORKDIR /work

# -----------------------------
# Default command
# -----------------------------
CMD ["/bin/bash"]

