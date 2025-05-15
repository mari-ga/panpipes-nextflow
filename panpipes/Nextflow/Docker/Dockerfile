FROM continuumio/miniconda3
# the following ARG conda_env variable must match the conda env name defined in environment.yml:
ARG conda_env=preprocessing

LABEL authors="M. Gonzales A." \
  maintainer="M.Gonzales <mariana.gonzales-andre@tum.de>" \
  description="Docker image for preprocessing modules - panpipes"

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf /var/cache/apt/*

# Install system dependencies - combining apt commands to reduce layers
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
    python3.11 \
    build-essential curl git procps \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Update conda
RUN conda update -n base -c defaults conda

# Add environment.yml file
COPY preprocessing.yml /tmp/preprocessing.yml

# Create conda environment
RUN conda env create -f /tmp/preprocessing.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "preprocessing", "/bin/bash", "-c"]

# Install all required Python packages
#RUN pip install scanpy pandas muon anndata spatialdata matplotlib seaborn glob2 \
#    && pip install argparse logging warnings

# Set PATH to use the conda environment by default
ENV PATH /opt/conda/envs/${conda_env}/bin:$PATH
ENV CONDA_DEFAULT_ENV ${conda_env}

# Activate the environment when container starts
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bff"]
CMD ["/bin/bash"]
