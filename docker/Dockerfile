FROM rocker/r-ver:4.4

MAINTAINER John Davis <cmap-soft@broadinstitute.org>

# Set environment variables to allow non-interactive installs
ENV DEBIAN_FRONTEND=noninteractive

# Switch to root user for system-wide installations
USER root

# 1. Install all system dependencies in a single layer
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    # R package dependencies
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libmariadb-dev-compat \
    libpq-dev \
    libssh2-1-dev \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    # Python build dependencies (Fixes the pyfarmhash error)
    python3-dev \
    build-essential \
    # General utilities
    jq \
    python3 \
    python3-pip \
    python3-venv \
    curl \
    # Clean up apt cache to reduce image size
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 2. Setup Python virtual environment and install Python packages
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --upgrade pip setuptools wheel \
    && /opt/venv/bin/pip install google-cloud-bigquery pyfarmhash pandas requests boto3 polars pyyaml

# 4. Make the virtual environment's python/pip the default
ENV PATH="/opt/venv/bin:${PATH}"

# 5. Copy R scripts into the container
COPY scripts/drc /scripts/drc
COPY scripts/biomarker /scripts/biomarker
COPY scripts/utils /scripts/utils
COPY docker/install_packages.R /src/install_packages.R

# 6. Install R dependencies
RUN Rscript /src/install_packages.R

# Label with git commit for build tracking
ARG GIT_COMMIT=unknown
LABEL org.opencontainers.image.revision=${GIT_COMMIT}