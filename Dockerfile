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
    # Python build dependencies
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

# 2. Setup Python virtual environment
RUN python3 -m venv /opt/venv
# Make the virtual environment's python/pip the default
ENV PATH="/opt/venv/bin:${PATH}"
# Upgrade pip inside the venv
RUN pip install --upgrade pip setuptools wheel

# 3. Install the local sushi-tools Python package and all its dependencies
WORKDIR /app
COPY pyproject.toml README.md ./
COPY sushilib ./sushilib  # Copy the new, clean source directory
RUN pip install -e .      # Install the package and all dependencies from pyproject.toml

# 4. Copy R scripts and other project scripts into the container
# This copies all your script folders (drc, biomarker, etc.)
COPY scripts /scripts
COPY docker/install_packages.R /src/install_packages.R

# 5. Install R dependencies
RUN Rscript /src/install_packages.R

