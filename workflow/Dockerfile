# Use Ubuntu as the base image
FROM docker.io/continuumio/miniconda3:latest

# Update the package list and install some basic tools
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    csh \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy application code to the container
COPY environment.yml .
COPY ./dnatraj ./dnatraj

# Enable libmamba as solver
RUN conda config --set solver libmamba

# Create new environment
RUN conda env create --file=environment.yml && conda clean -afy

# Install mdplus
RUN conda run -n glimps_env pip install --upgrade pip setuptools wheel mdplus 

# Install dnatraj
RUN cd dnatraj && conda run -n glimps_env pip install . 

# Install Switch Perl module using conda-cpan
RUN conda run -n glimps_env cpan Switch 

# Run a command 
# Code to run when container is started
RUN echo '#!/bin/sh' > entrypoint.sh && \
    echo 'conda run --no-capture-output -n glimps_env "$@"' >> entrypoint.sh && \
    chmod +x entrypoint.sh

# Serve the app
ENTRYPOINT ["/app/entrypoint.sh"]

