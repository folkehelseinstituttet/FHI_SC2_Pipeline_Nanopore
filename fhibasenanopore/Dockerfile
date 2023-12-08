FROM garcianacho/fhibasenanopore:31052022
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

USER root

# Remove the existing 'artic' environment
RUN /bin/bash -c "conda env remove -n artic"

# Install mamba
RUN conda install mamba -n base -c conda-forge

# Switch to user 'docker'
USER docker

# Move the old 'fieldbioinformatics' directory, clone the new one, and create a new 'artic' environment using mamba
RUN cd /home/docker \
    && mv fieldbioinformatics/ fieldbioinformatics_old \
    && git clone https://github.com/artic-network/fieldbioinformatics \
    && cd fieldbioinformatics \
    && mamba env create -f environment.yml

# Activate the 'artic' environment and install the Python package
RUN /bin/bash -c ". activate artic && \
    cd /home/docker/fieldbioinformatics && \
    python setup.py install"

# Switch back to user 'root' if needed
USER root
