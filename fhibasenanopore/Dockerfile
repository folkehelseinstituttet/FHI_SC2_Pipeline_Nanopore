FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"
LABEL version="26042022"

RUN apt-get update \
    && apt-get install -y --no-install-recommends cmake
USER docker
RUN cd /home/docker \
    && wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/docker/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install ivar \
    && conda install -c bioconda seqkit \
    && conda update -n base -c defaults conda \
    && conda install -c bioconda bedtools \
    && conda install -c bioconda nextalign \
    && conda create -n pangolin \
    && cd /home/docker \
    && git clone https://github.com/artic-network/fieldbioinformatics \
    && cd fieldbioinformatics \
    && conda env create -f environment.yml 

USER root
RUN Rscript -e "install.packages('seqinr')"
USER docker    
RUN /bin/bash -c ". activate pangolin && \
    conda install -c bioconda -c conda-forge -c defaults pangolin && \
    conda update -y pangolin"


RUN /bin/bash -c ". activate artic && \
    cd /home/docker/fieldbioinformatics && \
    python setup.py install"

USER root
RUN Rscript -e "install.packages(c('doSNOW', \
'progress','foreach','parallel', 'pdftools', 'doParallel', \
'BiocManager', 'ggplot2', 'xgboost', 'reshape2', 'httr','phangorn','ggpubr' ,'ape','readxl', 'geomnet','rvest', 'tidyverse','writexl','nnet', 'stringr', 'devtools'))"
RUN Rscript -e "devtools::install_github('sctyner/geomnet')"
RUN Rscript -e "devtools::install_version('rvcheck',version='0.1.8')"
RUN rm -rf /var/lib/apt/lists/* \
    && rm /usr/bin/gcc /usr/bin/gcc-ar /usr/bin/gcc-nm /usr/bin/gcc-ranlib \
    && ln /usr/bin/gcc-nm-9 /usr/bin/gcc-nm \
    && ln /usr/bin/gcc-ar-9 /usr/bin/gcc-ar \
    && ln /usr/bin/gcc-9 /usr/bin/gcc \
    && ln /usr/bin/gcc-ranlib-9 /usr/bin/gcc-ranlib \
    && ln -s /home/docker/miniconda3/lib/libcrypto.so.1.1 /home/docker/miniconda3/lib/libcrypto.so.1.0.0
RUN Rscript -e "BiocManager::install(c('msa','GenomicAlignments','ggtree'))"
RUN Rscript -e "install.packages(c('lubridate', 'entropy', 'keras', 'abind'))"

USER docker
RUN pip install seaborn \
    && pip install matplotlib \
    && pip install pysam \
    && conda install nextflow \
    && Rscript -e "keras::install_keras()"

RUN conda create -n nextclade
RUN /bin/bash -c ". activate nextclade && \
    conda install -c bioconda nextclade"
    
USER root