FROM garcianacho/fhibasenanopore:31052022
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

COPY CommonFiles/ /home/docker/CommonFiles/
COPY Scripts/ /home/docker/Scripts/
COPY Binaries/ /home/docker/Binaries/
RUN mkdir /Inference \
    && mkdir /Noise \
    && mkdir -p /Models/RecombinantModel \
    && chmod +x /home/docker/Binaries/* \
    && chmod +x /home/docker/Scripts/* \
    && chmod -R 777 /home/docker/CommonFiles/ \
    && chmod 777 /Inference \
    && chmod 777 /Noise \
    && chmod 777 /Models/* \
    && mv /home/docker/Binaries/* /usr/bin/
USER docker
WORKDIR /home/docker/Fastq

