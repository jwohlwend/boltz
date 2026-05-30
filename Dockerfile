FROM pytorch/pytorch:2.7.0-cuda12.6-cudnn9-runtime

LABEL author="Colby T. Ford <colby@tuple.xyz>"

## Environment Settings
ENV DEBIAN_FRONTEND=noninteractive

## Install Basic Dependencies
RUN apt-get clean && \
    apt-get update && \
    apt-get -y install \
        sudo \
        git \
        curl \
        wget \
        g++

## Install Boltz (+ dependencies)
# RUN pip install boltz==2.1.1 -U
RUN pip install boltz[cuda]==2.2.1 -U

## Download CCD data and model weights 
RUN mkdir ~/.boltz && \
    cd ~/.boltz && \
    wget https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar && \
    tar -xvf mols.tar && \
    rm mols.tar && \
    wget https://model-gateway.boltz.bio/boltz2_conf.ckpt && \
    wget https://model-gateway.boltz.bio/boltz2_aff.ckpt


CMD ["boltz"]