FROM nvidia/cuda:12.4.0-runtime-ubuntu22.04

LABEL Author="ZipBio" Version="1.0"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
    apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get update && \
    apt-get install -y python3.11 python3-pip && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3 1 && \
    pip install --upgrade pip && \
    pip install boltz -U
