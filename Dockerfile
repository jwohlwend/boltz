# Boltz-1 Dockerfile
ARG BASE_IMAGE=nvcr.io/nvidia/pytorch:25.10-py3

FROM ${BASE_IMAGE} AS boltz-base

# Install core apt packages.
RUN --mount=type=cache,id=apt-cache,target=/var/cache/apt,sharing=locked \
  --mount=type=cache,id=apt-lib,target=/var/lib/apt,sharing=locked \
  bash -c '\
  apt-get update -qy && \
  apt-get install -qyy \
  libsndfile1 \
  ffmpeg \
  git \
  curl \
  pre-commit \
  lsof \
  git-lfs \
  sudo && \
  apt-get upgrade -qyy \
  rsync && \
  rm -rf /tmp/* /var/tmp/*'

RUN apt-get install -y gnupg

RUN mkdir -p /workspace/boltz/

# Fix for duplicate triton installation due to pytorch_triton renaming
RUN bash -c ' \
  cd /usr/local/lib/python3*/dist-packages/ && \
  PTRITON_DIR=$(ls -d pytorch_triton-*.dist-info 2>/dev/null | head -n 1) && \
  if [ -n "$PTRITON_DIR" ]; then \
    VERSION=${PTRITON_DIR#pytorch_triton-} && \
    VERSION=${VERSION%.dist-info} && \
    NEW_DIR="triton-${VERSION}.dist-info" && \
    cp -r "$PTRITON_DIR" "$NEW_DIR" && \
    sed -i "s/Name: pytorch-triton/Name: triton/" "$NEW_DIR/METADATA" && \
    echo "Successfully aliased pytorch-triton to triton version $VERSION"; \
  fi'

ENV NVIDIA_TF32_OVERRIDE=0

RUN bash -c 'echo "ubuntu ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/ubuntu && \
    chmod 0440 /etc/sudoers.d/ubuntu'

FROM boltz-base AS dev

# Install boltz-1
COPY ./README.md /workspace/boltz/README.md
COPY ./pyproject.toml /workspace/boltz/pyproject.toml
COPY ./src /workspace/boltz/src
COPY ./tests /workspace/boltz/tests
COPY ./scripts /workspace/boltz/scripts
COPY ./examples /workspace/boltz/examples

WORKDIR /workspace/boltz
RUN bash -c 'find . -name __pycache__ -type d -print | xargs rm -rf'
RUN pip install --no-build-isolation --editable .[lint,test,cuda,dev]

ENV NVIDIA_TF32_OVERRIDE=0
