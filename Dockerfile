FROM nvidia/cuda:12.4.1-cudnn-runtime-ubuntu22.04 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
  git \
  build-essential \
  python3 \
  python3-venv \
  python3-dev \
  cuda-nvcc-12-4 \
  cuda-libraries-dev-12-4 \
  && python3 -m venv /opt/venv \
  && . /opt/venv/bin/activate \
  && pip install --no-cache-dir --upgrade pip

# Copy the local source code
COPY . /app/boltz/

# Install from local source
WORKDIR /app/boltz
RUN . /opt/venv/bin/activate && pip install -e .[cuda]

# Clean up
RUN apt-get purge -y git \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

FROM nvidia/cuda:12.4.1-cudnn-runtime-ubuntu22.04

COPY --from=builder /opt/venv /opt/venv
COPY --from=builder /app/boltz /app/boltz
COPY --from=builder /usr/local/cuda-12.4 /usr/local/cuda-12.4

# Install build-essential and python3-dev for Triton/TriFast runtime compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
  python3 \
  python3-dev \
  build-essential \
  cuda-nvcc-12-4 \
  cuda-libraries-dev-12-4 \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PATH="/opt/venv/bin:/usr/local/cuda-12.4/bin:$PATH" \
  LD_LIBRARY_PATH="/usr/local/cuda-12.4/lib64:$LD_LIBRARY_PATH" \
  LANG=C.UTF-8 \
  PYTHONUNBUFFERED=1 \
  CC=gcc

WORKDIR /app

COPY LICENSE /app/LICENSE
COPY README.md /app/README.md
COPY examples /app/examples
COPY scripts /app/scripts
COPY docs /app/docs

ENTRYPOINT ["boltz"]
