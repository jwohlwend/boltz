FROM nvidia/cuda:12.6.3-base-ubuntu24.04

# Install required packages
RUN apt-get update && apt-get install -y python3 python3-pip python3-venv

# Copy application files
COPY . /app
WORKDIR /app

# Create and activate virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install dependencies in virtual environment
RUN pip install -e .[cuda]

CMD ["/bin/bash"]