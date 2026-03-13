# Use Python 3.10 base image
FROM python:3.10-slim

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    pkg-config \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install mwa_hyperbeam
RUN pip install --no-cache-dir mwa_hyperbeam==0.10.4

WORKDIR /app

COPY hyperbeam_array.py /app
COPY hyperbeam_single.py /app

