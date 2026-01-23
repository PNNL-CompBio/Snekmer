# 1) Small base image with Python 3.11
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
      git bash build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# (Optional) let you switch branches at build time:
ARG SNEKMER_BRANCH=main

RUN git clone https://github.com/PNNL-CompBio/Snekmer.git \
 && cd Snekmer \
 && git checkout ${SNEKMER_BRANCH} \
 \
 # create the venv
 && python3.11 -m venv /opt/snekmer_env \
 \
 && /opt/snekmer_env/bin/pip install --upgrade pip setuptools wheel \
 \
 # install deps and the package
 && /opt/snekmer_env/bin/pip install -r requirements.txt \
 && /opt/snekmer_env/bin/pip install .

# 4) Entrypoint
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENV PATH="/opt/snekmer_env/bin:${PATH}"
WORKDIR /work

ENTRYPOINT ["entrypoint.sh"]

