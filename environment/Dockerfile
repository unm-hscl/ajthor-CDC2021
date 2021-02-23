# hash:sha256:0c92695f775c217e702d2124a6e32850b6878121292fc5f2d12b2f1a189fa70d
FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.4ubuntu1 \
        bzip2=1.0.6-8.1ubuntu0.2 \
        m4=1.4.18-1 \
    && rm -rf /var/lib/apt/lists/*

# Matlab installation script.
COPY setup_dependencies.m /

# post_install script.
COPY post_install.sh /
RUN /post_install.sh
