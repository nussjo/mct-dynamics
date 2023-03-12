FROM ubuntu:22.10

ARG DEBIAN_FRONTEND=noninteractive
ARG PIP_INDEX_URL

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        binutils \
        build-essential \
        curl \
        gcc \
        gcc-multilib \
        g++-multilib \
        gcc-12-multilib \
        g++-12-multilib \
        gdb \
        git \
        ninja-build \
        openjdk-11-jre-headless \
        python-is-python3 \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-tk \
        python3-venv \
        python3-wheel \
        ssh \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100

WORKDIR /opt/cuda

RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin && \
    mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600 && \
    wget https://developer.download.nvidia.com/compute/cuda/12.1.0/local_installers/cuda-repo-ubuntu2204-12-1-local_12.1.0-530.30.02-1_amd64.deb && \
    dpkg -i cuda-repo-ubuntu2204-12-1-local_12.1.0-530.30.02-1_amd64.deb && \
    cp /var/cuda-repo-ubuntu2204-12-1-local/cuda-*-keyring.gpg /usr/share/keyrings/ && \
    apt-get update && \
    apt-get -y install cuda

RUN useradd -ms /bin/bash nussjo

USER nussjo

RUN python -m pip install --upgrade \
        pip \
        cmake==3.25.0 \
        conan==1.57.0 \
        ninja==1.11.1

ENV PATH "/home/nussjo/.local/bin:$PATH"

WORKDIR /home/nussjo

RUN wget -L -O plantuml.jar https://github.com/plantuml/plantuml/releases/download/v1.2023.4/plantuml-1.2023.4.jar
