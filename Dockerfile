FROM continuumio/miniconda3:latest

ARG CMAKE_VERSION=3.29.2
ARG CMAKE_OS=linux
ARG CMAKE_ARCH=x86_64

RUN apt-get update && apt-get install -y \
    build-essential \
    ninja-build \
    wget \
    tar \
    && wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-${CMAKE_OS}-${CMAKE_ARCH}.tar.gz \
    && tar -xzvf cmake-${CMAKE_VERSION}-${CMAKE_OS}-${CMAKE_ARCH}.tar.gz \
    && mv cmake-${CMAKE_VERSION}-${CMAKE_OS}-${CMAKE_ARCH} /opt/cmake \
    && ln -s /opt/cmake/bin/* /usr/local/bin/ \
    && cmake --version \
    && rm cmake-${CMAKE_VERSION}-${CMAKE_OS}-${CMAKE_ARCH}.tar.gz \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y docker.io

WORKDIR /osmordred