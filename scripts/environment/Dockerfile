FROM debian:stretch

MAINTAINER Björn Dahlgren <bjodah@gmail.com>

ENV LANG C.UTF-8

# This dockerfile is designed to run on binder (mybinder.org)
RUN apt-get update && \
    apt-get --quiet --assume-yes install curl git g++-6 libgmp-dev binutils-dev bzip2 make cmake sudo \
    python3-dev python3-pip libgsl-dev liblapack-dev graphviz && \
    apt-get clean && \
    curl -LOs http://computation.llnl.gov/projects/sundials/download/sundials-3.1.2.tar.gz && \
    tar xzf sundials-3.1.2.tar.gz && mkdir build/ && cd build/ && \
    cmake -DCMAKE_PREFIX_PATH=/usr/local -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF -DEXAMPLES_ENABLE=OFF -DEXAMPLES_INSTALL=OFF \
    ../sundials*/ && make install && cd - && rm -r build/ sundials* && \
    python3 -m pip install --upgrade pip && \
    curl -LOs http://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.bz2 && \
    tar xjf boost_*.tar.bz2 && cd boost* && ./bootstrap.sh && ./b2 -j 2 --prefix=/usr/local/ install && cd -

# At this point the system should be able to pip-install the package and all of its dependencies. We'll do so
# when running the image using the ``host-jupyter-using-docker.sh`` script. Installed packages are cached.
