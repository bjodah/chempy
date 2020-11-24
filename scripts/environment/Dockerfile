FROM debian:bullseye

MAINTAINER Björn Dahlgren <bjodah@gmail.com>

ENV LANG C.UTF-8

RUN apt-get update && \
    apt-get --quiet --assume-yes install curl git g++-10 gfortran-10 libgmp-dev binutils-dev bzip2 make cmake sudo \
    python3-dev python3-pip libboost-dev libgsl-dev liblapack-dev libsuitesparse-dev graphviz && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /tmp/sundials-5.5.0-build && \
    curl -Ls https://github.com/LLNL/sundials/releases/download/v5.5.0/sundials-5.5.0.tar.gz | tar xz -C /tmp && \
    FC=gfortran-10 cmake \
        -S /tmp/sundials-5.5.0 \
        -B /tmp/sundials-5.5.0-build \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DEXAMPLES_ENABLE_C=OFF \
        -DEXAMPLES_INSTALL=OFF \
        -DENABLE_LAPACK=ON \
        -DSUNDIALS_INDEX_SIZE=32 \
        -DENABLE_KLU=ON \
        -DKLU_INCLUDE_DIR=/usr/include/suitesparse \
        -DKLU_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu && \
    cmake --build /tmp/sundials-5.5.0-build && \
    cmake --build /tmp/sundials-5.5.0-build --target install && \
    rm -r /tmp/sundials-5.5.0*/ && \
    python3 -m pip install --upgrade-strategy=eager --upgrade pip && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# http://computation.llnl.gov/projects/sundials/download/sundials-5.5.0.tar.gz

# At this point the system should be able to pip-install the package and all of its dependencies. We'll do so
# when running the image using the ``host-jupyter-using-docker.sh`` script. Installed packages are cached.
