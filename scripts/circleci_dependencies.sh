#!/bin/bash
if [[ ! -d $HOME/miniconda ]]; then
    curl -L --silent -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash miniconda.sh -b -p $HOME/miniconda
    rm miniconda.sh
    conda config --remove channels defaults
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --set show_channel_urls true
    conda config --set always_yes yes
    conda update --quiet --all
    conda install --quiet pip pytest numpy scipy matplotlib sympy quantities pyparsing
fi
