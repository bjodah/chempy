#!/bin/bash
echo "Running chempy tests with environment variables for the compiler pointing to prefix ..."
set -x
CPLUS_INCLUDE_PATH=$PREFIX/include LIBRARY_PATH=$PREFIX/lib LD_LIBRARY_PATH=$PREFIX/lib MPLBACKEND=agg LLAPACK=openblas python -m pytest --verbose -ra --pyargs chempy
