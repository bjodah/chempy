#!/bin/bash -e
# Usage
#   $ ./scripts/run_tests.sh
# or
#   $ ./scripts/run_tests.sh --cov chempy --cov-report html
CHEMPY_DEPRECATION_FILTER='ignore' ${PYTHON:-python3} -m pytest -ra --doctest-modules --flakes $@
MPLBACKEND=Agg ${PYTHON:-python3} -m doctest README.rst
rstcheck README.rst
