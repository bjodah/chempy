#!/bin/bash -e
# Usage
#   $ ./scripts/run_tests.sh
# or
#   $ ./scripts/run_tests.sh --cov pycvodes --cov-report html
${PYTHON:-python} -m pytest --doctest-modules --pep8 --flakes $@
${PYTHON:-python} -m doctest README.rst
