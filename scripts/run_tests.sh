#!/bin/bash -e
cd $(dirname $0)/..
python2 -m pytest --ignore build/ --ignore doc/ --ignore benchmarks/ --doctest-modules --pep8 --flakes $@
python3 -m pytest --ignore build/ --ignore doc/ --ignore benchmarks/
python -m doctest README.rst
