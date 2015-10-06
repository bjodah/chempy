#!/bin/bash -e
cd $(dirname $0)/..
python2 -m pytest --ignore build/ --ignore doc/ --doctest-modules --pep8 --flakes $@
python3 -m pytest --ignore build/ --ignore doc/
python -m doctest README.rst
