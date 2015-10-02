#!/bin/bash -ex
python2 -m pytest --pep8 --flakes --ignore doc/ --ignore examples/aqchem
python3 -m pytest --ignore build/ --ignore doc/ --ignore examples/aqchem
python -m doctest README.rst
