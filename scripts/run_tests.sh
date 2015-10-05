#!/bin/bash -ex
python2 -m pytest --pep8 --flakes --ignore doc/ $@
python3 -m pytest --ignore doc/ $@
python -m doctest README.rst
