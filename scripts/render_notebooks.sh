#!/bin/bash -e
cd examples/
jupyter nbconvert --debug --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb
../scripts/render_index.sh *.html
