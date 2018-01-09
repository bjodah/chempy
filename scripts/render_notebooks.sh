#!/bin/bash -e
for dir in . examples/; do
    cd $dir
    jupyter nbconvert --debug --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb
    cd -
done
#../scripts/render_index.sh *.html

cd examples/
../scripts/render_index.sh *.html
