#!/bin/bash -e
if [ -f index.ipynb ]; then
    sed -i.bak0 's/ipynb/html/' index.ipynb
    sed -i.bak1 's/filepath=index.html/filepath=index.ipynb/' index.ipynb  # mybinder link fix
fi
for dir in . examples/; do
    cd $dir
    jupyter nbconvert --debug --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb
    cd -
done

cd examples/
../scripts/render_index.sh *.html
