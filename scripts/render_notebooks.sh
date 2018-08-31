#!/bin/bash -e
if [ -f index.ipynb ]; then
    sed -i.bak0 's/ipynb/html/' index.ipynb
    sed -i.bak1 's/filepath=index.html/filepath=index.ipynb/' index.ipynb  # mybinder link fix
fi
SCRIPTS_PATH=$(unset CDPATH && cd "$(dirname "$0")" && echo $PWD)
for dir in . examples/; do
    cd $dir
    find . -maxdepth 1 -iname "*.ipynb" | xargs -P 2 "$SCRIPTS_PATH/render_notebook.sh"
    cd -
done

cd examples/
../scripts/render_index.sh *.html
