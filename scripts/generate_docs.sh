#!/bin/bash -xe
#
# Usage:
#
#    $ ./scripts/generate_docs.sh
#
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
AUTHOR=$(head -n 1 AUTHORS)
sphinx-apidoc --full --force -A "$AUTHOR" --module-first --doc-version=$(python setup.py --version) -F -o doc $PKG/ $PKG/tests
#sed -i 's/Contents/.. include:: ..\/README.rst\n\nContents/g' doc/index.rst
echo ".. include:: ../README.rst" >>doc/index.rst
sed -i "s/'sphinx.ext.viewcode',/'sphinx.ext.viewcode',\n    'numpydoc',/g" doc/conf.py
ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
( cd doc; PYTHONPATH=$ABS_REPO_PATH make html )
