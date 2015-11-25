#!/bin/bash -xe
#
# Usage:
#
#    $ ./scripts/generate_docs.sh
#
# Usage if doc/ is actually published in master branch on github:
#
#    $ ./scripts/generate_docs.sh my_github_username my_github_repo master
#
NARGS=$#
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
AUTHOR=$(head -n 1 AUTHORS)
sphinx-apidoc --full --force -A "$AUTHOR" --module-first --doc-version=$(python setup.py --version) -F -o doc $PKG/ $PKG/tests/
#sed -i 's/Contents/.. include:: ..\/README.rst\n\nContents/g' doc/index.rst
echo ".. include:: ../README.rst" >>doc/index.rst
sed -i "s/'sphinx.ext.viewcode',/'sphinx.ext.viewcode',\n    'sphinx.ext.autosummary',\n    'numpydoc',/g" doc/conf.py
sed -i "s/alabaster/sphinx_rtd_theme/g" doc/conf.py
if [[ $NARGS -eq 3 ]]; then
cat <<EOF>>doc/conf.py
context = {
    'conf_py_path': '/doc/',
    'github_user': '$1',
    'github_repo': '$2',
    'github_version': '$3',
    'display_github': True,
    'source_suffix': '.rst',
}

if 'html_context' in globals():
    html_context.update(context)
else:
    html_context = context
EOF
fi
echo "numpydoc_class_members_toctree = False" >>doc/conf.py
ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
( cd doc; PYTHONPATH=$ABS_REPO_PATH make html )
