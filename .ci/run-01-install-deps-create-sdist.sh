#!/bin/bash
set -euxo pipefail
. .ci/_source_common_env.sh
if [ ! -d $PYTHONUSERBASE ]; then mkdir -p $PYTHONUSERBASE; fi

python -m pip install --cache-dir $CACHE_ROOT/pip_cache --upgrade-strategy=eager --upgrade cython
python -m pip install --cache-dir $CACHE_ROOT/pip_cache -e .[all]
python -c "import pycvodes; import pyodesys; import pygslodeiv2"  # debug this CI config
git fetch -tq
#python setup.py sdist                    # test pip installable sdist (checks MANIFEST.in)
python -m build --sdist                    # test pip installable sdist (checks MANIFEST.in)
git archive -o dist/chempy-head.zip HEAD  # test pip installable zip (symlinks break)
