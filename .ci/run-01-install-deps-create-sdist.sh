#!/bin/bash
set -euxo pipefail
. .ci/_source_common_env.sh
if [ ! -d $PYTHONUSERBASE ]; then mkdir -p $PYTHONUSERBASE; fi

python3 -m pip install --cache-dir $CACHE_ROOT/pip_cache --user --upgrade-strategy=eager --upgrade cython
python3 -m pip install --cache-dir $CACHE_ROOT/pip_cache --user -e .[all]
python3 -c "import pycvodes; import pyodesys; import pygslodeiv2"  # debug this CI config
git fetch -tq
python3 setup.py sdist                    # test pip installable sdist (checks MANIFEST.in)
git archive -o dist/chempy-head.zip HEAD  # test pip installable zip (symlinks break)
