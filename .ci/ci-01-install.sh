#!/bin/bash
set -euxo pipefail
export CPATH=$SUNDBASE/include:$CPATH
export LIBRARY_PATH=$SUNDBASE/lib
export LD_LIBRARY_PATH=$SUNDBASE/lib
source /opt-3/cpython-v3.*-apt-deb/bin/activate
python3 -m pip install --cache-dir $CACHE_ROOT/pip_cache --user -e .[all]
python3 -c "import pycvodes; import pyodesys; import pygslodeiv2"  # debug this CI config
git fetch -tq
python3 setup.py sdist                    # test pip installable sdist (checks MANIFEST.in)
git archive -o dist/chempy-head.zip HEAD  # test pip installable zip (symlinks break)
mkdir -p deploy/public_html/branches/${CI_COMMIT_BRANCH}
cp dist/chempy-* deploy/public_html/branches/${CI_COMMIT_BRANCH}/
