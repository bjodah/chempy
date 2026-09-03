#!/bin/bash
set -euxo pipefail
. .ci/_source_common_env.sh
if [ ! -d $PYTHONUSERBASE ]; then mkdir -p $PYTHONUSERBASE; fi

python -m pip install --cache-dir $CACHE_ROOT/pip_cache --upgrade-strategy=eager --upgrade cython

# Build pykinsol against $PYKINSOL_SUNDIALS_ROOT (see _source_common_env.sh)
# and publish the wheel via $PIP_FIND_LINKS so every later "pip install" that
# resolves "pykinsol>=0.1.6" (below, and the sdist/git-archive smoke installs
# into fresh virtualenvs in run-02) picks it up pre-built instead of trying to
# compile it against the newer default $SUNDIALS_ROOT.
mkdir -p $PIP_FIND_LINKS
(
    export SUNDIALS_ROOT=$PYKINSOL_SUNDIALS_ROOT
    export CPATH=$SUNDIALS_ROOT/include:$SUITESPARSE_ROOT/include/suitesparse
    export LIBRARY_PATH=$SUNDIALS_ROOT/lib:$SUITESPARSE_ROOT/lib
    python -m pip wheel --no-deps --wheel-dir $PIP_FIND_LINKS --cache-dir $CACHE_ROOT/pip_cache "pykinsol>=0.1.6"
)

# Extract "anyode/anyode.hpp" et al. (see $ANYODE_INCLUDE_ROOT in
# _source_common_env.sh) from pycvodes' sdist, since the installed pycvodes
# package doesn't ship the vendored anyode headers it was built against, but
# pyodesys' native codegen needs them at test time.
if [ ! -d $ANYODE_INCLUDE_ROOT/anyode ]; then
    mkdir -p $ANYODE_INCLUDE_ROOT
    tmp_anyode_sdist=$(mktemp -d)
    python -m pip download --no-deps --no-binary pycvodes --cache-dir $CACHE_ROOT/pip_cache --dest $tmp_anyode_sdist "pycvodes>=0.14.5"
    tar xf $tmp_anyode_sdist/pycvodes-*.tar.gz -C $tmp_anyode_sdist
    cp -r $tmp_anyode_sdist/pycvodes-*/external/anyode/include/anyode $ANYODE_INCLUDE_ROOT/
    rm -rf $tmp_anyode_sdist
fi

python -m pip install --cache-dir $CACHE_ROOT/pip_cache -e .[all]
python -c "import pycvodes; import pyodesys; import pygslodeiv2"  # debug this CI config
git fetch -tq
#python setup.py sdist                    # test pip installable sdist (checks MANIFEST.in)
python -m build --sdist                    # test pip installable sdist (checks MANIFEST.in)
git archive -o dist/chempy-head.zip HEAD  # test pip installable zip (symlinks break)
