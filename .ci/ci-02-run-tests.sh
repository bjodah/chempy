#!/bin/bash
set -euxo pipefail
source /opt-3/cpython-v3.*-apt-deb/bin/activate
export CPATH=$SUNDBASE/include:${CPATH:-}
export LIBRARY_PATH=$SUNDBASE/lib
export LD_LIBRARY_PATH=$SUNDBASE/lib
PYTHON=python
INSTALL_PIP_FLAGS="--cache-dir $CACHE_ROOT/pip_cache ${INSTALL_PIP_FLAGS:-}"  # --user
$PYTHON -m pip install $INSTALL_PIP_FLAGS -e .[all]
./scripts/run_tests.sh --cov chempy --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
