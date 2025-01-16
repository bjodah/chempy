#!/bin/bash
set -euxo pipefail
export CPATH=$SUNDBASE/include:${CPATH:-}
export LIBRARY_PATH=$SUNDBASE/lib
export LD_LIBRARY_PATH=$SUNDBASE/lib
source /opt-3/cpython-v3.*-apt-deb/bin/activate
PYTHON=python
INSTALL_PIP_FLAGS="--cache-dir $CACHE_ROOT/pip_cache ${INSTALL_PIP_FLAGS:-}"  # --user
for pypkg in pyodeint pygslodeiv2 pycompilation pycodeexport pycvodes pykinsol sym pyodesys; do
    case $pypkg in
        sym)
            pypkg_fqn="git+https://github.com/bjodah/sym@jun21#egg=sym"
            ;;
        pyodeint)
            pypkg_fqn="git+https://github.com/bjodah/pyodeint@sep21#egg=pyodeint"
            ;;
        pygslodeiv2)
            pypkg_fqn="git+https://github.com/bjodah/pygslodeiv2@cython-except-plus#egg=pygslodeiv2"
            ;;
        pycompilation)
            pypkg_fqn="git+https://github.com/bjodah/pycompilation@use-importlib-rather-than-imp#egg=pycompilation"
            ;;
        pycodeexport)
            pypkg_fqn="git+https://github.com/bjodah/pycodeexport@qulify-extension-name-and-new-ci#egg=pycodeexport"
            ;;
        pycvodes)
            pypkg_fqn="git+https://github.com/bjodah/pycvodes@may21#egg=pycvodes"
            ;;
        pyodesys)
            pypkg_fqn="git+https://github.com/bjodah/pyodesys@bdf2#egg=pyodesys"
            ;;
        pykinsol)
            pypkg_fqn="git+https://github.com/bjodah/pykinsol@jan25#egg=pykinsol"
            ;;
        *)
            pypkg_fqn=$pypkg
            ;;
    esac
    $PYTHON -m pip install $INSTALL_PIP_FLAGS $pypkg_fqn
    #( cd /tmp; $PYTHON -m pytest -k "not pool_discontinuity_approx" --pyargs $pypkg )
done

$PYTHON -m pip install $INSTALL_PIP_FLAGS -e .[all]
$PYTHON -c "import pycvodes; import pyodesys; import pygslodeiv2"
#git fetch -tq
$PYTHON setup.py sdist                    # test pip installable sdist (checks MANIFEST.in)
git archive -o dist/chempy-head.zip HEAD  # test pip installable zip (symlinks break)
mkdir -p deploy/public_html/branches/${CI_COMMIT_BRANCH}
cp dist/chempy-* deploy/public_html/branches/${CI_COMMIT_BRANCH}/

set +e
[[ $($PYTHON setup.py --version) =~ ^[0-9]+.* ]]
./scripts/run_tests.sh --cov chempy --cov-report html
bash
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg

