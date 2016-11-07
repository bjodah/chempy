#!/bin/bash -xeu
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi
git archive -o /tmp/$PKG_NAME.zip HEAD  # test pip installable zip (symlinks break)
python3 -m pip install /tmp/$PKG_NAME.zip

python2 setup.py sdist  # test pip installable sdist (checks MANIFEST.in)
(cd dist/; python2 -m pip install $PKG_NAME-$(python2 ../setup.py --version).tar.gz)
(cd /; python2 -m pytest --pyargs $PKG_NAME)

PYTHON=python2 ./scripts/run_tests.sh
PYTHON=python3 ./scripts/run_tests.sh --cov $PKG_NAME --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
! grep "DO-NOT-MERGE!" -R . --exclude ci.sh

# Test package without any 3rd party libraries (only python stdlib):
pip install virtualenv
python -m virtualenv venv
set +u
(source ./venv/bin/activate; python -m pip install pytest; python -m pytest chempy)

python2 -m pip install .[all]
python3 -m pip install .[all]

./scripts/render_notebooks.sh
./scripts/generate_docs.sh
