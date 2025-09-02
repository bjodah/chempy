#!/bin/bash
set -euxo pipefail
. .ci/_source_common_env.sh

bash -c '[[ $(python3 setup.py --version) =~ ^[0-9]+.* ]]'
./scripts/run_tests.sh --cov chempy --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
cp -r htmlcov/ dist/*.* deploy/public_html/branches/${CI_COMMIT_BRANCH}/
./.ci/grep-for-merge-blocking-token.sh
export CHEMPY_DEPRECATION_FILTER=ignore
python3 -m virtualenv /tmp/test_sdist
python3 -m virtualenv /tmp/test_git_archive
cd deploy/public_html/branches/${CI_COMMIT_BRANCH}
unset CHEMPY_SKIP_NO_TESTS  # I can't get pip to install extras when using local file...
bash -c "\
source /tmp/test_sdist/bin/activate \
; pip install \
    --cache-dir $CACHE_ROOT/pip_cache \
    file://$(realpath $(eval ls chempy-*.tar.gz))#chempy[all] pytest \
; pytest --pyargs chempy"

bash -c "\
source /tmp/test_git_archive/bin/activate\
; pip install \
    --cache-dir $CACHE_ROOT/pip_cache \
    file://$(realpath chempy-head.zip)#chempy[all] pytest\
; pytest --pyargs chempy"
