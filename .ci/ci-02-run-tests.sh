#!/bin/bash
set -euxo pipefail
source /opt-3/cpython-v3.*-apt-deb/bin/activate
[[ $(python3 setup.py --version) =~ ^[0-9]+.* ]]
./scripts/run_tests.sh --cov chempy --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
