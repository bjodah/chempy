#!/bin/bash
set -euxo pipefail
. .ci/_source_common_env.sh

./scripts/render_notebooks.sh
mv index.html index.ipynb.html
cp -r index.* examples/ "deploy/public_html/branches/${CI_COMMIT_BRANCH}"
