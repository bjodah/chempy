#!/bin/bash -xe
# Usage:
#
#    $ ./scripts/release.sh v1.2.3
#
if [[ $1 != v* ]]; then
    echo "Argument does not start with 'v'"
    exit 1
fi
./scripts/check_clean_repo_on_master.sh
cd $(dirname $0)/..
# PKG will be name of the directory one level up containing "__init__.py" 
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
PKG_UPPER=$(echo $PKG | tr '[:lower:]' '[:upper:]')
./scripts/run_tests.sh
env ${PKG_UPPER}_RELEASE_VERSION=$1 python setup.py sdist

# Build conda recipe:
./scripts/build_conda_recipe.sh $1

# All went well

git tag -a $1 -m $1
git push
git push --tags
twine upload dist/${PKG}-${1#v}.tar.gz
echo "Remember to bump version (and commit and push)!"
