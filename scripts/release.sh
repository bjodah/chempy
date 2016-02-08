#!/bin/bash -xeu
# Usage:
#
#    $ ./scripts/release.sh v1.2.3 ~/anaconda2/bin
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
env ${PKG_UPPER}_RELEASE_VERSION=$1 ./scripts/generate_docs.sh
PATH=$2:$PATH ./scripts/build_conda_recipe.sh $1
# All went well
git tag -a $1 -m $1
git push
git push --tags
VERSION=${1#v}
twine upload dist/${PKG}-$VERSION.tar.gz
MD5=$(md5sum dist/${PKG}-${1#v}.tar.gz | cut -f1 -d' ')
cp -r conda-recipe/ dist/conda-recipe-${1#v}
sed -i -E -e "s/version:(.+)/version: $VERSION/" -e "s/path:(.+)/fn: $PKG-$VERSION.tar.gz\n    url: https:\/\/pypi.python.org\/packages\/source\/${PKG:0:1}\/$PKG\/$PKG-$VERSION.tar.gz#md5=$MD5\n    md5: $MD5/" dist/conda-recipe-${1#v}/meta.yaml
env ${PKG_UPPER}_RELEASE_VERSION=$1 python setup.py upload_sphinx

# Specific for this project:
scp -r dist/conda-recipe-${1#v}/ $PKG@hera:~
ssh $PKG@hera "conda-build conda-recipe-${1#v}/"

echo "Remember to build conda-recipe, bump version (and commit and push)!"
