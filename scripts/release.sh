#!/bin/bash -xeu
# Usage:
#
#    $ ./scripts/release.sh v1.2.3 GITHUB_USER GITHUB_REPO
#
# or, using some host specific settings:
#
#    $ export PYODESYS_CVODE_FLAGS="-isystem /opt/sundials-5.3.0-rel-klu-lapack/include -isystem /usr/include/suitesparse" PYODESYS_CVODE_LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,/opt/sundials-5.3.0-rel-klu-lapack/lib:/opt/openblas-0.3.9/lib -L/opt/sundials-5.3.0-rel-klu-lapack/lib -L/opt/openblas-0.3.9/lib -lopenblas"
#    $ ./scripts/release.sh ...

if [[ $1 != v* ]]; then
    echo "Argument does not start with 'v'"
    exit 1
fi
VERSION=${1#v}
find . -type f -iname "*.pyc" -exec rm {} +
find . -type f -iname "*.o" -exec rm {} +
find . -type f -iname "*.so" -exec rm {} +
find . -type d -name "__pycache__" -exec rmdir {} +
./scripts/check_clean_repo_on_master.sh
cd $(dirname $0)/..
# PKG will be name of the directory one level up containing "__init__.py" 
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
! grep --include "*.py" "will_be_missing_in='$VERSION'" -R $PKG/  # see deprecation()
PKG_UPPER=$(echo $PKG | tr '[:lower:]' '[:upper:]')
./scripts/run_tests.sh
env ${PKG_UPPER}_RELEASE_VERSION=v$VERSION python3 setup.py sdist
env ${PKG_UPPER}_RELEASE_VERSION=v$VERSION ./scripts/generate_docs.sh

# All went well, add a tag and push it.
git tag -a v$VERSION -m v$VERSION
git push
git push --tags
twine upload dist/${PKG}-$VERSION.tar.gz

set +x
echo ""
echo "    You may now create a new github release at with the tag \"v$VERSION\", here is a link:"
echo "        https://github.com/$2/${3:-$PKG}/releases/new "
echo "    name the release \"${PKG}-${VERSION}\", and don't foreget to manually attach the file:"
echo "        $(openssl sha256 $(pwd)/dist/${PKG}-${VERSION}.tar.gz)"
echo "    Then run:"
echo ""
echo "        $ ./scripts/post_release.sh $1 $2 myserver.university.edu"
echo ""
