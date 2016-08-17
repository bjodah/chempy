#!/bin/bash -xeu
# Usage:
#
#    $ ./scripts/post_release.sh v1.2.3 myserver githubuser
#
VERSION=${1#v}
SERVER=$2
GITHUBUSER=$3
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
PKG_UPPER=$(echo $PKG | tr '[:lower:]' '[:upper:]')
SDIST_FILE=dist/${PKG}-$VERSION.tar.gz
if [[ ! -f "$SDIST_FILE" ]]; then
    >&2 echo "Nonexistent file $SDIST_FILE"
    exit 1
fi
MD5=$(md5sum "$SDIST_FILE" | cut -f1 -d' ')
cp -r conda-recipe/ dist/conda-recipe-$VERSION
sed -i -E -e "s/version:(.+)/version: $VERSION/" -e "s/path:(.+)/fn: $PKG-$VERSION.tar.gz\n  url: https:\/\/github.com\/$GITHUBUSER\/$PKG\/releases\/download\/v$VERSION\/$PKG-$VERSION.tar.gz\n  md5: $MD5/" dist/conda-recipe-$VERSION/meta.yaml
                                                                                                  
env ${PKG_UPPER}_RELEASE_VERSION=v$VERSION python setup.py upload_sphinx

# Specific for this project:
ssh $PKG@$SERVER 'mkdir -p ~/public_html/conda-recipes/'
scp -r dist/conda-recipe-$VERSION/ "$PKG@$SERVER:~/public_html/conda-recipes/"
ssh $PKG@$SERVER 'mkdir -p ~/public_html/releases/'
scp "$SDIST_FILE" "$PKG@$SERVER:~/public_html/releases/"
for CONDA_PY in 2.7 3.4 3.5; do
    ssh $PKG@$SERVER "source /etc/profile; conda-build --python $CONDA_PY ~/public_html/conda-recipes/conda-recipe-$VERSION/"
done
