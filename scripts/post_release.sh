#!/bin/bash -xeu
# Usage:
#
#    $ ./scripts/post_release.sh v1.2.3 myserver githubuser
#
VERSION=${1#v}
SERVER=$2
GITHUBUSER=$3
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
MD5=$(md5sum dist/${PKG}-$VERSION.tar.gz | cut -f1 -d' ')
cp -r conda-recipe/ dist/conda-recipe-$VERSION
sed -i -E -e "s/version:(.+)/version: $VERSION/" -e "s/path:(.+)/fn: $PKG-$VERSION.tar.gz\n  url: https:\/\/github.com\/$GITHUBUSER\/$PKG\/releases\/download\/v$VERSION\/$PKG-$VERSION.tar.gz\n  md5: $MD5/" dist/conda-recipe-$VERSION/meta.yaml
                                                                                                  
env ${PKG_UPPER}_RELEASE_VERSION=v$VERSION python setup.py upload_sphinx

# Specific for this project:
scp -r dist/conda-recipe-$VERSION/ $PKG@$SERVER:~/public_html/conda-recipes/
scp dist/${PKG}-$VERSION.tar.gz $PKG@$SERVER:~/public_html/releases/
for CONDA_PY in 2.7 3.4 3.5; do
    ssh $PKG@$SERVER "source /etc/profile; conda-build --python $CONDA_PY ~/public_html/conda-recipes/conda-recipe-$VERSION/"
done
