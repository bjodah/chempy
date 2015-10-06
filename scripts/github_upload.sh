#!/bin/bash -eux
#
# Careful - this script rebases and pushes forcefully!
#
# Remember to set user and email in git:
#
#   $ git config --global user.name "First Lastname"
#   $ git config --global user.email "first.lastname@email.domain"
#

UPLOAD_DIR=$1
GITHUB_USER=$2
GITHUB_REPO=$3
OVERWRITE_UPLOAD_BRANCH=$4
WORKDIR=$(pwd)
TMPDIR=$(mktemp -d)
trap "rm -r ${TMPDIR}" EXIT SIGINT SIGTERM
cd $TMPDIR
git clone --quiet git@github.com:${GITHUB_USER}/${GITHUB_REPO} $TMPDIR > /dev/null
git checkout --orphan $OVERWRITE_UPLOAD_BRANCH
git rm -rf . > /dev/null
cd $WORKDIR
cp -r ${UPLOAD_DIR}/. $TMPDIR/
cd $TMPDIR
git add -f . > /dev/null
git commit -m "Lastest docs from successful drone build (hash: ${DRONE_COMMIT})"
git push -f origin $OVERWRITE_UPLOAD_BRANCH
