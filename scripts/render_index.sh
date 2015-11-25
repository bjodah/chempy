#!/bin/bash -e
#
# Usage (assuming: shopt -s extglob):
#
#  $ cd examples/ && ../scripts/render_index.sh !(index).html
#
mkdir -p thumbs
tmpdir=$(mktemp -d)
trap "rm -r $tmpdir" INT TERM EXIT
cat <<EOF>index.html
<!DOCTYPE html>
<html>
<head>
    <title>Notebook gallery</title>
</head>
<body>
EOF
for f in $@; do
    img=$(basename $f .html).png
    phantomjs $(unset CDPATH && cd "$(dirname "$0")" && echo $PWD)/rasterize.js $f $tmpdir/$img 1200px*900px
    convert $tmpdir/$img -resize 400x300 thumbs/$img
    cat <<EOF>>index.html
<p style='text-align: center'>
<a href='$f' style='ffont-ont-family: sans-serif'>
<img src='thumbs/$img'><br>
$f
</a></p>
EOF
done
cat <<EOF>>index.html
</body>
</html>
EOF
