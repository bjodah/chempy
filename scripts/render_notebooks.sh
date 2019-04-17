#!/bin/bash -e
function quiet_unless_fail {
    # suppresses function output unless exit status is != 0
    OUTPUT_FILE=$(tempfile)
    #/bin/rm --force /tmp/suppress.out 2>/dev/null
    EXECMD=${1+"$@"}
    $EXECMD > ${OUTPUT_FILE} 2>&1
    EXIT_CODE=$?
    if [ ${EXIT_CODE} -ne 0 ]; then
	cat ${OUTPUT_FILE}
	echo "The following command exited with exit status ${EXIT_CODE}: ${EXECMD}"
	/bin/rm ${OUTPUT_FILE}
	exit $?
    fi
    /bin/rm ${OUTPUT_FILE}
}

if [ -f index.ipynb ]; then
    sed -i.bak0 's/ipynb/html/' index.ipynb
    sed -i.bak1 's/filepath=index.html/filepath=index.ipynb/' index.ipynb  # mybinder link fix
fi
set +e
for dir in . examples/; do
    cd $dir
    for fname in *.ipynb; do
        echo "rendering ${fname}..."
        quiet_unless_fail jupyter nbconvert --debug --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 "${fname}" || exit 1
    done
    cd -
done
set -e
cd examples/
../scripts/render_index.sh *.html
