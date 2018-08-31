#!/bin/bash
jupyter nbconvert --debug --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 $1 2>&1 >$1.log
if [[ $? -ne 0 ]]; then
    2>&1 echo "Failed to render $1"
    cat $1.log
    exit 1
fi
