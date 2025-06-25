#!/usr/bin/env bash

PROJECT_DIR=${PWD}

if [ ! -d src/htslib ]; then
    mkdir ${PROJECT_DIR}/src/htslib
fi

(cd ${PROJECT_DIR}/htslib-1.21/ && \
./configure --prefix ${PROJECT_DIR}/src/htslib/ && \
make && \
make install) && \
rm -rf ${PROJECT_DIR}/htslib-1.21/
