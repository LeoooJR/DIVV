#!/usr/bin/env bash

PROJECT_DIR=$(PWD)

if [ ! -d htslib ]; then
    mkdir ${PROJECT_DIR}/htslib
fi

(cd ${PROJECT_DIR}/htslib-1.21/ && \
./configure --prefix ${PROJECT_DIR}/htslib/ && \
make && \
make install) && \
rm -rf ${PROJECT_DIR}/htslib-1.21/
