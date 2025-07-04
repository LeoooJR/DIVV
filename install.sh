#!/usr/bin/env bash

INSTALL_DIR=$1
CURRENT_DIR=$(pwd)

if [ ! -d ${INSTALL_DIR}/htslib ]; then
    mkdir -p ${INSTALL_DIR}/htslib
fi

(cd ${CURRENT_DIR}/htslib-1.21/ && \
./configure --prefix ${INSTALL_DIR}/htslib/ && \
make && \
make install) && \
rm -rf ${CURRENT_DIR}/htslib-1.21/
