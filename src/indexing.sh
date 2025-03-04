#!/bin/bash

which bgzip > /dev/null
bgzip=$?

which tabix > /dev/null
tabix=$?

which bcftools > /dev/null
bcftools=$?

indexing(){
    local input="$1"
    local output="$1.gz"

    if [[ $bgzip -eq 0  && $tabix -eq 0 ]]; then
        bgzip -c $input > $output && tabix -p vcf $output
    elif [ $bcftools -eq 0 ]; then
        bcftools view -O z -o $output $input && bcftools index $output
    else
        exit 1
    fi
    if [ $? -ne 0 ]; then
        exit 1
    fi
}

indexing $1

