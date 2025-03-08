#!/bin/bash

# Check if bgzip, tabix, and bcftools are installed
which bgzip > /dev/null
bgzip=$?

which tabix > /dev/null
tabix=$?

which bcftools > /dev/null
bcftools=$?

# Indexing VCF file
indexing(){
    local input="$1"
    local output="$1.gz"

    # If bgzip and tabix are installed
    if [[ $bgzip -eq 0  && $tabix -eq 0 ]]; then
        bgzip -c $input > $output && tabix -p vcf $output
    # If bcftools is installed
    elif [ $bcftools -eq 0 ]; then
        bcftools view -O z -o $output $input && bcftools index $output
    # If none of the tools are installed quit
    else
        exit 1
    fi
    # Check if indexing was successful
    if [ $? -ne 0 ]; then
        exit 1
    fi
}

indexing $1

