#!/bin/bash

# Fix corrupted BAM files.

. settings.conf

for BAM_FILE in $BAM_DIR; do
    $SAMTOOLS_BIN view -H $BAM_FILE > header.sam
    if [[ $(grep 'VN\W' header.sam) -eq 0 ]]; then
        sed -i header.sam 's/VN:\t/VN:0.5.7\t/'
        samtools reheader header.sam $(readlink -f $BAM_FILE) > $BAM_FILE
    fi
done
