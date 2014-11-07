#!/bin/bash

. ../settings.conf

IN_DIR=$WORK_DIR/01_fixbams
OUT_DIR=$WORK_DIR/10_pairedend
mkdir -p $OUT_DIR

for BAM_FILE in $(ls $IN_DIR/*.bam); do
    OUT_FILE=$(echo $OUT_DIR/$(basename $BAM_FILE) | sed s/bam$/histo/)
    if [[ ! -f $OUT_FILE ]]; then
        samtools view $BAM_FILE | tail -n +100000 | \
            $LUMPY_DIR/pairend_distro.py -r 100 -X 4 -N 10000 -o $OUT_FILE \
            > $(echo $OUT_FILE | sed s/histo/dat/)
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done 
