#!/bin/bash

. settings.conf

WORK_DIR=$(echo $WORK_DIR | sed s/dlbcl/colorectal/)

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/09_flagstat
mkdir -p $OUT_DIR

for BAM_FILE in $IN_DIR/*.bam; do
    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE | sed s/bam$/txt/)
    echo "source $HOME/.bash_profile; samtools flagstat $BAM_FILE > $OUT_FILE"
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name flagstat
fi
