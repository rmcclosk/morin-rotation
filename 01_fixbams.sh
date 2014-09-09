#!/bin/bash

# Fix exome BAM files with malformed header.

. settings.conf

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/01_fixbams
mkdir -p $OUT_DIR

for EXOME_BAM in $(ls $IN_DIR/*exome*.bam); do
    OUT_FILE=$OUT_DIR/$(basename $EXOME_BAM)
    if [[ ! -f $OUT_FILE ]]; then
        $SAMTOOLS_BIN view -H $EXOME_BAM > header.sam
        grep 'VN:\W' header.sam > /dev/null
        if [[ $? -eq 0 ]]; then
            sed 's/VN:\t/VN:0.5.7\t/' header.sam > $OUT_FILE.header
            echo "cp $(readlink -f $EXOME_BAM) $OUT_FILE && \
                  $SAMTOOLS_BIN reheader $OUT_FILE.header $OUT_FILE > $OUT_FILE.fixed" 
        else
            ln -s $(readlink -f $EXOME_BAM) $OUT_FILE
            ln -s $(readlink -f $EXOME_BAM.bai) $OUT_FILE.bai
        fi
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name fixbams
fi
