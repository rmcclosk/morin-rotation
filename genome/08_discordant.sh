#!/bin/bash

. ../settings.conf

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/08_discordant
mkdir -p $OUT_DIR

for BAM_FILE in $(ls $IN_DIR/*genome*.bam); do
    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE)
    if [[ ! -f $OUT_FILE ]]; then
        echo -n "source $HOME/.bash_profile; "
        echo -n "samtools view -u -F 0x0002 $BAM_FILE | "
        echo -n "samtools view -u -F 0x0100 - | "
        echo -n "samtools view -u -F 0x0004 - | "
        echo -n "samtools view -u -F 0x0008 - | "
        echo -n "samtools view -b -F 0x0400 - > "
        echo    "$OUT_FILE"
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name discordant
fi
