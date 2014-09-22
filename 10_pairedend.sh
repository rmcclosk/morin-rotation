#!/bin/bash

. settings.conf

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/10_pairedend
mkdir -p $OUT_DIR

for BAM_FILE in $(ls $IN_DIR/*genome*.bam); do
    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE)
    if [[ ! -f $OUT_FILE ]]; then
        echo -n "source $HOME/.bash_profile; "
        echo -n "samtools view -h $BAM_FILE | "
        echo -n "$LUMPY_DIR/extractSplitReads_BwaMem -i stdin | "
        echo -n "samtools view -Sb - "
        echo "> $OUT_FILE"
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name pairedend
fi
