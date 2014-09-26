#!/bin/bash

. settings.conf

WIN_SIZE=$1
if [[ "x$WIN_SIZE" == "x" ]]; then
    WIN_SIZE=1000
fi

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/04_readcounts/$WIN_SIZE
mkdir -p $OUT_DIR

GC_FILE=$OUT_DIR/$(basename $HUMAN_REF | sed 's/fa[sta]*$/gc.wig/')
if [[ ! -f $GC_FILE ]]; then
    echo "Making GC content file."
    gcCounter -w $WIN_SIZE $HUMAN_REF > $GC_FILE
fi

MAP_FILE=$OUT_DIR/$(basename $MAPPABILITY_FILE | sed s/bw/wig/)
if [[ ! -f $MAP_FILE ]]; then
    echo "Making mappability file."
    mapCounter -w $WIN_SIZE $MAPPABILITY_FILE > $MAP_FILE
fi

for BAM_FILE in $(ls $IN_DIR/*.bam); do
    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE | sed s/.bam$/.wig/)
    if [[ $BAM_FILE == *exome* ]]; then continue; fi

    if [[ ! -f $OUT_FILE ]]; then
        echo "source $HOME/.bash_profile; readCounter -w $WIN_SIZE -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $BAM_FILE > $OUT_FILE" >> jobs.txt
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name readCounter --qsub "-l h_vmem=1G -l mem_token=1G -l mem_free=1G"
fi
