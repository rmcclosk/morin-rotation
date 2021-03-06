#!/bin/bash

. ../settings.conf

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/06_unmapped
mkdir -p $OUT_DIR

for BAM_FILE in $(ls $IN_DIR/*genome*.bam); do
    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE | sed s/.bam$/.fastq/)
    if [[ ! -f $OUT_FILE ]]; then
        echo "$SAMTOOLS_BIN view $BAM_FILE | $LUMPY_DIR/split_unmapped_to_fasta.pl -b 20 > $OUT_FILE"
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name unmapped --qsub "-l h_vmem=1G -l mem_token=1G -l mem_free=1G"
fi
