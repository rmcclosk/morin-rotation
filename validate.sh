#!/bin/bash

. settings.conf

INDIR=$WORK_DIR/00_bams

for BAM_FILE in $(ls $INDIR/*genome*.bam); do
    OUTFILE=$WORKDIR/04_readcounts/$(basename $BAM_FILE | sed s/bam/wig/)
    if [[ ! -f $OUTFILE ]]; then
        echo "$JAVA_BIN -Xmx1048m -jar $PICARD_DIR/ValidateSamFile.jar \
              I=$BAM_FILE IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND"
    fi
done > jobs.txt

exit 0
if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name coverage --qsub "-l h_vmem=2G -l mem_token=2G -l mem_free=2G"
fi
