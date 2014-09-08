#!/bin/bash

. settings.conf

INDIR=$WORK_DIR/00_bams
OUTDIR=$WORK_DIR/01_coverage
mkdir -p $OUTDIR

for BAM_FILE in $(ls $INDIR/*.bam); do
    OUTFILE_STEM=$OUTDIR/$(basename $BAM_FILE | sed s/.bam$//)
    if [[ ! -f $OUTFILE_STEM.sample_summary ]]; then
        echo "$JAVA_BIN -jar $GATK_JAR -T DepthOfCoverage -omitBaseOutput \
            -omitLocusTable -R $HUMAN_REF -I $BAM_FILE -L $INTERVAL_LIST \
            -o $OUTFILE_STEM"
    else
        echo "Output files with stem $OUTFILE_STEM already exist" >&2
    fi
done
