#!/bin/bash

. settings.conf

INDIR=$WORK_DIR/01_fixbams
OUTDIR=$WORK_DIR/02_coverage
mkdir -p $OUTDIR

for BAM_FILE in $(ls $INDIR/*exome*.bam); do
    OUTFILE_STEM=$OUTDIR/$(basename $BAM_FILE | sed s/.bam$//)
    if [[ ! -f $OUTFILE_STEM.sample_summary ]]; then
        echo "$JAVA_BIN -Xmx2048m -jar $GATK_JAR -T DepthOfCoverage -omitBaseOutput \
            -omitLocusTable -R $HUMAN_REF -I $BAM_FILE -L $INTERVAL_LIST \
            -o $OUTFILE_STEM"
    else
        echo "Output files with stem $OUTFILE_STEM already exist" >&2
    fi
done > jobs.txt

exit 0
if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name coverage --qsub "-l h_vmem=3G -l mem_token=3G -l mem_free=3G"
fi
