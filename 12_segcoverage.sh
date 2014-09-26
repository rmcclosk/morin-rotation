#!/bin/bash

. settings.conf

BAM_DIR=$WORK_DIR/01_fixbams
OUT_DIR=$WORK_DIR/12_segcoverage
SEG_DIR=$WORK_DIR/05_hmmcopy/10000
mkdir -p $OUT_DIR

for SEG_FILE in $(ls $SEG_DIR/*/*_kmeans.seg); do
    OUTFILE_STEM=$OUT_DIR/$(basename $SEG_FILE | sed s/.bam$//)
    SAMPLE_ID=$(basename $SEG_FILE | cut -d '_' -f 1)
    BAM_FILE=$BAM_DIR/$(grep $SAMPLE_ID, $GENOME_METADATA | cut -d ',' -f 1,7 | tr ',' '_')_genome.GRCh37-lite.aln.bam
    if [[ ! -f $OUTFILE_STEM.sample_summary ]]; then
        tail -n +2 $SEG_FILE | tr -d '"' | cut -d ' ' -f 1-3 | sed -e 's/ /:/' -e 's/ /-/' > $OUTFILE_STEM.list
        echo "$JAVA_BIN -Xmx2048m -jar $GATK_JAR -T DepthOfCoverage -omitBaseOutput \
            -omitLocusTable -R $HUMAN_REF -I $BAM_FILE -L $OUTFILE_STEM.list \
            -o $OUTFILE_STEM"
    else
        echo "Output files with stem $OUTFILE_STEM already exist" >&2
    fi
done > jobs.txt

exit 0
if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name coverage --qsub "-l h_vmem=3G -l mem_token=3G -l mem_free=3G"
fi
