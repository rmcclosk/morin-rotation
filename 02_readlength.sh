#!/bin/bash

. settings.conf

OUT_DIR=$WORK_DIR/02_readlength
mkdir -p $OUT_DIR

rm -f jobs.txt
for BAM_FILE in $(ls $EXOME_BAM_DIR/*.bam); do
    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE | sed s/.bam$/.tbl/)
    if [[ ! -f $OUT_FILE ]]; then
        echo "$JAVA_BIN -Xmx512M -jar $GATK_JAR \
              -T ReadLengthDistribution \
              -I $BAM_FILE \
              -R $HUMAN_REF \
              -o $OUT_DIR/$(basename $BAM_FILE | sed s/.bam$/.tbl/)" >> jobs.txt
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done

mqsub --file jobs.txt --name readlength --chdir qsub-logs --qsub "-l s_vmem=1G -l h_vmem=1G -l mem_token=1G -l mem_free=1G"
