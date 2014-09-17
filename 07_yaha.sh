#!/bin/bash

. settings.conf

IN_DIR=$WORK_DIR/06_unmapped
OUT_DIR=$WORK_DIR/07_yaha
mkdir -p $OUT_DIR

MEM=20G

for FASTQ_FILE in $(ls $IN_DIR/*.fastq); do
    OUT_FILE=$OUT_DIR/$(basename $FASTQ_FILE | sed s/.fastq$/.bam/)
    if [[ ! -f $OUT_FILE ]]; then
        echo -n "source $HOME/.bash_profile; yaha -x $YAHA_REF_INDEX -q $FASTQ_FILE -osh stdout "
        echo -n "-M 15 -H 2000 -L 11 | samtools view -Sb - > "
        echo $OUT_FILE
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name yaha --qsub "-l h_vmem=$MEM -l mem_token=$MEM -l mem_free=$MEM"
fi
