#!/bin/bash

. ../settings.conf

IN_DIR=$WORK_DIR/06_unmapped
OUT_DIR=$WORK_DIR/07_yaha
mkdir -p $OUT_DIR

MEM=2G
NCPUS=12

for FASTQ_FILE in $(ls $IN_DIR/*.fastq); do
    OUT_FILE=$OUT_DIR/$(basename $FASTQ_FILE | sed s/.fastq$/.bam/)
    if [[ ! -f $OUT_FILE ]]; then
        echo -n "source $HOME/.bash_profile; yaha -t 4 -x $YAHA_REF_INDEX -q $FASTQ_FILE -osh stdout "
        echo -n "-M 15 -H 2000 -L 11 | samtools view -Sb - > "
        echo $OUT_FILE
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name yaha --qsub "-pe ncpus $NCPUS -l h_vmem=$MEM -l mem_token=$MEM -l mem_free=$MEM"
fi

for BAM_FILE in $(ls $OUT_DIR/*.bam); do
    samtools view -H $BAM_FILE | grep SO:coordinate > /dev/null
    if [[ $? -ne 0 ]]; then
        echo "source $HOME/.bash_profile; samtools sort -O bam -T $BAM_FILE -o $BAM_FILE.sort $BAM_FILE && mv $BAM_FILE.sort $BAM_FILE"
    else
        echo "BAM file $BAM_FILE is already sorted" >&2
    fi
done > jobs.txt

NCPUS=1
MEM=4G
if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name sort --qsub "-pe ncpus $NCPUS -l h_vmem=$MEM -l mem_token=$MEM -l mem_free=$MEM"
fi
