#!/bin/bash

. ../settings.conf

IN_DIR=$WORK_DIR/01_fixbams
OUT_DIR=$WORK_DIR/08_discordant
mkdir -p $OUT_DIR

for BAM_FILE in $(ls $IN_DIR/*.bam); do
    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE)
    if [[ ! -f $OUT_FILE ]]; then
        echo -n "source $HOME/.bash_profile; "
        echo -n "samtools view -u -F 0x0002 $BAM_FILE | "
        echo -n "samtools view -u -F 0x0100 - | "
        echo -n "samtools view -u -F 0x0004 - | "
        echo -n "samtools view -u -F 0x0008 - | "
        echo -n "samtools view -b -F 0x0400 - > "
        echo    "$OUT_FILE"
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name discordant
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
