#!/bin/bash

# Fix exome BAM files with malformed header.

. settings.conf

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/01_fixbams
mkdir -p $OUT_DIR

for EXOME_BAM in $(ls $IN_DIR/*.bam); do
    if [[ $EXOME_BAM == *genome* ]]; then continue; fi

    OUT_FILE=$OUT_DIR/$(basename $EXOME_BAM)
    if [[ ! -f $OUT_FILE ]]; then
        samtools view -H $EXOME_BAM > header.sam

        # Check for missing version number.
        grep 'VN:\W' header.sam > /dev/null
        GOOD_VN=$?

        grep 'RG' header.sam > /dev/null
        MISSING_RG=$?

        if [[ $GOOD_VN -eq 0 ]]; then
            sed -i 's/VN:\t/VN:0.5.7\t/' header.sam 
        fi
        if [[ $MISSING_RG -eq 1 ]]; then
            echo "@RG\tID:1" >> header.sam
        fi

        if [[ $GOOD_VN -eq 0 || $MISSING_RG -eq 1 ]]; then
            cp header.sam $OUT_FILE.header
            echo -n "source $HOME/.bash_profile; "
            echo -n "samtools reheader $OUT_FILE.header "
            echo "$(readlink -f $EXOME_BAM) > $OUT_FILE" 
        else
            ln -s $(readlink -f $EXOME_BAM) $OUT_FILE
            ln -s $(readlink -f $EXOME_BAM.bai) $OUT_FILE.bai
        fi
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt
exit 0

if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name fixbams
fi
