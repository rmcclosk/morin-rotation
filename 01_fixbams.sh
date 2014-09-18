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
        grep "VN:   " header.sam > /dev/null
        if [[ $? -eq 0 ]]; then
            echo "Replacing empty version flag with 0.5.7" >&2
            sed -i 's/VN:/VN:0.5.7/' header.sam
            samtools reheader header.sam $EXOME_BAM > $OUT_FILE
        fi
    
        # Check for inconsistency in read groups between body and header.
        BODY_RG=$(samtools view $EXOME_BAM | head -n 100 | grep -o "RG:Z:[^  ]*" | sort | uniq | cut -d ':' -f 3)
        HEAD_RG=$(grep RG header.sam | cut -f 2 | cut -d ':' -f 2)
        if [[ "x$HEAD_RG" -ne "x$BODY_RG" ]]; then
            echo -n "source $HOME/.bash_profile; "
            echo -n "java -Xmx2048m -jar $PICARD_DIR/AddOrReplaceReadGroups.jar "
            echo -n "INPUT=$EXOME_BAM OUTPUT=$OUT_FILE.regrouped.bam RGLB=N/A "
            echo -n "RGPL=Ilumina RGPU=N/A RGSM=N/A"
        fi
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

exit 0
if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name fixbams
fi
