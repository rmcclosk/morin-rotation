#!/bin/bash

# Fix BAM files with malformed header.

. settings.conf

IN_DIR=$WORK_DIR/00_bams
OUT_DIR=$WORK_DIR/01_fixbams
mkdir -p $OUT_DIR

for BAM_FILE in $(ls $IN_DIR/*.bam); do

    OUT_FILE=$OUT_DIR/$(basename $BAM_FILE)
    if [[ ! -f $OUT_FILE ]]; then
        samtools view -H $BAM_FILE > header.sam
        BODY_RG=$(samtools view $BAM_FILE | head -n 100 | grep -o "RG:Z:[^  ]*" | cut -f 1 | cut -d ':' -f 3 | sort | uniq | tr -d $'\n')
        HEAD_RG=$(grep RG header.sam | cut -f 2 | cut -d ':' -f 2 | sort | tr -d $'\n')

        # Check for missing version number.
        grep "VN:   " header.sam > /dev/null
        if [[ $? -eq 0 ]]; then
            echo "Replacing empty version flag with 0.5.7" >&2
            sed -i 's/VN:/VN:0.5.7/' header.sam
            samtools reheader header.sam $BAM_FILE > $OUT_FILE

        # Check for inconsistency in read groups between body and header.
        elif [[ "x$HEAD_RG" != "x$BODY_RG" || "x$HEAD_RG" == "x" ]]; then
            echo "Going to fix inconsistent read groups" >&2
            echo -n "source $HOME/.bash_profile; "
            echo -n "java -Xmx2048m -jar $PICARD_DIR/AddOrReplaceReadGroups.jar "
            echo -n "VALIDATION_STRINGENCY=SILENT "
            echo -n "INPUT=$BAM_FILE OUTPUT=$OUT_FILE RGLB=N/A "
            echo    "RGPL=Ilumina RGPU=N/A RGSM=N/A"
        else
            ln -s $BAM_FILE $OUT_FILE
            ln -s $BAM_FILE.bai $OUT_FILE.bai
        fi
    else
        echo "Output file $OUT_FILE already exists" >&2
    fi
done > jobs.txt

exit 0
if [[ -s jobs.txt ]]; then
    mqsub --file jobs.txt --chdir qsub-logs --name fixbams --qsub "-l mem_free=3G -l mem_token=3G -l h_vmem=3G"
fi
