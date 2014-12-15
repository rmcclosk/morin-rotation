#!/bin/bash

. settings.conf

BAM_DIR=$WORK_DIR/01_fixbams
OUT_DIR=$WORK_DIR/TITAN

IFS=$'\n'
for LINE in $(tail -n +2 $METADATA | grep -v cell-line); do
    NORMAL=$(echo $LINE | cut -f 2)
    TUMOR=$(echo $LINE | cut -f 3)
    if [[ ! -d $OUT_DIR/$TUMOR ]]; then
        if [[ ! -f $BAM_DIR/$TUMOR.bam ]] || [[ ! -f $BAM_DIR/$NORMAL.bam ]]; then
            continue
        fi
        echo -en "$TUMOR\t$TUMOR\t$BAM_DIR/$TUMOR.bam\t" > titan_input.txt
        echo -e  "$NORMAL\t$NORMAL\t$BAM_DIR/$NORMAL.bam" >> titan_input.txt
        $TITAN_DIR/TitanRunner.py -i $HOME/morin-rotation/titan_input.txt \
                                         --project-name $TUMOR \
                                         --project-path $OUT_DIR \
                                         -c $HOME/morin-rotation/titan_config.cfg
        if [[ $? -eq 0 ]]; then
            mv $OUT_DIR/${TUMOR}_R* $OUT_DIR/$TUMOR
        fi
    fi
done
