#!/bin/bash

set -e

USER=rmccloskey
HOST=jango.bcgsc.ca
DB=colorectal
PASS=rmccloskey

REF=/extscratch/morinlab/shared/rmccloskey/GRCh37-lite.fa
BAM_DIR=/extscratch/morinlab/shared/rmccloskey/colorectal_exome/01_fixbams
for MAF in mafs/*; do
    LIBNAME=$(basename $MAF | grep -o 'A[0-9]*' | tail -n 1)
    QUERY="SELECT id FROM library WHERE library_name = '$LIBNAME'"
    LIBID=$(mysql $DB -u$USER -p$PASS -s --column-names=FALSE -e "$QUERY")
    NORM=$(basename $MAF | cut -d '_' -f 1 | cut -d '-' -f 1-2)
    TUM=$(basename $MAF | cut -d '_' -f 3 | cut -d '-' -f 1-4)
    if [[ "$NORM" == HT29* ]]; then
        NORM=HT29
        TUM=R3
    elif [[ "$NORM" == "46-WB" ]]; then 
        NORM="046-WB"
    fi
    echo $TUM $NORM
    ./MAF_to_database.py $MAF $LIBID $HOST $DB $USER $PASS --normal-bam $BAM_DIR/$NORM.bam --tumour-bam $BAM_DIR/$TUM.bam --reference $REF
done
