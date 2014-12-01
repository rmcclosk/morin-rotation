#!/bin/bash

set -e

USER=rmccloskey
HOST=jango.bcgsc.ca
DB=colorectal
PASS=rmccloskey

REF=/extscratch/morinlab/shared/rmccloskey/GRCh37-lite.fa
BAM_DIR=/extscratch/morinlab/shared/rmccloskey/colorectal_exome/01_fixbams
for MAF in ../MAF/by-sample/*; do
    SAMPLE=$(basename $MAF | cut -d '.' -f 1)
    QUERY="SELECT library.id FROM library JOIN sample ON library.sample_id = sample.id WHERE sample.sample_id = '$SAMPLE' AND library.library_type = 'exome'"
    LIBID=$(mysql $DB -u$USER -p$PASS -s --column-names=FALSE -e "$QUERY")
    TUM=$(basename $MAF | cut -d '.' -f 1)
    NORM=$(echo $TUM | cut -d '-' -f 3)-BF
    if [[ "$TUM" == R3* ]]; then
        NORM=HT29
        TUM=R3
    elif [[ "$NORM" == "046-BF" ]]; then 
        NORM="046-WB"
    fi
    ./MAF_to_database.py $MAF $LIBID $HOST $DB $USER $PASS --normal-bam $BAM_DIR/$NORM.bam --tumour-bam $BAM_DIR/$TUM.bam --reference $REF
done
