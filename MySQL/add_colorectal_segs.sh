#!/bin/bash

set -e

USER=rmccloskey
HOST=jango.bcgsc.ca
DB=colorectal
PASS=rmccloskey

for SEGS in ../TITAN/titan/*_segs.txt.bz2; do
    SAMPLE=$(basename $SEGS | cut -d '_' -f 1)
    echo $SAMPLE
    #bunzip2 $SEGS
    QUERY="SELECT id FROM library WHERE sample = '$LIBNAME'"
    #LIBID=$(mysql $DB -u$USER -p$PASS -s --column-names=FALSE -e "$QUERY")
    ./titan_segs_to_database.py $SEGS somatic exome 
done
