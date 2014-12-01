#!/bin/bash

USER=rmccloskey
HOST=jango.bcgsc.ca
DB=colorectal
PASS=rmccloskey

for SEGS in ../TITAN/titan/*_segs.txt.bz2; do
    # get sample id
    SAMPLE=$(basename $SEGS | cut -d '_' -f 1)
    QUERY="SELECT id FROM sample WHERE sample_id = '$SAMPLE'"
    SAMPLE_ID=$(mysql $DB -u$USER -p$PASS -s --column-names=FALSE -e "$QUERY")

    # get library id
    QUERY="SELECT id FROM library WHERE sample_id = '$SAMPLE_ID' and library_type = 'exome'"
    LIBRARY_ID=$(mysql $DB -u$USER -p$PASS -s --column-names=FALSE -e "$QUERY")
    bunzip2 $SEGS
    ./titan_seg_to_database.py $(echo $SEGS | sed s/'.bz2'//) $LIBRARY_ID $DB $HOST $USER $PASS
    bzip2 $(echo $SEGS | sed s/'.bz2'//)
done
