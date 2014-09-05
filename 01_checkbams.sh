#!/bin/bash

. settings.conf

OUT_FILE=$WORK_DIR/01_checkbams.csv

if [[ -f $OUT_FILE ]]; then
	echo "Output file $OUT_FILE exists."
	exit 1
fi

rm -f jobs.txt
for BAM_FILE in $BAM_DIR/*.bam; do
	COMMAND='export _JAVA_OPTIONS="-Xmx512M"'
	COMMAND="$COMMAND; $JAVA_BIN -jar $PICARD_DIR/ValidateSamFile.jar \
		IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND \
		VALIDATE_INDEX=false INPUT=$BAM_FILE MAX_OUTPUT=1"
	COMMAND="$COMMAND; echo $(basename $BAM_FILE),\$? >> $OUT_FILE"
	echo $COMMAND >> jobs.txt
	break
done

mqsub --file jobs.txt --chdir qsub-logs --qsub "-l h_vmem=1G -l s_vmem=1G -l mem_token=1G -l mem_free=1G" 
