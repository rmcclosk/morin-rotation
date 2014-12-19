#!/bin/bash

# Fix BAM file with malformed header. If there are no problems, just link the
# starting file to the end file. Otherwise, fix the problems and create a new
# BAM file. The input BAM file must have an index.
# Usage: ./fixbam.sh in.bam out.bam

BAM_FILE=$1
OUT_FILE=$2

HEADER=$(tempfile)
samtools view -H $BAM_FILE > HEADER
BODY_RG=$(samtools view $BAM_FILE | head -n 100 | grep -o "RG:Z:[^  ]*" | cut -f 1 | cut -d ':' -f 3 | sort | uniq | tr -d $'\n')
HEAD_RG=$(grep RG $HEADER | cut -f 2 | cut -d ':' -f 2 | sort | tr -d $'\n')

# Check for missing version number.
grep "VN:   " header.sam > /dev/null
HAS_VN=$?

# Check for wrong chromosome notation.
grep "chr" header.sam > /dev/null
NO_CHRS=$?

# Header problems.
if [[ $HAS_VN -eq 0 ]] || [[ $NO_CHRS -eq 0 ]]; then
    if [[ $NO_CHRS -eq 0 ]]; then
        echo "Removing chr notation" >&2
        sed -i 's/chrM	/chrMT	/g' header.sam
        sed -i 's/chr//g' header.sam
    fi

    if [[ $HAS_VN -eq 0 ]]; then
        echo "Replacing empty version flag with 0.5.7" >&2
        sed -i 's/VN:/VN:0.5.7/' header.sam
    fi
    samtools reheader header.sam $BAM_FILE > $OUT_FILE

# Check for inconsistency in read groups between body and header.
elif [[ "x$HEAD_RG" != "x$BODY_RG" || "x$HEAD_RG" == "x" ]]; then
    echo "Going to fix inconsistent read groups" 
    java -Xmx2048m -jar $PICARD_DIR/AddOrReplaceReadGroups.jar \
        VALIDATION_STRINGENCY=SILENT \
        INPUT=$BAM_FILE OUTPUT=$OUT_FILE RGLB=N/A \
        RGPL=Ilumina RGPU=N/A RGSM=N/A

else
    ln -s $BAM_FILE $OUT_FILE
    ln -s $BAM_FILE.bai $OUT_FILE.bai
fi
