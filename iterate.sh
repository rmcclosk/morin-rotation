#!/bin/bash
set -e

ITER=1
while [[ $ITER -le 100 ]]; do
    cd ExomeCNV
    make -j 10
    cp segments/*.seg.bz2 ../ABSOLUTE/segments
    mv segments segments-$ITER
    mkdir segments
    
    cd ../ABSOLUTE
    make -j 10
    IFS=$'\n'
    for LINE in $(cat summary.PP-calls_tab.txt); do
        echo 1-`echo $LINE | cut -f 4` | bc > ../ExomeCNV/admixture/`echo $LINE | cut -f 1`.txt
    done
    mv summary.PP-calls_tab.txt summary.$ITER.PP-calls_tab.txt
    rm absolute/* -rf
    ITER=$(($ITER+1))

    cd ..
done
