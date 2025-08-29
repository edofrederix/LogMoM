#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase
rm -rf \
    0 \
    constant/polyMesh \
    __pycache__ \
    *.pdf *.txt \
    system/probeFields \
    constant/phaseProperties \
    constant/FPT

find $CASEPATH -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    if [ -f "$OUT" ]; then
        rm $OUT
    fi

done

rm -f code/populationBalance/libpopulationBalance.so
