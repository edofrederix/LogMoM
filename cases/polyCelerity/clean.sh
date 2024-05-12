#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase
rm -rf \
    0 \
    constant/phaseProperties \
    *.txt \
    constant/FPT \
    constant/momentumTransport.air[0-9] \
    constant/momentumTransport.air[0-9][0-9] \
    constant/thermophysicalProperties.air[0-9] \
    constant/thermophysicalProperties.air[0-9][0-9] \
    system/functions \
    system/sampleFields

find $CASEPATH -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    if [ -f "$OUT" ]; then
        rm $OUT
    fi

done
