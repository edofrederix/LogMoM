#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase
rm -rf 0
rm -rf postProcessing/*
rm -f *.pdf
rm -f *.png
rm -f *.txt
rm -f *.obj
rm -f *.msh
rm -rf processor*
rm -rf __pycache__

find $CASEPATH -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    if [ -f "$OUT" ]; then
        rm $OUT
    fi

done
