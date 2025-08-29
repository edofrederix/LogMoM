#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase
rm -rf \
    0 \
    constant/polyMesh \
    system/blockMeshDict \
    __pycache__ \
    *.pdf [df].txt \
    system/sampleFields \
    system/averagingFields \
    constant/phaseProperties \
    constant/momentumTransfer \
    constant/FPT \
    system/sixDoFRigidBodyState \
    constant/dynamicMeshDict

find $CASEPATH -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    if [ -f "$OUT" ]; then
        rm $OUT
    fi

done
