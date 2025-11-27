#!/bin/bash

MESH=$1

if [[ ! "$MESH" =~ ^[1-9][0-9]+?$ ]]; then

    echo "Invalid mesh size (should be an integer)"
    exit

fi

source $FOAM_SRC/../bin/tools/RunFunctions

runApplication python3 PILULE.py $MESH -nopopup

runApplication gmshToFoam PILULE.msh2

runApplication wmake -s deform
runApplication ./deform/deform

# Set boundary types

WALLS=(wall_pipe wall_cylinder)

for WALL in "${WALLS[@]}"; do

    runApplication -append foamDictionary \
        -entry entry0/$WALL/type -set wall \
        constant/polyMesh/boundary

done

cp -r constant/polyMesh .
tar czf mesh${MESH}.tar.gz polyMesh
rm -r polyMesh

cp -r 0.orig 0
