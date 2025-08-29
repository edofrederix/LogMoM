#!/bin/bash

MESH=$1

if [[ ! "$MESH" =~ ^[1-9][0-9]+?$ ]]; then

    echo "Invalid mesh size (should be an integer)"
    exit

fi

source $FOAM_SRC/../bin/tools/RunFunctions

runApplication python3 BWRX300.py $MESH -nopopup

runApplication gmshToFoam BWRX300.msh2

runApplication wmake -s deform
runApplication ./deform/deform -overwrite

# Set boundary types

WALLS=(baffles plenum pipes sides)

for WALL in "${WALLS[@]}"; do

    PATCH="wall_${WALL}"

    runApplication -append foamDictionary \
        -entry entry0/$PATCH/type -set wall \
        constant/polyMesh/boundary

done

# Extract symmetry boundaries

runApplication topoSet
runApplication createPatch -overwrite

# Set defaultFaces as the chimney wall

sed -i 's/defaultFaces/wall_chimney/g' constant/polyMesh/boundary

runApplication -append foamDictionary \
    -entry entry0/wall_chimney/type -set wall \
    constant/polyMesh/boundary

cp -r constant/polyMesh .
tar czf mesh${MESH}.tar.gz polyMesh
rm -r polyMesh

cp -r 0.orig 0
