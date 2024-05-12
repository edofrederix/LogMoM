#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

SIGMA=0.5
DSM=2e-3
MESH=1
POLYCELERITY=true

if [ -f "mesh$MESH.tar.gz" ]; then
    tar xzf mesh$MESH.tar.gz -C constant
else
    echo "Mesh file mesh$MESH.tar.gz not found"
    exit
fi

if [ "$POLYCELERITY" == "true" ]; then

    FLUX0="phi0"
    FLUX2="phi2"
    PHASESYSTEM="basicPolyPhaseSystem"
    PHASEMODEL="pureIsothermalPolyPhaseModel"

else

    FLUX0="phi"
    FLUX2="phi"
    PHASESYSTEM="basicMultiphaseSystem"
    PHASEMODEL="pureIsothermalPhaseModel"

fi

cp -r 0.orig 0

VARS="\
    -DVARSIGMA=$SIGMA \
    -DVARDSM=$DSM \
    -DVARFLUX0=$FLUX0 \
    -DVARFLUX2=$FLUX2 \
    -DVARPHASESYSTEM=$PHASESYSTEM \
    -DVARPHASEMODEL=$PHASEMODEL \
    "

find -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    m4 $VARS $IN > $OUT

done

rm -f 0/*.m4

# runApplication mapFields -mapMethod cellPointInterpolate -sourceTime 'latestTime' ../somePrecursorCase
# runApplication setFields

runApplication setLogNormal air $SIGMA $DSM

runApplication decomposePar
