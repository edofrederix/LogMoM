#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

MODE=${1:-logmom}
MESH=${2:-4}
CASE=${3:-1}
FSI=${4:-false}

TURBMODELAIR=${5:-kOmegaSST}
TURBMODELWATER=${6:-kOmegaSSTSato}

SIGMA=0.5
DSM=2e-3
NSECTIONS=16

if [[ ! "$MESH" =~ ^[1-9][0-9]+?S?N?$ ]]; then

    echo "Invalid mesh"
    exit

fi

if [ -f "mesh$MESH.tar.gz" ]; then
    tar xzf mesh$MESH.tar.gz -C constant
elif [ -f "../PILULEMesh/mesh$MESH.tar.gz" ]; then
    tar xzf ../PILULEMesh/mesh$MESH.tar.gz -C constant
else
    echo "Mesh file mesh$MESH.tar.gz not found."
    echo "Generate package with PILULEMesh case."
    exit
fi

if [[ ! "$CASE" =~ ^([0-9]|1[0-4])$ ]]; then

    echo "Invalid case (should be 0-14)"
    exit

fi

if [[ ! "$MODE" =~ ^(logmom|fpt)$ ]]; then

    echo "Invalid mode (should be logmom or fpt)"
    exit
fi

if [[ ! "$FSI" =~ ^(true|false)$ ]]; then

    echo "Invalid FSI mode (should be true or false)"
    exit
fi

CASEDATA=$(head -n $(($CASE + 1)) caseData.txt | tail -n 1)

JL=$(echo $CASEDATA | cut -d " " -f 4)
JG=$(echo $CASEDATA | cut -d " " -f 5)
ALPHAG=$(echo $CASEDATA | cut -d " " -f 6)

ALPHAL=$(echo "print(1.0-$ALPHAG)" | python)

if [ ! "$CASE" == "0" ]; then
    UG=$(echo "print($JG/$ALPHAG)" | python)
    UL=$(echo "print($JL/$ALPHAL)" | python)
else
    UG=$JL
    UL=$JL
    ALPHAG=1e-5
    ALPHAL=0.99999
fi

cp -r 0.orig 0

if [ "$MODE" == "logmom" ]; then

    cp system/sampleFields.LogMoM system/sampleFields
    cp system/averagingFields.LogMoM system/averagingFields
    cp constant/phaseProperties.LogMoM constant/phaseProperties
    cp constant/momentumTransfer.LogMoM constant/momentumTransfer

else

    cp system/sampleFields.FPT system/sampleFields
    cp system/averagingFields.FPT system/averagingFields
    cp constant/phaseProperties.FPT constant/phaseProperties
    cp constant/momentumTransfer.FPT constant/momentumTransfer

    python3 sizeGroups.py $DSM $SIGMA $NSECTIONS

    rm -f constant/FPT/sizeGroups
    mkdir -p constant/FPT

    for ((I=0; I<$NSECTIONS; I++)); do

        J=$((I+1))

        D=($(cat d.txt | head -n $J | tail -n 1))
        F=($(cat f.txt | head -n $J | tail -n 1))

        echo "{dSph $D;}" >> constant/FPT/sizeGroups

        m4 -DVARFIELD=$F -DVARFNAME=f$I 0/f.m4 > 0/f$I.air

        echo f$I.air >> system/sampleFields

    done

fi

if [ "$FSI" == "true" ]; then

    UAIRCYLINDERBC="movingWallSlipVelocity"
    UWATERCYLINDERBC="movingWallVelocity"

else

    UAIRCYLINDERBC="slip"
    UWATERCYLINDERBC="noSlip"

fi

rm 0/f.m4

VARS="\
    -DVARSIGMA=$SIGMA \
    -DVARDSM=$DSM \
    -DVARALPHAG=$ALPHAG \
    -DVARALPHAL=$ALPHAL \
    -DVARUG=$UG \
    -DVARUL=$UL \
    -DVARTURBMODELAIR=$TURBMODELAIR \
    -DVARTURBMODELWATER=$TURBMODELWATER \
    -DVARUAIRCYLINDERBC=$UAIRCYLINDERBC \
    -DVARUWATERCYLINDERBC=$UWATERCYLINDERBC \
    "

find -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    m4 $VARS $IN > $OUT

done

if [ "$TURBMODELAIR" == "laminar" ]; then

    sed -i 's/simulationType  RAS/simulationType  laminar/g' \
        constant/momentumTransport.air

fi

if [ "$TURBMODELWATER" == "laminar" ]; then

    sed -i 's/simulationType  RAS/simulationType  laminar/g' \
        constant/momentumTransport.water

fi

rm -f 0/*.m4

if [ "$FSI" == "false" ]; then

    rm constant/dynamicMeshDict
    echo > system/sixDoFRigidBodyState
    rm 0/pointDisplacement

fi


if [ "$MODE" == "logmom" ]; then

    runApplication setLogNormal air $SIGMA $DSM

fi

runApplication decomposePar
