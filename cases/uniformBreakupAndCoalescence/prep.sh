#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

MODEL=${1:-logmom}
BREAKUPMODEL=${2-LehrMilliesMewes}
COALESCENCEMODEL=${34-LehrMilliesMewes}
DSM=${4:-0.002}
SIGMA=${5:-0.5}
ALPHA=${6:-0.1}
NSECTIONS=${7:-32}

# ------------------------------------------------------------------------------

MATH="import math as m"

ALPHAWATER=$(echo "print(1.0-$ALPHA)" | python3)

if [[ ! "$MODEL" =~ ^(logmom|fpt)$ ]]; then

    echo "Invalid mode (should be logmom or fpt)"
    exit
fi

if [[ ! "$BREAKUPMODEL" =~ ^[a-zA-Z]+$ ]]; then

    echo "Invalid break-up model"
    exit

fi

if [[ ! "$COALESCENCEMODEL" =~ ^[a-zA-Z]+$ ]]; then

    echo "Invalid coalescence model model"
    exit

fi

if [[ ! "$DSM" =~ ^[0-9.e-]+$ ]]; then

    echo "Invalid Sauter mean diameter (should be float)"
    exit

fi

if [[ ! "$SIGMA" =~ ^[0-9.e-]+$ ]]; then

    echo "Invalid distribution width (should be float)"
    exit

fi

if [[ ! "$ALPHA" =~ ^[0-9.e-]+$ ]]; then

    echo "Invalid void fraction (should be float)"
    exit

fi

if [[ ! "$NSECTIONS" =~ ^[0-9]+$ ]]; then

    echo "Invalid number of sections (should be int)"
    exit

fi

cp -r 0.org 0

m4 -DVARALPHA=$ALPHA 0/alpha.air.m4 > 0/alpha.air
m4 -DVARALPHA=$ALPHAWATER 0/alpha.water.m4 > 0/alpha.water

if [ "$MODEL" == "logmom" ]; then

    cp system/probeFields.LogMoM system/probeFields

    m4 \
        -DVARBREAKUPMODEL=$BREAKUPMODEL \
        -DVARCOALESCENCEMODEL=$COALESCENCEMODEL \
        constant/phaseProperties.LogMoM.m4 > constant/phaseProperties

    echo 1 > properties.txt

else

    cp system/probeFields.FPT system/probeFields

    m4 \
        -DVARBREAKUPMODEL=$BREAKUPMODEL \
        -DVARCOALESCENCEMODEL=$COALESCENCEMODEL \
        constant/phaseProperties.FPT.m4 > constant/phaseProperties

    python3 sizeGroups.py $DSM $SIGMA $NSECTIONS

    rm -f constant/FPT/sizeGroups
    mkdir -p constant/FPT

    for ((I=0; I<$NSECTIONS; I++)); do

        J=$((I+1))

        D=($(cat d.txt | head -n $J | tail -n 1))
        F=($(cat f.txt | head -n $J | tail -n 1))

        echo "{dSph $D;}" >> constant/FPT/sizeGroups

        m4 -DVARFIELD=$F -DVARFNAME=f$I 0/f.m4 > 0/f$I.air

        echo f$I.air >> system/probeFields

    done

    echo 2 > properties.txt

fi

echo $DSM >> properties.txt
echo $SIGMA >> properties.txt
echo $ALPHA >> properties.txt

if [ "$MODEL" == "fpt" ]; then

    echo $NSECTIONS >> properties.txt

fi

rm -f 0/*.m4

wmake -s code/populationBalance

runApplication blockMesh

if [ "$MODEL" == "logmom" ]; then

    runApplication setLogNormal air $SIGMA $DSM

fi
