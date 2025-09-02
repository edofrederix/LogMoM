#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

MODE=${1:-logmom}
MESH=${2:-2}
DT=${3:-10}
DSM=${4:-0.002}
HEATMODEL=${5-NuHZDR}
TURBMODEL_G=${6-kOmegaSST}
TURBMODEL_L=${7-kOmegaSSTSato}

SIGMA=0.5
NSECTIONS=16

# Check input

if [[ ! "$MESH" =~ ^[1-9][0-9]+?S?N?$ ]]; then

    echo "Invalid mesh"
    exit

fi

if [ -f "mesh$MESH.tar.gz" ]; then
    tar xzf mesh$MESH.tar.gz -C constant
elif [ -f "../BWRX300Mesh/mesh$MESH.tar.gz" ]; then
    tar xzf ../BWRX300Mesh/mesh$MESH.tar.gz -C constant
else
    echo "Mesh file mesh$MESH.tar.gz not found."
    echo "Generate package with BWRX300Mesh case."
    exit
fi

if [[ ! "$DT" =~ ^[0-9]+\.?[0-9]*$ ]]; then

    echo "Invalid sub-cooling temperature difference"
    exit

fi

TLOW=$(echo "print(560.15 - $DT)" | python)
TMIN=$(echo "print($TLOW-1.0)" | python)
TMAX=561.15

if [[ ! "$DSM" =~ ^[0-9.e-]+$ ]]; then

    echo "Invalid Sauter mean diameter (should be float)"
    exit

fi

if [[ ! "$MODE" =~ ^(logmom|fpt)$ ]]; then

    echo "Invalid mode (should be logmom or fpt)"
    exit

fi

case $HEATMODEL in

    RanzMarshall|Tomiyama|Hughmark|ChenMayinger|NuHZDR)
        ;;
    *)
        echo "Invalid interfacial heat transfer model (should be RanzMarshall, Tomiyama, Hughmark, ChenMayinger, or NuHZDR)"
        exit
        ;;
esac

case $TURBMODEL_G in

    kEpsilon|kOmegaSST)
        ;;
    *)
        echo "Invalid turbulence model for steam"
        exit
        ;;
esac

case $TURBMODEL_L in

    kEpsilon|kOmegaSST|kOmegaSSTSato)
        ;;
    *)
        echo "Invalid turbulence model for steam"
        exit
        ;;
esac

##

cp -r 0.org 0

VARS=" \
    -DVARTLOW=$TLOW \
    -DVARTMIN=$TMIN \
    -DVARTMAX=$TMAX \
    -DVARSIGMA=$SIGMA \
    -DVARDSM=$DSM \
    -DVARHEATMODEL=$HEATMODEL \
    -DVARTURBMODEL_L=$TURBMODEL_L \
    -DVARTURBMODEL_G=$TURBMODEL_G"

if [ "$MODE" == "logmom" ]; then

    cp system/sampleFields.LogMoM system/sampleFields
    cp system/functions.LogMoM system/functions

    m4 $VARS 0/kappai.steam.m4 > 0/kappai.steam
    m4 $VARS 0/lambda.steam.m4 > 0/lambda.steam

    cp constant/phaseProperties.LogMoM constant/phaseProperties

    cp constant/momentumTransfer.LogMoM constant/momentumTransfer

    m4 $VARS constant/heatTransfer.LogMoM.m4 > constant/heatTransfer

else

    cp system/sampleFields.FPT system/sampleFields
    cp system/functions.FPT system/functions

    python3 sizeGroups.py $DSM $SIGMA $NSECTIONS

    rm -f constant/FPT/sizeGroups
    mkdir -p constant/FPT

    for ((I=0; I<$NSECTIONS; I++)); do

        J=$((I+1))

        D=($(cat d.txt | head -n $J | tail -n 1))
        F=($(cat f.txt | head -n $J | tail -n 1))

        echo "{dSph $D;}" >> constant/FPT/sizeGroups

        m4 -DVARFIELD=$F -DVARFNAME=f$I 0/f.m4 > 0/f$I.steam

        echo f$I.steam >> system/sampleFields

    done

    cp constant/phaseProperties.FPT constant/phaseProperties

    cp constant/momentumTransfer.FPT constant/momentumTransfer

    m4 $VARS constant/heatTransfer.FPT.m4 > constant/heatTransfer

fi

m4 $VARS 0/T.steam.m4 > 0/T.steam
m4 $VARS 0/T.water.m4 > 0/T.water

m4 $VARS constant/momentumTransport.steam.m4 > constant/momentumTransport.steam
m4 $VARS constant/momentumTransport.water.m4 > constant/momentumTransport.water

m4 $VARS system/fvConstraints.m4 > system/fvConstraints

rm -f 0/*.m4

if [ "$MODE" == "logmom" ]; then

    runApplication setLogNormal steam $SIGMA $DSM

fi

runApplication decomposePar
