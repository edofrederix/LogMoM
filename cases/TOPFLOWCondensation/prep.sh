#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

# Case setup
# 1: P = 10 bar, Jl = 1.017, Jg = 0.219, Tsub = 3.9, Dori = 1 mm
# 2: P = 10 bar, Jl = 1.017, Jg = 0.219, Tsub = 5.0, Dori = 1 mm
# 3: P = 20 bar, Jl = 1.017, Jg = 0.219, Tsub = 3.7, Dori = 1 mm
# 4: P = 20 bar, Jl = 1.017, Jg = 0.219, Tsub = 6.0, Dori = 1 mm
# 5: P = 20 bar, Jl = 1.017, Jg = 0.219, Tsub = 6.0, Dori = 4 mm
# 6: P = 40 bar, Jl = 1.017, Jg = 0.219, Tsub = 5.0, Dori = 1 mm

MODE=${1:-logmom}
CASE=${2:-1}
HEATMODEL=${3-NuHZDR}
TURBMODEL=${4-kOmegaSSTSato}

# Check input

if [[ ! "$CASE" =~ ^[1-6]$ ]]; then

    echo "Invalid case (should be 1-6)"
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

case $TURBMODEL in

    kEpsilon|kOmegaSST|kOmegaSSTSato)
        ;;
    *)
        echo "Invalid turbulence model"
        exit
        ;;
esac

##

SIGMA=0.5
NSECTIONS=16

# Water velocities
VWATERS=(1.05450 1.03663 1.08603 1.04195 1.08517 1.09046)

# Average bubble diameter at inlet
DBUBS=(0.0117414 0.0081615 0.0099606 0.0063471 0.016852 0.0073394)

# Flow Properties
PRESSURES=(1e6 1e6 2e6 2e6 2e6 4e6)

TVAP=(453.03 453.03 485.53 485.53 485.53 523.50)
TLIQ=(449.13 448.03 481.83 479.53 479.53 518.50)

CPVAP=(2711.38 2711.38 3191.04 3191.04 3191.04 4020.32)
CPLIQ=(4404.48 4404.48 4565.5 4565.5 4565.5 4873.91)

HVAP=(2777110 2777110 2798290 2798290 2798290 2800820)
HLIQ=(762515 762515 908498 908498 908498 1087490)

RHOVAP=(5.145 5.145 10.0417 10.0417 10.0417 20.09)
RHOLIQ=(887.129 887.129 849.798 849.798 849.798 798.368)

MUVAP=(15.0213e-6 15.0213e-6 16.1441e-6 16.1441e-6 16.1441e-6 17.5079e-6)
MULIQ=(150.244e-6 150.244e-6 126.102e-6 126.102e-6 126.102e-6 105.951e-6)

PRVAP=(1.11807 1.11807 1.21009 1.21009 1.21009 1.37056)
PRLIQ=(0.982735 0.982735 0.878538 0.878538 0.878538 0.831835)

TR=(453.03 453.03 485.53 485.53 485.53 523.50)

STEN=(0.0422174 0.0422174 0.0348321 0.0348321 0.0348321 0.0259597)

mkdir -p constant/boundaryData/inlet/0

cp data/Case$CASE/Alpha/alpha.steam constant/boundaryData/inlet/0/alpha.steam
cp data/Case$CASE/Alpha/alpha.water constant/boundaryData/inlet/0/alpha.water
cp data/Case$CASE/Vg/U.steam constant/boundaryData/inlet/0/U.steam

cp -r 0.org 0

C=$(($CASE-1))

m4 -DVARVWATER=${VWATERS[$C]} 0/U.water.m4 > 0/U.water
m4 -DVARTLIQ=${TLIQ[$C]} 0/T.water.m4 > 0/T.water
m4 -DVARTVAP=${TVAP[$C]} 0/T.steam.m4 > 0/T.steam
m4 -DVARPRESSURE=${PRESSURES[$C]} 0/p.m4 > 0/p
m4 -DVARPRESSURE=${PRESSURES[$C]} 0/p_rgh.m4 > 0/p_rgh

m4 -DVARSIGMA=$SIGMA -DVARDSM=${DBUBS[$C]} 0/A.steam.m4 > 0/A.steam
m4 -DVARSIGMA=$SIGMA -DVARDSM=${DBUBS[$C]} 0/N.steam.m4 > 0/N.steam

m4 -DVARTR=${TR[$C]} constant/saturationTemperature.m4 > constant/saturationTemperature

m4  -DVARRHOLIQ=${RHOLIQ[$C]} \
    -DVARCPLIQ=${CPLIQ[$C]} \
    -DVARTR=${TR[$C]} \
    -DVARHLIQ=${HLIQ[$C]} \
    -DVARMULIQ=${MULIQ[$C]} \
    -DVARPRLIQ=${PRLIQ[$C]} \
    constant/thermophysicalProperties.water.m4 > \
    constant/thermophysicalProperties.water

m4  -DVARCPVAP=${CPVAP[$C]} \
    -DVARTR=${TR[$C]} \
    -DVARHVAP=${HVAP[$C]} \
    -DVARRHOVAP=${RHOVAP[$C]} \
    -DVARMUVAP=${MUVAP[$C]} \
    -DVARPRVAP=${PRVAP[$C]} \
    constant/thermophysicalProperties.steam.m4 > \
    constant/thermophysicalProperties.steam

if [ "$MODE" == "logmom" ]; then

    cp system/sampleFields.LogMoM system/sampleFields
    cp system/functions.LogMoM system/functions

    m4  -DVARSTEN=${STEN[$C]} \
        -DVARPRESSURES=${PRESSURES[$C]} \
        -DVARTR=${TR[$C]} \
        -DVARHEATMODEL=$HEATMODEL \
        constant/phaseProperties.LogMoM.m4 > constant/phaseProperties

    cp constant/momentumTransfer.LogMoM constant/momentumTransfer

    m4  -DVARHEATMODEL=$HEATMODEL \
        constant/heatTransfer.LogMoM.m4 > constant/heatTransfer

else

    cp system/sampleFields.FPT system/sampleFields
    cp system/functions.FPT system/functions

    m4  -DVARSTEN=${STEN[$C]} \
        -DVARPRESSURES=${PRESSURES[$C]} \
        -DVARTR=${TR[$C]} \
        -DVARHEATMODEL=$HEATMODEL \
        constant/phaseProperties.FPT.m4 > constant/phaseProperties

    cp constant/momentumTransfer.FPT constant/momentumTransfer

    m4  -DVARHEATMODEL=$HEATMODEL \
        constant/heatTransfer.FPT.m4 > constant/heatTransfer

    python3 sizeGroups.py ${DBUBS[$C]} $SIGMA $NSECTIONS

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

fi

m4 -DVARTURBMODEL=$TURBMODEL constant/momentumTransport.water.m4 \
    > constant/momentumTransport.water

rm -f 0/*.m4

runApplication blockMesh

if [ "$MODE" == "logmom" ]; then

    runApplication setLogNormal steam $SIGMA ${DBUBS[$C]}

fi

runApplication decomposePar
