
#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

ALPHAINLET=0.01

MODE=$1
DSM=$2
SIGMA=$3

if [[ ! "$MODE" =~ ^(logmom|fpt)$ ]]; then

    echo "Invalid mode (should be logmom or fpt)"
    exit
fi

if [[ ! "$DSM" =~ ^[0-9]+(\.[0-9]+?)?$ ]]; then

    echo "Invalid value set for dsm (should be float)"
    exit 1

fi

if [[ ! "$SIGMA" =~ ^[0-9]+(\.[0-9]+?)?$ ]]; then

    echo "Invalid value set for sigma (should be float)"
    exit 1

fi

echo "Mode = $MODE, dsm = $DSM, sigma = $SIGMA"

# ------------------------------------------------------------------------------

cp -r 0.org 0

if [ "$MODE" == "logmom" ]; then

    POLYCELERITY=$4

    if [[ ! "$POLYCELERITY" =~ ^(true|false)$ ]]; then

        echo "Specify if simulation is poly-celeric (true or false)"
        exit 1

    fi

    cp system/sampleFields.LogMoM system/sampleFields
    cp system/functions.LogMoM system/functions

    if [[ "$POLYCELERITY" == "true" ]]; then

        FLUX0="phi0"
        FLUX2="phi2"
        PHASESYSTEM="basicPolyPhaseSystem"
        PHASEMODEL="pureIsothermalPolyPhaseModel"

    else

        FLUX0="phi"
        FLUX2="phi"
        PHASESYSTEM="basicMultiphaseSystem"
        PHASEMODEL="pureIsothermalPhaseModel"

        sed -i 's/U0\.air//g' system/sampleFields
        sed -i 's/U2\.air//g' system/sampleFields

        sed -i 's/phi2/phi/g' system/functions
        sed -i 's/phi0/phi/g' system/functions

    fi

    m4 -DVARPHASESYSTEM=$PHASESYSTEM -DVARPHASEMODEL=$PHASEMODEL \
        constant/phaseProperties.LogMoM.m4 > constant/phaseProperties

    m4 -DVARPHASENAME=air -DVARINLET=$ALPHAINLET 0/alpha.air.m4 > 0/alpha.air
    m4 -DVARPHASENAME=air 0/U.air.m4 > 0/U.air
    m4 -DVARPHASENAME=air 0/T.air.m4 > 0/T.air

    m4 -DVARSIGMA=$SIGMA -DVARDSM=$DSM -DVARFLUX2=$FLUX2 0/A.air.m4 > 0/A.air
    m4 -DVARSIGMA=$SIGMA -DVARDSM=$DSM -DVARFLUX0=$FLUX0 0/N.air.m4 > 0/N.air

else

    NGROUPS=$4
    NSECTIONSPERGROUP=$5

    if [[ ! "$NGROUPS" =~ ^[0-9]+$ ]]; then

        echo "Invalid number of groups"
        exit 1

    fi

    if [[ ! "$NSECTIONSPERGROUP" =~ ^[0-9]+$ ]]; then

        echo "Invalid number of sections per group"
        exit 1

    fi

    cp constant/phaseProperties.FPT constant/phaseProperties
    cp system/sampleFields.FPT system/sampleFields
    cp system/functions.FPT system/functions

    python3 sizeGroups.py $DSM $SIGMA $NGROUPS $NSECTIONSPERGROUP $ALPHAINLET

    rm -fr constant/FPT
    mkdir -p constant/FPT

    K=0

    A=($(cat alpha.txt))
    SUMF=($(cat sumf.txt))

    for ((I=1; I<=$NGROUPS; I++)); do

        ALPHAINLET=${A[$I-1]}
        SUMFI=${SUMF[$I-1]}

        N=$((I*NSECTIONSPERGROUP))
        D=($(cat d.txt | head -n $N | tail -n $NSECTIONSPERGROUP))
        F=($(cat f.txt | head -n $N | tail -n $NSECTIONSPERGROUP))

        PHASENAME="air$I"
        PHASEPAIR="${PHASENAME}_dispersedIn_water"

        for ((J=0; J<$NSECTIONSPERGROUP; J++)); do

            DJ=${D[$J]}
            FJ=${F[$J]}

            K=$((K+1))

            echo "f$K {dSph $DJ; value $FJ;}" \
                >> constant/FPT/sizeGroups.$PHASENAME
            m4 -DVARPHASENAME=$PHASENAME -DVARFIELD=1.0 \
                -DVARINLET=$FJ -DVARFNAME=f$K 0/f.air.m4 > 0/f$K.$PHASENAME

            m4  -DVARFFIELD=f$K.$PHASENAME \
                -DVARALPHAFLUX=alphaPhi.$PHASENAME \
                -DVARFUNCTIONNAME=outletAlphaPhiF${K} \
                -DVARPATCHNAME=outlet \
                system/templates/functionF.m4 >> system/functions

            echo f$K.$PHASENAME >> system/sampleFields

        done

        VARS="-DVARPHASENAME=$PHASENAME -DVARPHASEPAIR=$PHASEPAIR"

        echo $PHASENAME >> constant/FPT/phaseNames
        m4 $VARS constant/templates/phase.m4 >> constant/FPT/phases
        m4 $VARS constant/templates/drag.m4 >> constant/FPT/drags
        m4 $VARS constant/templates/lift.m4 >> constant/FPT/lifts
        m4 $VARS constant/templates/surfaceTension.m4 \
            >> constant/FPT/surfaceTensions
        m4 $VARS constant/templates/turbulentDispersion.m4 \
            >> constant/FPT/turbulentDispersions
        m4 $VARS constant/templates/virtualMass.m4 >> constant/FPT/virtualMasses

        sed -i 's/@include/#include/g' constant/FPT/phases

        m4 -DVARPHASENAME=$PHASENAME -DVARINLET=$ALPHAINLET 0/alpha.air.m4 \
            > 0/alpha.$PHASENAME
        m4 -DVARPHASENAME=$PHASENAME  0/U.air.m4 > 0/U.$PHASENAME
        m4 -DVARPHASENAME=$PHASENAME  0/T.air.m4 > 0/T.$PHASENAME

        m4 -DVARPHASENAME=$PHASENAME -DVARFIELD=1.0 -DVARINLET=1.0 \
            -DVARFNAME=f 0/f.air.m4 > 0/f.$PHASENAME

        cp constant/momentumTransport.air constant/momentumTransport.$PHASENAME
        cp constant/thermophysicalProperties.air \
            constant/thermophysicalProperties.$PHASENAME

        m4  -DVARFLUX=phi.$PHASENAME \
            -DVARFUNCTIONNAME=outletPhi${I} \
            -DVARPATCHNAME=outlet \
            system/templates/functionPhi.m4 >> system/functions

        echo U.$PHASENAME >> system/sampleFields


    done

    echo water >> constant/FPT/phaseNames

fi

rm -f 0/*.m4

runApplication blockMesh

if [ "$MODE" == "logmom" ]; then

    runApplication setLogNormal air $SIGMA $DSM

fi

runApplication decomposePar
