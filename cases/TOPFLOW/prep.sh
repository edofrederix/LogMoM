
#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

D=0.195         # Diameter of the pipe [m]
L=7.5855        # Length of the pipe = D*(40-1.1) in [m]
POUT=120000     # Outlet pressure [Pa]

CASE=$1
MESH=$2
MODE=$3

case $CASE in

    A)

        JL=1
        JG=0.037
        DSM=0.00114
        SIGMA=0.3062
        ALPHAIN=0.031828

    ;;

    B)

        JL=1
        JG=0.22
        DSM=0.003717
        SIGMA=0.666149
        ALPHAIN=0.10225

    ;;

    C)

        JL=1
        JG=0.53
        DSM=0.003289
        SIGMA=0.73136
        ALPHAIN=0.17311

    ;;

    *)

        echo "Invalid case specified, should be [A-C]"
        exit

    ;;

esac

if [[ ! "$MESH" =~ ^[1-9][0-9]?$ ]]; then

    echo "Invalid mesh size. Specify the number of cells on the radius"
    exit

fi

if [[ ! "$MODE" =~ ^(logmom|fpt)$ ]]; then

    echo "Invalid mode (should be logmom or fpt)"
    exit
fi

echo "Mode = $MODE, dsm = $DSM, sigma = $SIGMA"

echo $CASE > case.txt
echo $MESH > mesh.txt
echo $MODE > mode.txt
echo $DSM > dsm.txt
echo $SIGMA > sigma.txt

# ------------------------------------------------------------------------------

ALPHAWATERIN=$(echo "print(1.0-$ALPHAIN)" | python3)

UAIRIN=$(python3 U.py $POUT $L $ALPHAIN $JG)
UWATERIN=$(echo "print($JL/$ALPHAWATERIN)" | python3)

cp -r 0.org 0

if [ "$MODE" == "logmom" ]; then

    POLY=$4

    cp system/sampleFields.LogMoM system/sampleFields
    cp system/functions.LogMoM system/functions

    case $POLY in

        true)

            m4 -DVARPHINAME="phi2.air" -DVARSIGMA=$SIGMA -DVARDSM=$DSM \
                0/A.air.m4 > 0/A.air
            m4 -DVARPHINAME="phi0.air" -DVARSIGMA=$SIGMA -DVARDSM=$DSM \
                0/N.air.m4 > 0/N.air

            m4 -DVARPHASESYSTEM=basicPolyPhaseSystem \
                -DVARPHASEMODEL=pureIsothermalPolyPhaseModel \
                constant/phaseProperties.LogMoM.m4 > constant/phaseProperties

        ;;

        false)

            m4 -DVARPHINAME="phi.air" -DVARSIGMA=$SIGMA -DVARDSM=$DSM \
                0/A.air.m4 > 0/A.air
            m4 -DVARPHINAME="phi.air" -DVARSIGMA=$SIGMA -DVARDSM=$DSM \
                0/N.air.m4 > 0/N.air

            m4  -DVARPHASESYSTEM=basicMultiphaseSystem \
                -DVARPHASEMODEL=pureIsothermalPhaseModel \
                constant/phaseProperties.LogMoM.m4 > constant/phaseProperties

            sed -i 's/U0\.air//g' system/sampleFields
            sed -i 's/U2\.air//g' system/sampleFields

            sed -i 's/phi2/phi/g' system/functions
            sed -i 's/phi0/phi/g' system/functions

        ;;

        *)

            echo "Specify if simulation is poly-celeric (true or false)"
            exit
        ;;

    esac

    m4 -DVARPHASENAME=air -DVARALPHAIN=$ALPHAIN -DVARCASE=$CASE \
        0/alpha.air.m4 > 0/alpha.air
    m4 -DVARALPHAWATERIN=$ALPHAWATERIN -DVARCASE=$CASE \
        0/alpha.water.m4 > 0/alpha.water

    m4 -DVARPHASENAME=air -DVARUAIRIN=$UAIRIN 0/U.air.m4 > 0/U.air
    m4 -DVARPHASENAME=air -DVARUAIRIN=$UAIRIN 0/U0.air.m4 > 0/U0.air
    m4 -DVARPHASENAME=air -DVARUAIRIN=$UAIRIN 0/U2.air.m4 > 0/U2.air

    m4 -DVARUWATERIN=$UWATERIN 0/U.water.m4 > 0/U.water

    m4 -DVARPHASENAME=air 0/T.air.m4 > 0/T.air

else

    cp constant/phaseProperties.FPT constant/phaseProperties
    cp system/sampleFields.FPT system/sampleFields
    cp system/functions.FPT system/functions

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

    python3 sizeGroups.py $DSM $SIGMA $NGROUPS $NSECTIONSPERGROUP $ALPHAIN

    rm -fr constant/FPT
    mkdir -p constant/FPT

    K=0

    A=($(cat alpha.txt))
    SUMF=($(cat sumf.txt))

    for ((I=1; I<=$NGROUPS; I++)); do

        ALPHAIN=${A[$I-1]}
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

            m4  -DVARFFIELD=f$K.$PHASENAME \
                -DVARALPHAFLUX=alphaPhi.$PHASENAME \
                -DVARFUNCTIONNAME=inletAlphaPhiF${K} \
                -DVARPATCHNAME=inlet \
                system/templates/functionF.m4 >> system/functions

            echo f$K.$PHASENAME >> system/sampleFields

        done

        VARS="-DVARPHASENAME=$PHASENAME -DVARPHASEPAIR=$PHASEPAIR"

        echo $PHASENAME >> constant/FPT/phaseNames
        m4 $VARS constant/templates/phase.m4 >> constant/FPT/phases
        m4 $VARS constant/templates/drag.m4 >> constant/FPT/drags
        m4 $VARS constant/templates/lift.m4 >> constant/FPT/lifts
        m4 $VARS constant/templates/wallLubrication.m4 \
            >> constant/FPT/wallLubrications
        m4 $VARS constant/templates/surfaceTension.m4 \
            >> constant/FPT/surfaceTensions
        m4 $VARS constant/templates/turbulentDispersion.m4 \
            >> constant/FPT/turbulentDispersions
        m4 $VARS constant/templates/virtualMass.m4 >> constant/FPT/virtualMasses

        sed -i 's/@include/#include/g' constant/FPT/phases

        m4 -DVARPHASENAME=$PHASENAME -DVARALPHAIN=$ALPHAIN -DVARCASE=$CASE \
            0/alpha.air.m4 > 0/alpha.$PHASENAME
        m4 -DVARPHASENAME=$PHASENAME  -DVARUAIRIN=$UAIRIN 0/U.air.m4 \
            > 0/U.$PHASENAME
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

        m4  -DVARFLUX=phi.$PHASENAME \
            -DVARFUNCTIONNAME=inletPhi${I} \
            -DVARPATCHNAME=inlet \
            system/templates/functionPhi.m4 >> system/functions

        echo U.$PHASENAME >> system/sampleFields

    done

    m4 -DVARUWATERIN=$UWATERIN 0/U.water.m4 > 0/U.water
    m4 -DVARALPHAWATERIN=$ALPHAWATERIN -DVARCASE=$CASE 0/alpha.water.m4 \
        > 0/alpha.water

    echo water >> constant/FPT/phaseNames

fi

rm -f 0/*.m4

wmake -s TOPFLOWAlphaInlet

NR=$MESH
NZ=$(echo "print(int(round(2.0*$L/$D*$MESH/6.0)))" | python)

m4 -DVARNR=$NR -DVARNZ=$NZ system/blockMeshDict.m4 > system/blockMeshDict

runApplication blockMesh

if [ "$MODE" == "logmom" ]; then

    runApplication setLogNormal air $SIGMA $DSM

fi

runApplication decomposePar
