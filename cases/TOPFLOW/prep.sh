
#!/bin/bash

# Parameters -------------------------------------------------------------------

# Simulation type (either SECTIONAL or LOGMOM)

SIMTYPE=LOGMOM

# Case setup (1 -> Jg = 0.036 m/s, 2 -> Jg = 0.22 m/s, see Prasser et al.
# (2007))

SETUP=2

# Length of the pipe scaled scaled by diameter (i.e., L/D). Setup 1 works with
# 2.5, 8, 12.7, 23.2 or 40. Setup 2 only works with L/D = 40

LENGTH=40

# Number of classes (only relevant for the sectional simulation type)

CLASSES=30

# ------------------------------------------------------------------------------

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

if [[ $SETUP == 1 ]]; then

    SIGMA=0.16;
    DSM=6.8E-3;

    NX=45;

    case $LENGTH in

        2.5)
            L=0.4875
            T=4
            NZ=50
            UAIR=1.137
            ;;
        8)
            L=1.56
            T=6
            NZ=150
            UAIR=1.052
            ;;
        12.7)
            L=2.4765
            T=7
            NZ=250
            UAIR=0.988
            ;;
        23.2)
            L=4.524
            T=10
            NZ=450
            UAIR=0.871
            ;;
        40)
            L=7.8
            T=12
            NZ=750
            UAIR=0.732
            ;;
        *)
            echo "Invalid length (2.5, 8, 12.7, 23.2 or 40)"
            exit
            ;;
    esac

fi

if [[ $SETUP == 2 ]]; then

    SIGMA=1.024;
    DSM=21.35E-3;

    NX=30;

    L=7.8
    T=12
    NZ=500
    UAIR=1

fi

if [[ $SETUP != 1 && $SETUP != 2 ]]; then

    echo "Invalid setup (1 or 2)"
    exit

fi

if [[ $CLASSES -gt 100 || $CLASSES -lt 3 ]]; then

    echo "Invalid number of classes (3 =< CLASSES < 100)"
    exit

fi

if [[ $SIMTYPE != "SECTIONAL" && $SIMTYPE != "LOGMOM" ]]; then

    echo "Invalid simulation type (LOGMOM or SECTIONAL)"
    exit

fi

cp constant/phaseProperties.$SIMTYPE constant/phaseProperties
cp system/sampleFields.$SIMTYPE system/sampleFields

cp -r 0.orig.$SETUP 0

VARS="\
    -DVARNX=$NX \
    -DVARL=$L \
    -DVART=$T \
    -DVARNZ=$NZ \
    -DVARUAIR=$UAIR \
    -DVARSIGMA=$SIGMA \
    -DVARDSM=$DSM \
    "

find -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    m4 $VARS $IN > $OUT

done

rm -f 0/*.m4

runApplication blockMesh

runApplication setFields

if [[ $SIMTYPE == "SECTIONAL" ]]; then

    python sizeGroups.py $CLASSES $SETUP

fi

if [[ $SIMTYPE == "LOGMOM" ]]; then

    runApplication setLogNormal air $SIGMA $DSM

fi
