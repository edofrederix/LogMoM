#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

CASE=$1

if [[ ! "$CASE" =~ ^R[0-4]?$ ]]; then

    echo "Invalid case. Case should be R0, R1, R2, R3 or R4."
    exit

fi

NGH=5

ALPHA=0.1
SIGMA=1.0
DSM=1e-3

# ------------------------------------------------------------------------------

MATH="import math as m"

# The zeroth, second and third moment

M0=$(echo "$MATH; print($ALPHA*6.0/m.pi/$DSM**3.0*m.exp(3.0*$SIGMA**2.0))" | python3)
M2=$(echo "$MATH; print($ALPHA*6.0/m.pi/$DSM)" | python3)
M3=$(echo "$MATH; print($ALPHA*6.0/m.pi)" | python3)

# Read the case

case $CASE in

    R0)
        BREAKRATE1=1.0
        BREAKRATE2=""

        R1=0
        R2=""

        CASENUM=0

        ;;
    R1)
        BREAKRATE1=$(echo "print(($M0/$M3)**(1.0/3.0))" | python3)
        BREAKRATE2=""

        R1=1
        R2=""

        CASENUM=1

        ;;
    R2)
        BREAKRATE1=$(echo "print($M0/$M3)" | python3)
        BREAKRATE2=""

        R1=3
        R2=""

        CASENUM=2

        ;;
    R3)
        BREAKRATE1=1.0
        BREAKRATE2=$(echo "print(($M0/$M3)**(1.0/3.0))" | python3)

        R1=0
        R2=1

        CASENUM=3

        ;;
    R4)
        BREAKRATE1=1.0
        BREAKRATE2=$(echo "print($M0/$M3)" | python3)

        R1=0
        R2=3

        CASENUM=4

        ;;
    *)
        echo "Invalid case (should be R0, R1, R2, R3 or R4)"
        exit
        ;;
esac

# Read the number of quadrature nodes

case $NGH in

    5|10|20)
        ;;
    *)
        echo "Invalid number of Gauss-Hermite quadrature nodes (should be 5, 10 or 20)"
        exit
        ;;
esac

echo $M0 > properties.txt
echo $M2 >> properties.txt
echo $M3 >> properties.txt
echo $CASENUM >> properties.txt

VARS="\
    -DVARALPHA=$ALPHA \
    -DVARBREAKRATE1=$BREAKRATE1 \
    -DVARBREAKRATE2=$BREAKRATE2 \
    -DVARR1=$R1 \
    -DVARR2=$R2 \
    -DVARNGH=$NGH \
    "

find -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    m4 $VARS $IN > $OUT

done

runApplication blockMesh

cp -r 0.org 0

rm 0/*.m4

runApplication setLogNormal bubbles $SIGMA $DSM
