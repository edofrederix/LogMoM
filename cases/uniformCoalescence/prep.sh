#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

# Case type (S1, S2, S4, CR, GSR or FMR)

CASE=S1

# Number of quadrature points (5, 10 or 20)

NGH=20

# ------------------------------------------------------------------------------

ALPHA=0.1
SIGMA=1.0
DSM=1e-3

MATH="import math as m"

# The zeroth, second and third moment

N0=$(echo "$MATH; print($ALPHA*6.0/m.pi/$DSM**3.0*m.exp(3.0*$SIGMA**2.0))" | python3)
S0=$(echo "$MATH; print($ALPHA*6.0/m.pi/$DSM)" | python3)
A0=$(echo "$MATH; print($ALPHA*6.0/m.pi)" | python3)

# Set the case

case $CASE in

    S1)
        COARATE1=$(echo "print(1.0/$N0)" | python3)
        COARATE2=""
        COARATE3=""
        COARATE4=""

        P1=0
        P2=""
        P3=""
        P4=""
        Q1=0
        Q2=""
        Q3=""
        Q4=""

        CASENUM=1

        ;;
    S2)
        COARATE1=$(echo "print(1.0/$A0)" | python3)
        COARATE2=""
        COARATE3=""
        COARATE4=""

        P1=0
        P2=""
        P3=""
        P4=""
        Q1=3
        Q2=""
        Q3=""
        Q4=""

        CASENUM=2

        ;;
    S4)
        COARATE1=$(echo "print(1.0/$N0)" | python3)
        COARATE2=$(echo "print(1.0/$A0)" | python3)
        COARATE3=""
        COARATE4=""

        P1=0
        P2=0
        P3=""
        P4=""
        Q1=0
        Q2=3
        Q3=""
        Q4=""

        CASENUM=4

        ;;

    CR)
        COARATE1=$(echo "print(1.0/$N0)" | python3)
        COARATE2=$(echo "print(1.0/$N0)" | python3)
        COARATE3=""
        COARATE4=""

        P1=0
        P2=1
        P3=""
        P4=""
        Q1=0
        Q2=-1
        Q3=""
        Q4=""

        CASENUM=5

        ;;

    GSR)
        COARATE1=$(echo "print(1.0/$N0)" | python3)
        COARATE2=$(echo "print(1.0/$N0)" | python3)
        COARATE3=$(echo "print(1.0/$N0*($A0/$N0)**(1.0/3.0))" | python3)
        COARATE4=$(echo "print(1.0/$N0*($A0/$N0)**(1.0/3.0))" | python3)

        P1=0
        P2=1
        P3=0
        P4=1
        Q1=0
        Q2=-1
        Q3=-1
        Q4=-2

        CASENUM=6

        ;;

    FMR)
        COARATE1=$(echo "print(1.0/($N0**(5.0/6.0)*$A0**(1.0/6.0)))" | python3)
        COARATE2=$(echo "print(2.0/($N0**(5.0/6.0)*$A0**(1.0/6.0)))" | python3)
        COARATE3=$(echo "print(1.0/($N0**(5.0/6.0)*$A0**(1.0/6.0)))" | python3)
        COARATE4=""

        P1=0.5
        P2=1.0
        P3=2.0
        P4=""
        Q1=0
        Q2=-0.5
        Q3=-1.5
        Q4=""

        CASENUM=7

        ;;
    *)
        echo "Invalid case (should be S1, S2, S4, CR, GSR or FMR)"
        exit
        ;;
esac

# Set the number of quadrature nodes

case $NGH in

    5|10|20)
        ;;
    *)
        echo "Invalid number of Gauss-Hermite quadrature nodes (should be 5, 10 or 20)"
        exit
        ;;
esac

echo $N0 > properties.txt
echo $S0 >> properties.txt
echo $A0 >> properties.txt
echo $CASENUM >> properties.txt

VARS="\
    -DVARALPHA=$ALPHA \
    -DVARCOARATE1=$COARATE1 \
    -DVARCOARATE2=$COARATE2 \
    -DVARCOARATE3=$COARATE3 \
    -DVARCOARATE4=$COARATE4 \
    -DVARP1=$P1 \
    -DVARP2=$P2 \
    -DVARP3=$P3 \
    -DVARP4=$P4 \
    -DVARQ1=$Q1 \
    -DVARQ2=$Q2 \
    -DVARQ3=$Q3 \
    -DVARQ4=$Q4 \
    -DVARNGH=$NGH \
    "

find -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    m4 $VARS $IN > $OUT

done

runApplication blockMesh

cp -r 0.orig 0

rm 0/*.m4

runApplication setLogNormal bubbles $SIGMA $DSM
