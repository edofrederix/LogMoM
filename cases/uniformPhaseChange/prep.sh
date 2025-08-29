#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

# Parameters -------------------------------------------------------------------

# Case types
#
# Case 1: constant Nusselt number
# Case 2: quadratic Nusselt number
#

CASE=${1:-1}

ALPHA=0.01
SIGMA=0.5
DSM=1e-2

MATH="import math as m"

# The zeroth, second and third moment

M3=$(echo "$MATH; print($ALPHA*6.0/m.pi)" | python3)
M2=$(echo "$MATH; print($M3/$DSM)" | python3)
M0=$(echo "$MATH; print($M3/$DSM**3.0*m.exp(3.0*$SIGMA**2.0))" | python3)

ALPHAWATER=$(echo "$MATH; print(1-$ALPHA)" | python3)

# Set the case

case $CASE in

    1)
        NUCORR=constantNu

        ;;
    2)
        NUCORR=sqrNu

        ;;

    *)
        echo "Invalid case specified (should be 1 or 2)"
        exit
        ;;

esac

NU=1
D0=$DSM

# ------------------------------------------------------------------------------

echo $CASE > properties.txt
echo $DSM >> properties.txt
echo $SIGMA >> properties.txt
echo $ALPHA >> properties.txt
echo $NU >> properties.txt
echo $D0 >> properties.txt

VARS="\
    -DVARALPHASTEAM=$ALPHA \
    -DVARALPHAWATER=$ALPHAWATER \
    -DVARNUCORR=$NUCORR \
    -DVARNU=$NU \
    -DVARD0=$D0 \
    "

find -name *.m4 | while read IN; do

    OUT=$(echo $IN | rev | cut -c 4- | rev)

    m4 $VARS $IN > $OUT

done

runApplication blockMesh

cp -r 0.org 0

rm 0/*.m4

runApplication setLogNormal steam $SIGMA $DSM
