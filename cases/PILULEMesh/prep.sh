#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

MESH=$1         # Number of cells per pipe radius

D1=0.034        # Diameter of the fluid pipe minus twice the gap size
D2=0.012        # Diameter of the cross pipe
S=0.001         # Gap size
F=0.5           # Ratio of the center block size and the pipe diameter
H=0.4           # Ratio of the center block vertex and the pipe diameter
THETA0=45       # Angle of the corner points at the zero plane
LI=0.5          # Inlet length
LO=0.25         # Outlet length
GL=5            # Inlet grading

##

if [[ ! "$MESH" =~ ^[1-9][0-9]+?$ ]]; then

    echo "Invalid mesh size (should be an integer)"
    exit

fi

MATH="import math as m"

R1=$(echo "print($D1/2)" | python)
R2=$(echo "print($D2/2)" | python)
R3=$(echo "print($R1+$S)" | python)

Z2=$(echo "print($D1/2*0.85)" | python)
Z3=$LI
Z32=$(echo "print($Z3-$Z2)" | python)

Z2_NO=-$LI
Z3_NO=$LO

FR1=$(echo "print($F*$R1)" | python)
FR2=$(echo "print($F*$R2)" | python)
FR3=$(echo "print($F*$R3)" | python)
HR1=$(echo "print($H*$R1)" | python)
HR2=$(echo "print($H*$R2)" | python)
HR3=$(echo "print($H*$R3)" | python)

R1MR2=$(echo "print(($R1**2 - $R2**2)**(1/2))" | python)
R3MR2=$(echo "print(($R3**2 - $R2**2)**(1/2))" | python)
R1MFR2=$(echo "print(($R1**2 - $FR2**2)**(1/2))" | python)
R3MFR2=$(echo "print(($R3**2 - $FR2**2)**(1/2))" | python)
RN=$(echo "print(0.33*$R1+0.67*$R2)" | python)
RM=$(echo "print(0.4*$R1+0.6*$R2)" | python)

R1BYSQRT2=$(echo "print($R1/2**0.5)" | python)
R2BYSQRT2=$(echo "print($R2/2**0.5)" | python)
R3BYSQRT2=$(echo "print($R3/2**0.5)" | python)

R1SINPI9=$(echo "$MATH; print($R1*m.sin(m.pi/9))" | python)
R1COSPI9=$(echo "$MATH; print($R1*m.cos(m.pi/9))" | python)

R1SINPI2O9=$(echo "$MATH; print($R1*m.sin(2*m.pi/9))" | python)
R1COSPI2O9=$(echo "$MATH; print($R1*m.cos(2*m.pi/9))" | python)

R1SINTHETA0=$(echo "$MATH; print($R1*m.sin($THETA0/180*m.pi))" | python)
R1COSTHETA0=$(echo "$MATH; print($R1*m.cos($THETA0/180*m.pi))" | python)

R3SINPI9=$(echo "$MATH; print($R3*m.sin(m.pi/9))" | python)
R3COSPI9=$(echo "$MATH; print($R3*m.cos(m.pi/9))" | python)

R3SINPI2O9=$(echo "$MATH; print($R3*m.sin(2*m.pi/9))" | python)
R3COSPI2O9=$(echo "$MATH; print($R3*m.cos(2*m.pi/9))" | python)

R3SINTHETA0=$(echo "$MATH; print($R3*m.sin($THETA0/180*m.pi))" | python)
R3COSTHETA0=$(echo "$MATH; print($R3*m.cos($THETA0/180*m.pi))" | python)

THETA50=$(echo "$MATH; print($THETA0 - m.atan($R2/$R1MR2)/m.pi*180)" | python)
THETA65=$(echo "print(90-$THETA0)" | python)

NR1=$(echo "$MATH; print(int(round(1.8*($R1-$FR1)/$R1*$MESH)))" | python)
NB=$(echo "$MATH; print(int(round($Z2/$R1*$MESH)))" | python)
NBH=$(echo "$MATH; print(int(round($NB/2)))" | python)
NC=$(echo "$MATH; print(int(round(($Z2-$R2)/$R1*$MESH)))" | python)
NS=$(echo "$MATH; print(max(int(round($S/($R1-$FR1)*$NR1)),1))" | python)

echo $NS

Q=$(echo "$MATH; print(($Z32-$GL*$R1/$MESH)/($Z32-$R1/$MESH))" | python)
NL=$(echo "$MATH; print(int(round(m.log(1.0/$GL)/m.log($Q)+1)))" | python)

NL_NO=$(echo "print(int(round(($LI+$LO)/$R1*$MESH/$GL)))" | python)

vars () {
    echo \
        -DVARR1=$R1 \
        -DVARR2=$R2 \
        -DVARR3=$R3 \
        -DVARZ2=$Z2 \
        -DVARZ3=$Z3 \
        -DVARFR1=$FR1 \
        -DVARFR2=$FR2 \
        -DVARFR3=$FR3 \
        -DVARHR1=$HR1 \
        -DVARHR2=$HR2 \
        -DVARHR3=$HR3 \
        -DVARR1MR2=$R1MR2 \
        -DVARR3MR2=$R3MR2 \
        -DVARR1MFR2=$R1MFR2 \
        -DVARR3MFR2=$R3MFR2 \
        -DVARRN=$RN \
        -DVARRM=$RM \
        -DVARR1BYSQRT2=$R1BYSQRT2 \
        -DVARR2BYSQRT2=$R2BYSQRT2 \
        -DVARR3BYSQRT2=$R3BYSQRT2 \
        -DVARR1SINPI9=$R1SINPI9 \
        -DVARR1COSPI9=$R1COSPI9 \
        -DVARR1SINPI2O9=$R1SINPI2O9 \
        -DVARR1COSPI2O9=$R1COSPI2O9 \
        -DVARR1SINTHETA0=$R1SINTHETA0 \
        -DVARR1COSTHETA0=$R1COSTHETA0 \
        -DVARR3SINPI9=$R3SINPI9 \
        -DVARR3COSPI9=$R3COSPI9 \
        -DVARR3SINPI2O9=$R3SINPI2O9 \
        -DVARR3COSPI2O9=$R3COSPI2O9 \
        -DVARR3SINTHETA0=$R3SINTHETA0 \
        -DVARR3COSTHETA0=$R3COSTHETA0 \
        -DVARTHETA50=$THETA50 \
        -DVARTHETA65=$THETA65 \
        -DVARNR1=$NR1 \
        -DVARNB=$NB \
        -DVARNBH=$NBH \
        -DVARNC=$NC \
        -DVARNS=$NS \
        -DVARNL=$NL \
        -DVARGL=$GL \
        -DVARZ2_NO=$Z2_NO \
        -DVARZ3_NO=$Z3_NO \
        -DVARNL_NO=$NL_NO
}

# Create inlet mesh

m4 $(vars) system/blockMeshDict.m4 > system/blockMeshDict

runApplication -overwrite -suffix 1 blockMesh

runApplication -overwrite -suffix 1 \
    mirrorMesh -overwrite -dict system/mirrorMeshDict.y
runApplication -overwrite -suffix 2 \
    mirrorMesh -overwrite -dict system/mirrorMeshDict.x

runApplication -overwrite -suffix 1 transformPoints Ry=180

rm -fr inlet; mkdir inlet
cp -r system constant inlet

# Create outlet mesh

GL=$(echo "print($GL*$LO/$LI)" | python)
Z3=$LO
Z32=$(echo "print($Z3-$Z2)" | python)
Q=$(echo "$MATH; print(($Z32-$GL*$R1/$MESH)/($Z32-$R1/$MESH))" | python)
NL=$(echo "$MATH; print(int(round(m.log(1.0/$GL)/m.log($Q)+1)))" | python)

m4 $(vars) system/blockMeshDict.m4 > system/blockMeshDict

runApplication -overwrite -suffix 2 blockMesh

runApplication -overwrite -suffix 1 \
    mirrorMesh -overwrite -dict system/mirrorMeshDict.y
runApplication -overwrite -suffix 2 \
    mirrorMesh -overwrite -dict system/mirrorMeshDict.x

sed -i 's/inlet/outlet/g' constant/polyMesh/boundary
sed -i 's/defaultFaces/defaultFaces2/g' constant/polyMesh/boundary

# Merge

runApplication -overwrite -suffix 1 \
    mergeMeshes -overwrite -addCases '("inlet")'

runApplication -overwrite -suffix 1 \
    stitchMesh -overwrite '((defaultFaces defaultFaces2))'

# Pack

echo "Packing"
cp -r constant/polyMesh .
tar czf mesh$MESH.tar.gz polyMesh

rm -r polyMesh

# Create quarter symmetry mesh

runApplication -overwrite -suffix 1 topoSet
runApplication -overwrite -suffix 1 subsetMesh -overwrite extraction
sed -i 's/oldInternalFaces/symm/g' constant/polyMesh/boundary
sed -i 's/internal/symmetry/g' constant/polyMesh/boundary

# Pack

echo "Packing"
cp -r constant/polyMesh .
tar czf mesh${MESH}S.tar.gz polyMesh

rm -r polyMesh

# Create mesh without obstacle

m4 $(vars) system/blockMeshDictNoObstacle.m4 > system/blockMeshDict

runApplication -overwrite -suffix 3 blockMesh

sed -i 's/defaultFaces/symm/g' constant/polyMesh/boundary
sed -i 's/empty/symmetry/g' constant/polyMesh/boundary

# Pack

echo "Packing"
cp -r constant/polyMesh .
tar czf mesh${MESH}SN.tar.gz polyMesh

rm -r polyMesh

# Create full mesh without obstacle

runApplication -overwrite -suffix 3 \
    mirrorMesh -overwrite -dict system/mirrorMeshDict.y
runApplication -overwrite -suffix 4 \
    mirrorMesh -overwrite -dict system/mirrorMeshDict.x

# Pack

echo "Packing"
cp -r constant/polyMesh .
tar czf mesh${MESH}N.tar.gz polyMesh

rm -r polyMesh

##

echo Done
