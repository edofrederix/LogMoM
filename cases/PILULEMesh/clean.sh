#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -rf 0 polyMesh inlet system/blockMeshDict
