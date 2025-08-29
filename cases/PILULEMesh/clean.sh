#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -rf \
    system/blockMeshDict \
    system/createPatchDict \
    inlet \
    pipe
