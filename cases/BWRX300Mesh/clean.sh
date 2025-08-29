#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

wclean deform

rm -fr BWRX300.msh2 polyMesh 0
