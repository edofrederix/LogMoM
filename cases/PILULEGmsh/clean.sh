#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -fr PILULE.msh2 polyMesh 0 deform/deform
