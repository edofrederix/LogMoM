#!/bin/sh

# Run from this directory
cd "${0%/*}" || exit 1

# Source run functions
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# Prepare case
cp -r 0.org 0
runApplication blockMesh
runApplication topoSet
runApplication decomposePar
