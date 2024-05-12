#!/bin/sh

# Run from this directory
cd "${0%/*}" || exit 1

# Source run functions
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# Run case
runParallel    $(getApplication)
runApplication reconstructPar
