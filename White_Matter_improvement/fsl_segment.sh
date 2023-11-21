#!/bin/bash
# -----------------------------------------------------------
# Bash script which applies a 3 class FAST segmentation on  
# white matter volumes from reference models
# 
# REQUIRES FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FAST)
# Note: Since FSL tools work only in linux based OS, this script 
#       has been launch from WSL2
#
# 2022-03-23 Andrés le Boeuf
# andres.le.boeuf@estudiantat.upc.edu
# -----------------------------------------------------------

if [ "$1" == "-h" ] ; then
    echo "Entrance variables: "
    echo "\$1 = output image base name"
    exit 1
fi

INPUT_IMAGE=$1

GETF0="/usr/local/fsl/bin/fast"

$GETF0 -t 2 -n 3 -H 0.1 -I 4 -l 20.0 -o $1 
