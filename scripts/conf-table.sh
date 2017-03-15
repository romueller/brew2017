#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8

########## CONFIGURATION ##########

# Adjust to your machine

# Binaries
PYTHON="python2.7" # python command / path to binary


# Adjust only if you know what you are doing

# Scripts
BUILD_CONFUSION_TABLE="scripts/confusion_table.py"
COMPUTE_CLUSTER_METRICS="scripts/CValidate.pl"

###################################

# Inputs
CLUSTERING_RESULTS=$1
SHUFFLING=$2
METHOD=$3
THRESHOLD=$4
TAXONOMIC_ASSIGNMENTS=$5
CONFUSION_TABLE=$6
METRICS=$7

# Build confusion table; derive and store metric values
${PYTHON} "${BUILD_CONFUSION_TABLE}" -t "${TAXONOMIC_ASSIGNMENTS}" -s "${CLUSTERING_RESULTS}" > "${CONFUSION_TABLE}"
METRICS_DATA=$(perl "${COMPUTE_CLUSTER_METRICS}" --cfile="${CONFUSION_TABLE}")
echo "${METHOD},${SHUFFLING},${THRESHOLD},${METRICS_DATA}" >> "${METRICS}"

