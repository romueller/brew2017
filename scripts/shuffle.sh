#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8

PYTHON=$1
FASTA_SHUFFLER=$2
FASTA=$3
NUM_REPS=$4

${PYTHON} ${FASTA_SHUFFLER} -i ${FASTA} -n ${NUM_REPS}
cp ${FASTA} ${FASTA%.fasta}_shuffled_$(echo ${NUM_REPS} | sed 's/./0/g').fasta
