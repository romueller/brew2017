#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8

# Inputs
INPUT_TABLE=$1
LONG_TABLE=$2

# Produce a long table version
sort -t "," -k1,1d -k2,2n -k3,3n ${INPUT_TABLE} | \
    awk 'BEGIN {
             FS = ","
             OFS = ","
             split("recall precision NMI rand adjustedrand", a, " ")
         }
         {for (i=1; i<=5; i++) print $1, $2, $3, a[i], $(i+3)}' > ${LONG_TABLE}
