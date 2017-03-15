#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8

DATA_SET=$1

wget http://sbr2.sb-roscoff.fr/download/externe/de/fmahe/${DATA_SET}.fasta.bz2
bzip2 -d ${DATA_SET}.fasta.bz2
mv ${DATA_SET}.fasta data/