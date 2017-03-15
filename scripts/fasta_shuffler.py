#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Reproduced from:
# Mah� F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8

"""
    Parse a fasta file and output the entries in a shuffled order.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.de>"
__date__ = "2012/05/14"
__version__ = "$Revision: 2.0"

import os
import sys
from random import shuffle
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser

#**********************************************************************#
#                                                                      #
#                            Functions                                 #
#                                                                      #
#**********************************************************************#

def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Parse a fasta file and output the entries in a shuffled order."""
    
    parser = OptionParser(
        usage="usage: %prog --input_file FILENAME --number_of_shufflings INT",
        description = desc,
        version = "%prog version 2.0")

    parser.add_option(
        "-i", "--input_file",
        metavar = "<FILENAME>",
        action = "store",
        dest = "input_file",
        help = "set <FILENAME> as input file.")

    parser.add_option(
        "-n", "--number_of_shufflings",
        metavar = "<INT>",
        action = "store",
        dest = "number_of_shuffles",
        type = "int",
        default = 10,
        help = "create <INT> shuffled outputs (default is 10).")

    (options, args) = parser.parse_args()
    
    return options.input_file, options.number_of_shuffles

#**********************************************************************#
#                                                                      #
#                              Body                                    #
#                                                                      #
#**********************************************************************#

if __name__ == '__main__':
    
    # Parse command line options.
    input_file, number_of_shuffles = option_parse()
    basename = os.path.splitext(os.path.abspath(input_file))[0]
    extension = os.path.splitext(input_file)[1]

    # Build a list of fasta entries
    input_format = "fasta"
    with open(input_file, "rU") as input_file:
        records = SeqIO.parse(input_file, input_format)
        records_list = [(record.description, str(record.seq))
                        for record in records]

    # Is this a real fasta file? (this test needs improvement)
    if len(records_list) < 2:
        print("Error: something is wrong with the input file.", file=sys.stderr)
        sys.exit(-1)

    # Create and shuffle a list of indices, and output the fasta
    # entries using that list, repeat n times
    index = range(0, len(records_list))
    zeros = len(str(number_of_shuffles))
    for n in range(0, number_of_shuffles):
        shuffle(index)
        output_file = basename + "_shuffled_" + str(n).zfill(zeros) + extension
        with open(output_file, "w") as output_file:
            for i in index:
                print(">", records_list[i][0], "\n", records_list[i][1], sep="",
                      file=output_file)
        
sys.exit(0)
