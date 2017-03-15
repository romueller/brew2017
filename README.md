# On improving the quality of OTU assignments by generalising Swarm's fastidious clustering approach (BREW 2017)

Scripts for the evaluation of GeFaST's extended fastidious-clustering capabilities.

The analysis of the even resp. uneven data set is launched via targets in Makefile. 
See the makefile for more details on how to configure (mandatory) and run the analyses.

The whole analysis from FASTA file to metric values and plots is organised via the makefile.
Only testing the significance of differences is handled separately (see check-significance.R in the scripts folder). 


## Required software
 * [GeFaST](https://github.com/romueller/gefast) (version 0.7.0)
 * [Swarm](https://github.com/torognes/swarm) (version 1.2.3 and version 2.1.12)
 * [Usearch](http://www.drive5.com/usearch/download.html) (version 7)
 * Python (version 2.7 or higher)
 * Perl
 * make
 * R (with packages ggplot2 and reshape2)
 * Common Unix tools (awk, bzip2, grep, sed, sha1sum, sort, wget)
