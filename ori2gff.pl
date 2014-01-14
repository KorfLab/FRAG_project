#!/usr/bin/perl
#
# ori2gff.pl
#
# A script to take origin of replication data from Nature paper and convert to GFF
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings;

my $usage = "Usage: $0 <file with Ori data>\n";
die $usage unless @ARGV == 1;
# Raw data looks like this
#Supplementary Table 1.
#N       ORIGIN NAME     CHR     START   END     SIZE    MIDPOINT        %G+C (MIDPOINT+/-100bp) GENOMIC FEATURE TYPE OF GENOMIC FEATURE DISTANCE TO GENOMIC FEATURE     CLOSEST GENOMIC FEATURE SEQUENCE (MIDPOINT+/-100bp)
#1       ori1-0010       chr1    48318   49144   827     48731   43,28   AT1G01090       gene     ---     ---    CCTTTGCTGAGGGCATGGACATGGTCACGGTAGGTACTAACGACAGAGTCAGACTTGGTAAGGAGCTTGATAAAGCCAGTAGAAACAGCCTCTTGGCCATTGTACAAGTGAACAAAACCAAACATCTTGCCTCGGTAATACATTTGAGCACACATGTCTTCGAAAGATCTACCTAGTATCATATCTTCATACAACTCCAAT
#2       ori1-0020       chr1    117311  118096  786     117704  52,74   AT1G01300       gene     ---     ---    GACCCGATATTCGACCCGAGAAAATCCAAGACCTATGCCACAATCCCCTGTTCTTCACCTCACTGCCGCCGATTAGACTCCGCCGGATGCAACACCCGTCGTAAGACTTGTCTCTACCAAGTCTCTTACGGAGATGGTTCTTTCACCGTCGGCGATTTCTCCACCGAAACGCTAACTTTCCGGCGAAATCGCGTTAAAGGC

# assume single file is specified as input
while(<>){
	next if m/^Supplementary/;
	next if m/^N/;
	
	my @data = split;
	print "$data[2]\t$data[1]\torigin_of_replication\t$data[3]\t$data[4]\t";
	print ".\t+\t.\t$data[8]\n";

}