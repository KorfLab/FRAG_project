#!/usr/bin/perl
#
# bed2gff.pl
#
# A script to convert simple BED file to desired GFF format
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage: $0 <filename>\n" unless (@ARGV == 1);

# print GFF3 output header

# add an ID for last column
my $id = 1;

while(<>){
	# only want DHS sites that are from WT leaf
	next unless m/wtleaf/;
	my ($chr, $start, $end) = split;
	# Chr1    PRICE   chromosome_breakpoint   5310181 5310181 .       -       .       ID=breakpoint0018;Parent=block0009;Name=03c2;Note="paired with breakpoint0089"
	my $formatted_id = sprintf("%05s", $id);
	print "$chr\tGEO\tDNAseI_hypersensitive_site\t$start\t$end\t.\t.\t.\tID=DHS_$formatted_id\n";
	$id++;
}