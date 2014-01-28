#!/usr/bin/perl
#
# find_gff_overlaps.pl
#
# A simple script to find overlaps between a specified pair of coordinates and a named
# GFF file. This will ultimately be combined and updated with find_overlaps.pl
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use Getopt::Long;

######################################################
#              Command-line options                  #
######################################################

my ($chromosome, $start, $end, $help);

my $usage = "
$0 --chr <chr number> --start <start coord> --end <end coord> <gff file>

optional arguments:
        --help
";

GetOptions (
        "start=i" => \$start,
        "end=i"   => \$end,
        "chr=s"   => \$chromosome,
        "help"    => \$help,
);

die $usage unless @ARGV == 1;
my ($gff) = @ARGV;

open(my $in, "<", "$gff") or die "Can't read from $gff\n";
while(my $line = <$in>){
	my ($chr, undef, undef, $s, $e, undef, undef, undef, undef) = split(/\t/, $line);
	
	# only want to consider matches to target chromosome 
	next unless $chr eq "Chr$chromosome";
	
	#     |------------------------|
	#           s|---|e

	# do they overlap?
	if ($s >= $start and $s <= $end or
		$e >= $start and $e <= $end){
		print $line;
	}
	
}
close($in);


exit(0);




