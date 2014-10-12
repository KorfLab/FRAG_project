#!/usr/bin/perl
#
# chromatin_state2gff.pl
#
# A script to convert simple BED file to desired GFF format
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';


# add an ID for last column
my $id = 1;

foreach my $file_number (1..9){
	open(my $in, "<", "state$file_number.txt") or die "Can't open file\n";
	while(my $line = <$in>){
		next if ($line !~ m/^\d/);

		my ($chr, $start, $end) = split(/\s+/, $line);
		# Chr1    PRICE   chromosome_breakpoint   5310181 5310181 .       -       .       ID=breakpoint0018;Parent=block0009;Name=03c2;Note="paired with breakpoint0089"
		my $formatted_id = sprintf("%05s", $id);
		print "Chr${chr}\tSequeira-Mendes_2014_paper\topen_chromatin_state\t$start\t$end\t.\t.\t.\tID=chromatin_region_$formatted_id;Note=\"state$file_number\"\n";
		$id++;
	}
	close($in);
}
__END__
