#!/usr/bin/perl
#
# dehanitize.pl
#
# A script to conert Han's block coordinate nomenclature into a 1 line per entry system
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';


# Han represents blocks like this:

# 01a1,Chr1,1
# 01a2,Chr1,205500
# 01b1,Chr1,32500
# 01b2,Chr1,86500

# I.e. duplicated block 01a has a start (01a1) and end (01a2), each of which has a 
# chromosomal coordinate. In turn block 01a has a region that is triplicated. This is
# 01b. 

# rather than use two lines for each block could instead just do more of a GFF
# style representation

# Chr1	BLOCK0001	01a	1	205500	2x	.
# Chr1	BLOCK0002	01b	32500	86500 3x	"Contained in BLOCK0001"

my $usage = "Usage: $0 <duplicated blocks file from Han";
die $usage unless (@ARGV == 1);

my ($input) = @ARGV;

# will want to store some block details for later
my %han_id_to_block_id;
my $block_number = "0";

open(my $in, "<", "$input") or die "Can't read from $input\n";

# read lines in pairs
while(my $first_line = <$in>){

	# skip potential header line
	next if ($first_line =~ m/^Block/);
	
	# skip unbalanced lines
#	next if $first_line =~ m/03c/;
	next if $first_line =~ m/21a1/;

	my $second_line = <$in>;
	chomp($first_line, $second_line);
	
#	print "\n1) $first_line\n2) $second_line\n";
	
	# split two lines into separate, unique, variables
	my ($han_id, $chr, $start) = split(/\t/, $first_line);
	my (undef, undef, $end)    = split(/\t/, $second_line);
	

	# trim trailing digit from Han's ID
	$han_id =~ s/[12]$//;

	# remove underscores from Han's ID (if present)
	$han_id =~ s/_//g;
	
	# Will want to store a block ID which includes
	# a block number that increases for each new unique block
	my $block_id;
	$block_number++;
	$block_id = sprintf("BLOCK%04s", $block_number);

	$han_id_to_block_id{$han_id} = $block_id;
	
	# will want to store copy number of block (2x vs 3x etc)
	my $copy;

	# will also append a comment, to indicate which parent block a 3x block has
	my $comment = ".";


	# is this a duplicated block (ends in 'a')?
	if ($han_id =~ m/a$/){
		$copy = "2x";
	} else{
		$copy = "3x";
		
		# work out parent block ID and add to comment
		my $trimmed_han_id = $han_id;
		$trimmed_han_id =~ s/[b-z]$/a/;
		
		$comment = "Contained in $han_id_to_block_id{$trimmed_han_id}";
	}
	
	
	print "$chr\t$block_id\t$han_id\t$start\t$end\t$copy\t$comment";

	print "\n";
}

close($in);


exit(0);