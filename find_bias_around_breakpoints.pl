#!/usr/bin/perl
#
# find_bias_around_breakpoints.pl
#
# A script to take breakpoint regions and see whether certain chromosome features
# are biased in frequency in the immediate vicinity
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use FRAG;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($breakpoint_gff, $feature_gff,  $help, $verbose, $flip) = (undef, undef, undef, undef, 0);
my ($range, $bin_size, $step_size)  = (25000, 500, 50);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <GFF file of target feature>

Optional arguments:
--flip <whether to flip 'right' blocks to make results consistent for inside/outside blocks>
--range <max distance from breakpoint, default = $range>
--bin_size <interval to use, default = $bin_size>
--step_size <interval to step over, default = $step_size>
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"breakpoint_gff=s" => \$breakpoint_gff,
	"feature_gff=s"    => \$feature_gff,
	"range=i"          => \$range,
	"bin_size=i"       => \$bin_size,
	"step_size=i"      => \$step_size,
	"flip"             => \$flip,
	"help"             => \$help,
	"verbose"          => \$verbose,
);


die $usage if ($help);
die $usage if (not defined $breakpoint_gff);
die $usage if (not defined $feature_gff);

####################################
# Set up chromosome data
####################################

# Use a hash to store sizes of chromosomes 
# keys will be 'chr1', 'chr2' etc.
# want the shorter 'tailswap' length for Chr1
my %chr_sizes = FRAG::get_chromosome_lengths('tailswap'); 

# where does the tail swap begin on Chr4?
my $pre_tailswap_length = FRAG::get_chr4_pre_tailswap_length();

# read main block data to get left & right coordinates of each block
my %blocks = FRAG::read_block_data($breakpoint_gff);

# read and store coordinates of feature from a single-feature GFF file
my %features = FRAG::read_feature_data($feature_gff); 

# we will need to generate a fake string corresponding to the length of each chromosome
# will end up representing each chromosome as a string of dashes
my %chr_seqs = FRAG::initalize_virtual_chromosome_sequences();



####################################
#
# Main loop
#
####################################

print "Start\tEnd\tMidpoint\tBreakpoint_bp\tNon_breakpoint_bp\t";
print "Feature_bp_inside\t%Inside\t";
print "Feature_bp_outside\t%Outside\tRatio\n";

# will want to store actual ratios from real data in a hash
# key will be feature, value will be ratio. Do similar thing for differences
my %main_results;



##########################################
# Main loop over different positions around breakpoint
##########################################

my ($min, $max,) = (-$range, $range);

for(my $i = $min; $i + $bin_size <= $max; $i += $step_size){
	
	my ($s, $e) = ($i, $i + $bin_size);
	my $mid = int($s + (($e - $s) / 2));
	warn "Analysizing all windows at $mid bp ($s to $e) from all breakpoints\n" if ($verbose); 

	# create two copies of virtual chromosome sequences
	# both copies will be masked with breakpoint regions
	# but only 2nd copy will be additionally masked with target features
	# can then compare both sets
	my %tmp_chr_seqs_1 = %chr_seqs;
	my %tmp_chr_seqs_2 = %chr_seqs;

	# now need to mask this coordinate range around each breakpoint
	EDGE: foreach my $edge (sort {$a <=> $b} keys %blocks){
	
		my $bp_chr = $blocks{$edge}{'chr'};
		my $left   = $blocks{$edge}{'left'};
		my $right  = $blocks{$edge}{'right'};

		# may not be able to deal with ends of some blocks if they are
		# first or last on chromosome (or if $bin_size is really large)
		next EDGE if ($left  - abs($s) < 1);
		next EDGE if ($left  + abs($e) > $chr_sizes{$bp_chr});
		next EDGE if ($right - abs($s) < 1);
		next EDGE if ($right + abs($e) > $chr_sizes{$bp_chr});


		# what coordinate do we start masking breakpoint regions
		my $left_mask_start  = $left + $s;
		my $right_mask_start = $right + $s;
		my $left_mask_end    = $left_mask_start + $bin_size - 1;
		my $right_mask_end   = $right_mask_start + $bin_size - 1;

		# if we are in $flip mode, then we need to deal with the 'right' edge differently
		# e.g. position -1,000 to -900 should be changed to be +900 to +1,000
		# this ensures that we are consistent with looking at regions that are 
		# either inside or outside blocks
		if ($flip){
			$right_mask_start = $right - $e;
			$right_mask_end   = $right_mask_start + $bin_size - 1;
		} 

		# now mask where breakpoint regions are in chromosome
		substr($tmp_chr_seqs_1{$bp_chr}, $left_mask_start,  $bin_size) = ("B" x $bin_size);
		substr($tmp_chr_seqs_2{$bp_chr}, $left_mask_start,  $bin_size) = ("B" x $bin_size);
		substr($tmp_chr_seqs_1{$bp_chr}, $right_mask_start, $bin_size) = ("B" x $bin_size);
		substr($tmp_chr_seqs_2{$bp_chr}, $right_mask_start, $bin_size) = ("B" x $bin_size);						

		warn "\t$s to $e\t$edge\t$bp_chr\tBLOCK:$left-$right\t$left_mask_start-$left_mask_end\t$right_mask_start-$right_mask_end\n" if ($verbose);
	}	
	
	# now mask these sequences with genomic features (on top of masking with breakpoint regions)
	foreach my $feat (sort {$a <=> $b} keys %features){ 

		my $chr = $features{$feat}{'chr'};
		my $s   = $features{$feat}{'start'};
		my $e   = $features{$feat}{'end'};

		# mask where feature occurs in virtual chromosome sequence
		my $length = $e - $s + 1;
		substr($tmp_chr_seqs_2{$chr}, $s, $length) = ("o" x $length); 
 	}


	# now compare patterns in original virtual sequence (just blocks)
	# and tmp virtual sequence (masked with feature)
	my $original_breakpoint_bp      = $tmp_chr_seqs_1{'Chr1'} =~ tr/B/B/;
	$original_breakpoint_bp        += $tmp_chr_seqs_1{'Chr4'} =~ tr/B/B/;

	my $original_non_breakpoint_bp  = $tmp_chr_seqs_1{'Chr1'} =~ tr/-/-/;
	$original_non_breakpoint_bp    += $tmp_chr_seqs_1{'Chr4'} =~ tr/-/-/;

	my $remaining_breakpoint_bp     = $tmp_chr_seqs_2{'Chr1'} =~ tr/B/B/;
	$remaining_breakpoint_bp       += $tmp_chr_seqs_2{'Chr4'} =~ tr/B/B/;

	my $remaining_non_breakpoint_bp = $tmp_chr_seqs_2{'Chr1'} =~ tr/-/-/;
	$remaining_non_breakpoint_bp   += $tmp_chr_seqs_2{'Chr4'} =~ tr/-/-/;

	my $feature_overlapping_breakpoint_regions     = $original_breakpoint_bp - $remaining_breakpoint_bp;
	my $feature_overlapping_non_breakpoint_regions = $original_non_breakpoint_bp - $remaining_non_breakpoint_bp;

	# likely to have some zero counts for bases overlapping junction regions, so handle accordingly
	my ($percent_overlapping_breakpoint_regions, $percent_overlapping_non_breakpoint_regions);

	if ($feature_overlapping_breakpoint_regions == 0){
		$percent_overlapping_breakpoint_regions = "0";
	} else{
		$percent_overlapping_breakpoint_regions = ($feature_overlapping_breakpoint_regions / $original_breakpoint_bp) * 100;
	}

	if ($feature_overlapping_non_breakpoint_regions == 0){
		$percent_overlapping_non_breakpoint_regions = "0";
	} else{
		$percent_overlapping_non_breakpoint_regions = ($feature_overlapping_non_breakpoint_regions / $original_non_breakpoint_bp) * 100;
	}

	my $ratio = sprintf("%.4f", $percent_overlapping_breakpoint_regions / $percent_overlapping_non_breakpoint_regions);
	$percent_overlapping_breakpoint_regions     = sprintf("%.2f", $percent_overlapping_breakpoint_regions);
	$percent_overlapping_non_breakpoint_regions = sprintf("%.2f", $percent_overlapping_non_breakpoint_regions);
	
	print "$s\t$e\t$mid\t"; 
	print "$original_breakpoint_bp\t";
	print "$original_non_breakpoint_bp\t";
	print "$feature_overlapping_breakpoint_regions\t";
	print "$percent_overlapping_breakpoint_regions\t";
	print "$feature_overlapping_non_breakpoint_regions\t";
	print "$percent_overlapping_non_breakpoint_regions\t";
	print "$ratio\t";
	print "\n";

}




exit;




__END__

tmp1 -----BBBBB---------BBBBB-----BBBBB------------------BBBBB----------
tmp2 -----BBBBB---------BBBBB-----BBBBB------------------BBBBB----------

tmp1 -----BBBBB---------BBBBB-----BBBBB------------------BBBBB----------
tmp2 -----BBooB-------oooooBB-----ooooo------------------BBBBB--ooooo---