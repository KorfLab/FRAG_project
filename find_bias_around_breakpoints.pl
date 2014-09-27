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

###################################################### 
#              Command-line options                  #
###################################################### 

my ($breakpoint_gff, $feature_gff,  $help, $verbose, $mode) = (undef, undef, undef, undef, undef);
my ($range, $bin_size, $step_size)  = (50000, 2500, 500);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <GFF file of target feature>

Optional arguments:
--mode <which edge of block to look at, left or right>
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
	"help"             => \$help,
	"mode=s"           => \$mode,
	"verbose"          => \$verbose,
);


die $usage if ($help);
die $usage if (not defined $breakpoint_gff);
die $usage if (not defined $feature_gff);
die $usage if (not defined $mode);


####################################
# Set up chromosome data
####################################

# Need a hash to store sizes of chromosomes 
my %chr_sizes;
$chr_sizes{'Chr1'} = 30427671;
$chr_sizes{'Chr1'} = 28315915; # a fudge to deal with tail swap, coordinate from Han
$chr_sizes{'Chr2'} = 19698289;
$chr_sizes{'Chr3'} = 23459830;
$chr_sizes{'Chr4'} = 18585056;
$chr_sizes{'Chr5'} = 26975502;

# where does the tail swap begin on Chr4?
my $pre_tailswap_length = 16541500; # coordinate from Han

# will end up representing each chromosome as a string of dashes
my %chr_seqs;


####################################
#
# Main loop
#
####################################

print "Start\tEnd\tBreakpoint_region_bp\tNon_breakpoint_region_bp\t";
print "Feature_bp_inside\t%Inside\t";
print "Feature_bp_outside\t%Outside\n";

# will want to store actual ratios from real data in a hash
# key will be feature, value will be ratio. Do similar thing for differences
my %main_results;

# we will need to generate a fake 
# string corresponding to the length of each chromosome
initalize_virtual_chromosome_sequences();

# read main breakpoint data and feature data and store in hashes
my %breakpoints;
read_breakpoint_data();

my %features;
read_feature_data(); 

##########################################
# Main loop over different positions around breakpoint
##########################################

my ($min, $max,) = (-$range, $range);

for(my $i = $min; $i + $bin_size <= $max; $i += $step_size){
	
	my ($s, $e) = ($i, $i + $bin_size);
#		print "$s $e\n"; 

	# create copies of virtual chromosome sequences
	# one to count breakpoint bp before masking with features
	# one to count breakpoint bp after masking with feature
	my %tmp_chr_seqs_1 = %chr_seqs;
	my %tmp_chr_seqs_2 = %chr_seqs;

	# now need to mask this coordinate range around each breakpoint
	EDGE: foreach my $edge (sort {$a <=> $b} keys %breakpoints){
	
		my $bp_chr = $breakpoints{$edge}{'chr'};
		my $left   = $breakpoints{$edge}{'left'};
		my $right  = $breakpoints{$edge}{'right'};

		# may not be able to deal with ends of some blocks if they are
		# first or last on chromosome
		next EDGE if ($mode eq 'left'  and ($left  - abs($s) < 1));
		next EDGE if ($mode eq 'right' and ($right + abs($e) > $chr_sizes{$bp_chr}));

#		print "$edge\t$bp_chr\t$left\t$right\t";
		
		# mask where breakpoint regions are in chromosome
		if ($mode eq 'left'){
			substr($tmp_chr_seqs_1{$bp_chr}, $left  + $s, $bin_size) = ("B" x $bin_size);
			substr($tmp_chr_seqs_2{$bp_chr}, $left  + $s, $bin_size) = ("B" x $bin_size);
#			print $left + $s, "\t", $left + $s + $bin_size, "\t";
		} else {
			substr($tmp_chr_seqs_1{$bp_chr}, $right + $s, $bin_size) = ("B" x $bin_size);
			substr($tmp_chr_seqs_2{$bp_chr}, $right + $s, $bin_size) = ("B" x $bin_size);
#			print $right  + $s, "\t", $right + $s + $bin_size, "\t";

		}
#		print "\n";
	}	
	
	# now mask features against breakpoint regions that have already been masked
	
	foreach my $feat (sort {$a <=> $b} keys %features){ 

		my $chr = $features{$feat}{'chr'};
		my $s   = $features{$feat}{'start'};
		my $e   = $features{$feat}{'end'};

		# mask where feature occurs in virtual chromosome sequence
		my $length = $e - $s + 1;
		substr($tmp_chr_seqs_2{$chr}, $s, $length) = ("o" x $length); 
 	}

	# now compare patterns in original virtual sequence (just breakpoints)
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

	my $tmp_ratio = sprintf("%.4f", $percent_overlapping_breakpoint_regions / $percent_overlapping_non_breakpoint_regions);
	$percent_overlapping_breakpoint_regions     = sprintf("%.2f", $percent_overlapping_breakpoint_regions);
	$percent_overlapping_non_breakpoint_regions = sprintf("%.2f", $percent_overlapping_non_breakpoint_regions);
	
	print "$s $e\t"; 
	print "$original_breakpoint_bp\t";
	print "$original_non_breakpoint_bp\t";
	print "$feature_overlapping_breakpoint_regions\t";
	print "$percent_overlapping_breakpoint_regions\t";
	print "$feature_overlapping_non_breakpoint_regions\t";
	print "$percent_overlapping_non_breakpoint_regions\t";
	print "$tmp_ratio\t";



	print "\n";

}




exit;


####################################
#
# S U B R O U T I N E S
#
####################################

# create fake chromosome sequences as strings of '-'
sub initalize_virtual_chromosome_sequences{
	foreach my $chr (qw(Chr1 Chr2 Chr3 Chr4 Chr5)){
		$chr_seqs{$chr} = '-' x $chr_sizes{$chr};
	
		# have to do something different for Chromosome 4 as we are only interested in the 
		# tail swap region. So mask out first part of chromosome with a different character
		if($chr eq 'Chr4'){
			substr($chr_seqs{$chr}, 0, $pre_tailswap_length) = ("x" x $pre_tailswap_length);		
		}
	}
}


###############################################################
# read junction coordinates and represent in virtual sequences
###############################################################

sub read_breakpoint_data{
	open (my $in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

	my $block_counter = 0;
	while(my $line = <$in>){
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

		# skip comments
		next if ($line =~ m/^#/);

		# only interested in 2x or 3x blocks
		next unless ($feature eq 'copy_number_gain');

		$block_counter++;
		
		# 3 things to store for each block
		$breakpoints{$block_counter}{'chr'}   = $chr;
		$breakpoints{$block_counter}{'left'}  = $s;
		$breakpoints{$block_counter}{'right'} = $e;
	}

	close($in);
}



###############################################################
# read feature data from input GFF file 
###############################################################

sub read_feature_data{

	my $feature_count = 0;
	
	open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

	while(my $line = <$in>){
		chomp($line);
		# skip GFF header lines
		next if ($line =~ m/^#/);
		
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

		#  skip if not chr1 or chr4?
		next unless (($chr eq 'Chr1') or ($chr eq 'Chr4'));

		# skip tailswap regions of Chr1 and Chr4
		next if ($chr eq 'Chr1' and $s > $chr_sizes{'Chr1'});
		next if ($chr eq 'Chr4' and $s < $pre_tailswap_length);

		$feature_count++;
		
		# add data to hash
		$features{$feature_count}{'chr'}   = $chr;
		$features{$feature_count}{'start'} = $s;
		$features{$feature_count}{'end'}   = $e;

	}
	close($in);
}
