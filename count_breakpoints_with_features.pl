#!/usr/bin/perl
#
# count_breakpoints_with_features.pl
#
# A script to find overlaps between various GFF features in A. thaliana
# and our predicted regions of duplication/triplication in our mutant line
# that exhibits chromosome shattering
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use List::Util qw(shuffle);
use Getopt::Long;
use FRAG;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($breakpoint_gff, $feature_gff, $bp, $help, $verbose) = (undef, undef, 100, undef, undef);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>

Optional arguments:
--bp <how many bp to extract from inside/outside region (default = $bp)>
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"breakpoint_gff=s" => \$breakpoint_gff,
	"feature_gff=s"    => \$feature_gff,
	"bp=i"             => \$bp,
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


# get details of breakpoints coordinates
my %breakpoints = FRAG::read_breakpoint_data($breakpoint_gff);


# get list of all target GFF features we are interested in
my @target_gff_features = sort (FRAG::get_list_of_GFF_features($feature_gff));

####################################
#
# Main loop
#
####################################

open (my $in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

# hash to keep track of which breakpoints overlap which features
my %breakpoints_containing_features;
	

BREAKPOINT: foreach my $breakpoint (sort {$a <=> $b} keys %breakpoints){
	print "Breakpoint: $breakpoint\n" if ($verbose);
	my $chrA = $breakpoints{$breakpoint}{chr};
	my $qs   = $breakpoints{$breakpoint}{pos};
	my $qe   = $breakpoints{$breakpoint}{pos};


	# now loop over feature GFF file
	open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

	OUTER: while(my $gff_line = <$in>){
		my ($chrB, undef, $feature, $ss, $se, undef, undef, undef, $comment) = split(/\t/, $gff_line);	

		# skip comments
		next if ($gff_line =~ m/^#/);

		# only look at features on same chromosome as current breakpoint
		next unless ($chrB eq $chrA);

		# skip tailswap regions of Chr1 and Chr4
		next if ($chrB eq 'Chr1' and $ss > $chr_sizes{'Chr1'});
		next if ($chrB eq 'Chr4' and $se < $pre_tailswap_length);


		#  want to separate out chromatin state data (if it exists)
		if($feature eq 'open_chromatin_state'){
			my ($state) = $comment =~ m/Note=\"state(\d+)\"/;
			$feature = "open_chromatin_state_${state}";
		}


		# only want to look features on our list
		next unless ($feature ~~ @target_gff_features);

		# can now ask where current GFF feature overlaps current breakpoint
		my $overlap = check_for_overlap($qs, $qe, $ss, $se);
		if ($overlap){

			# flag that for this combination of feature and breakpoint we have see an overlap
			$breakpoints_containing_features{$feature}{$breakpoint}++;

			my $overlap_count = keys (%{$breakpoints_containing_features{$feature}});
			print "\t$feature: $overlap_count breakpoints overlap\n" if ($verbose and $breakpoints_containing_features{$feature}{$breakpoint} <2);
			# once we know a breakpoint overlaps at least
		}
	}
}

close($in);

my $number_of_breakpoints = keys %breakpoints;

foreach my $feature (@target_gff_features){
	my ($overlap_count, $percent) = (0, 0.00);

	$overlap_count = keys (%{$breakpoints_containing_features{$feature}});
	$percent = sprintf("%.1f", ($overlap_count / $number_of_breakpoints) * 100);

	print "FINAL: $overlap_count/$number_of_breakpoints breakpoint regions ($percent%) overlap with $feature\n";

}


exit;


sub check_for_overlap{

	my ($start1, $end1, $start2, $end2) = @_;

	# modify coordinates by breakpoint by $bp
	$start1 = $start1 - ($bp / 2);
	$end1   = $end1   + ($bp / 2);

	# LHS overlap of first window (or same coords)
	if(($start1 <= $start2) and ($end1 >= $start2)){
		return(1);
	}
	# RHS overlap of first window (or same coords)
	elsif(($start1 >= $start2) and ($start1 <= $end2)){
		return(1);
	}
	# First window completely encloses second window
	elsif(($start1 <= $start2) and ($end1 >= $end2)){
		return(1);
	}
	# First window is completely enclosed by second window
	elsif(($start1 > $start2) and ($end1 < $end2)){
		return(1);
	} else{
		return(0); # no overlap
	}
}