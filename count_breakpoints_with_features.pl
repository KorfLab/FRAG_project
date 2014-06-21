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

###################################################### 
#              Command-line options                  #
###################################################### 

my ($breakpoint_gff, $feature_gff, $bp, $help, $verbose) = (undef, undef, 1000, undef, undef);

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



####################################
#
# Main loop
#
####################################


my @gff_features = qw(gene DNA_replication_origin DNAseI_hypersensitive_site);

foreach my $gff_feature (@gff_features){

	open (my $in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

	my $breakpoint_count = 0;
	my $breakpoints_containing_features = 0;
		
	OUTER: while(my $breakpoint_data = <$in>){
		my ($chrA, undef, $featureA, $qs, $qe, undef, undef, undef, $comment) = split(/\t/, $breakpoint_data);	

		# skip comments
		next if ($breakpoint_data =~ m/^#/);

		# want chromosome breakpoints, but ignore any which are effectively the ends of
		# the chromosomes
		next unless ($featureA eq 'chromosome_breakpoint' and $comment !~ m/telomeric end/);
		
		$breakpoint_count++;

		# now main GFF file
		open (my $in2, "<", $feature_gff) or die "Can't read $feature_gff\n";

		while(my $gff_line = <$in2>){
			my ($chrB, undef, $featureB, $ss, $se, undef, undef, undef, undef) = split(/\t/, $gff_line);	

			# skip comments
			next if ($gff_line =~ m/^#/);

			# only want to look at one feature at a time
			next unless $featureB eq $gff_feature;

			# only look at features on same chromosome as current breakpoint
			next unless ($chrB eq $chrA);

			# skip tailswap regions of Chr1 and Chr4
			next if ($chrB eq 'Chr1' and $ss > $chr_sizes{'Chr1'});
			next if ($chrB eq 'Chr4' and $se < $pre_tailswap_length);

			# can now ask where current GFF feature overlaps current breakpoint
			my $overlap = check_for_overlap($qs,$qe,$ss,$se);
			if ($overlap){
				$breakpoints_containing_features++;
				print "\tProcessed $breakpoint_count breakpoints, found $breakpoints_containing_features overlap with $gff_feature\n" if ($verbose);
				next OUTER;
			}
	
		}
		close($in2);

		print "\tProcessed $breakpoint_count breakpoints, found $breakpoints_containing_features overlap with $gff_feature\n" if ($verbose);

	}
	my $percent = sprintf("%.1f", ($breakpoints_containing_features / $breakpoint_count) * 100);
	print "FINAL: Processed $breakpoint_count breakpoint regions, found $breakpoints_containing_features ($percent%) overlap with $gff_feature\n";

	close($in);
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