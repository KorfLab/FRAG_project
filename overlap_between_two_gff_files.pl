#!/usr/bin/perl
#
# overlap_between_two_gff_files.pl
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

my ($breakpoint_gff, $feature_gff, $bp, $shuffles, $help, $verbose) = (undef, undef, 1000, 0, undef, undef);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>

Optional arguments:
--bp <how many bp to extract from inside/outside region (default = $bp)>
--shuffles <how many shuffling iterations to perform (default = $shuffles)>
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"breakpoint_gff=s" => \$breakpoint_gff,
	"feature_gff=s"    => \$feature_gff,
	"bp=i"             => \$bp,
	"shuffles=i"       => \$shuffles,
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

# for the main run and each shuffle, we will need to generate a fake 
# string corresponding to the length of each chromosome
# will end up representing each chromosome as a string of dashes
my %chr_seqs = FRAG::initalize_virtual_chromosome_sequences();

# get details of breakpoints coordinates
my %breakpoints = FRAG::read_breakpoint_data($breakpoint_gff);



####################################
#
# Main loop
#
####################################

print "Run\tReal_ratio\tBp\tFeature\tBreakpoint_region_bp\tNon_breakpoint_region_bp\t";
print "Feature_bp_inside\t%Inside\t";
print "Feature_bp_outside\t%Outside\t";
print "Shuffled_ratio\t";
print "Above\tSame\tBelow\n";

# will want to store actual ratios from real data in a hash
# key will be feature, value will be ratio. Do similar thing for differences
my %main_results;

# need to do one master run and then a number of shuffles to test significance
# (as denoted by $shuffles)
for (my $i = 0; $i <= $shuffles; $i++){

	# need to reset our master set of sequences
	# we'll use tailswap mode to treat Chr4 differently
	%chr_seqs = FRAG::initalize_virtual_chromosome_sequences('tailswap');

	if ($i == 0){
		warn "Main run with unshuffled data\n";
	} else{
		warn "Shuffle $i\n";	
		# need to shuffle location of junctions
		shuffle_breakpoints();
	}

	# now mask $bp base pairs around each breakpoints (for sequences stored in %chr_seqs)
	%chr_seqs = FRAG::mask_breakpoint_regions($bp, \%chr_seqs, \%breakpoints);
  
 
	##########################################
	# Main loop over each possible feature
	##########################################

	# track all results by final ratio and by difference
	my %tmp_results;
	
	
	# now loop over all GFF feature sin our input file	
	foreach my $gff_feature (sort (FRAG::get_list_of_GFF_features($feature_gff))){
		warn "\tProcessing $gff_feature data\n" if ($verbose);
		
		# create copies of virtual chromosome sequences
		my %tmp_chr_seqs = %chr_seqs;
 
		# count how many times we see this feature
		my $feature_count = 0;
	 
		open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

		LINE: while(my $line = <$in>){
			chomp($line);
			# skip GFF header lines
			next if ($line =~ m/^#/);
			
			my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

			# only want to look at one feature at a time
			if ($gff_feature =~ m/open_chromatin_state/){
				my ($number) = $gff_feature =~ m/_(\d)/;

				next LINE unless ($feature =~ m/open_chromatin_state/ and $comment =~ m/state$number/);
			} else{
				next LINE unless $feature eq $gff_feature;
			}	
	
			# temporarily skip if not chr1 or chr4?
			next unless (($chr eq 'Chr1') or ($chr eq 'Chr4'));
	
			# skip tailswap regions of Chr1 and Chr4
			next if ($chr eq 'Chr1' and $s > $chr_sizes{'Chr1'});
			next if ($chr eq 'Chr4' and $s < $pre_tailswap_length);
	
			# mask where feature occurs in virtual chromosome sequence
			my $length = $e - $s + 1;
			substr($tmp_chr_seqs{$chr}, $s, $length) = ("o" x $length);
		
			$feature_count++;
		}
		close($in);
	
		# possible that we will get here without seeing any of the desired features
		# (depending on what input GFF is used), so skip if $feature_count == 0
		next if ($feature_count == 0);
	
	
		# now compare patterns in original virtual sequence (just breakpoints)
		# and tmp virtual sequence (masked with feature)
		my $original_breakpoint_bp      = $chr_seqs{'Chr1'} =~ tr/B/B/;
		$original_breakpoint_bp        += $chr_seqs{'Chr4'} =~ tr/B/B/;

		my $original_non_breakpoint_bp  = $chr_seqs{'Chr1'} =~ tr/-/-/;
		$original_non_breakpoint_bp    += $chr_seqs{'Chr4'} =~ tr/-/-/;

		my $remaining_breakpoint_bp     = $tmp_chr_seqs{'Chr1'} =~ tr/B/B/;
		$remaining_breakpoint_bp       += $tmp_chr_seqs{'Chr4'} =~ tr/B/B/;

		my $remaining_non_breakpoint_bp = $tmp_chr_seqs{'Chr1'} =~ tr/-/-/;
		$remaining_non_breakpoint_bp   += $tmp_chr_seqs{'Chr4'} =~ tr/-/-/;
	
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
		
		# store all of the current results in tmp hash
		$tmp_results{$gff_feature}{total_breakpoint_bp}       = $original_breakpoint_bp;
		$tmp_results{$gff_feature}{total_non_breakpoint_bp}   = $original_non_breakpoint_bp;
		$tmp_results{$gff_feature}{feature_breakpoint_bp}     = $feature_overlapping_breakpoint_regions;
		$tmp_results{$gff_feature}{feature_non_breakpoint_bp} = $feature_overlapping_non_breakpoint_regions;
		$tmp_results{$gff_feature}{feature_breakpoint_pc}     = $percent_overlapping_breakpoint_regions;
		$tmp_results{$gff_feature}{feature_non_breakpoint_pc} = $percent_overlapping_non_breakpoint_regions;
		$tmp_results{$gff_feature}{ratio}                     = $tmp_ratio;

		# do a couple of things differently now based on whether this is the
		# first run (main results) or a shuffled result

		# assess significance
		my $stats_ratio;
				
		if ($i == 0){
			$main_results{$gff_feature}{ratio_over}  = 0;
			$main_results{$gff_feature}{ratio_under} = 0;
			$main_results{$gff_feature}{ratio_same}  = 0;
		
			$main_results{$gff_feature}{ratio} = $tmp_ratio;
			$stats_ratio = 1;

		} else{		
			# if we are shuffling (i.e. $i > 0), we want to count whether the current ratio 
			# (from shuffling) beats real ratio in unshuffled data
			if($tmp_ratio > $main_results{$gff_feature}{ratio}){
				$main_results{$gff_feature}{ratio_over}++;
			} elsif($tmp_ratio < $main_results{$gff_feature}{ratio}){
				$main_results{$gff_feature}{ratio_under}++;
			} else{
				$main_results{$gff_feature}{ratio_same}++;
			}
			# do something different for stat testing depending on whether main 
			# ratio was above or below 1 (or zero)
			if($main_results{$gff_feature}{ratio} >= 1){
				$stats_ratio = $main_results{$gff_feature}{ratio_over} / $i;				
			} elsif($main_results{$gff_feature}{ratio} == 0){
				$stats_ratio = $main_results{$gff_feature}{ratio_same} / $i;							
			} else{
				$stats_ratio = $main_results{$gff_feature}{ratio_under} / $i;							
			}
		}

		print "$i\t";
		print "$main_results{$gff_feature}{ratio}\t";
		print "$bp\t";
		print "$gff_feature\t";
		print "$original_breakpoint_bp\t";
		print "$original_non_breakpoint_bp\t";
		print "$feature_overlapping_breakpoint_regions\t";
		print "$percent_overlapping_breakpoint_regions\t";
		print "$feature_overlapping_non_breakpoint_regions\t";
		print "$percent_overlapping_non_breakpoint_regions\t";
		print "$tmp_ratio\t";
		print "$main_results{$gff_feature}{ratio_over}\t";
		print "$main_results{$gff_feature}{ratio_same}\t";
		print "$main_results{$gff_feature}{ratio_under}\t";

		if ($main_results{$gff_feature}{ratio} == 0){
		}  elsif ($stats_ratio == 0){
			print "*****";
		} elsif ($stats_ratio < 0.0001){
			print "****";
		} elsif ($stats_ratio < 0.001){
			print "***";
		} elsif ($stats_ratio < 0.01){
			print "**";
		} elsif ($stats_ratio < 0.05){
			print "*";
		}

		print "\n";

	}
	
	
}


exit;


####################################
#
# S U B R O U T I N E S
#
####################################





###############################################################
# randomly shuffle junction coordinates
###############################################################

sub shuffle_breakpoints{

	foreach my $breakpoint (keys %breakpoints){
		
		# choose random position anywhere between length of sequence representing
		# chr1 (pre-tailswap) + region on chr4 after tailswap
		# have to factor in size of $bp as will need to extract a window either side 
		# of breakpoint
		
		my $max_coord = $chr_sizes{'Chr1'} + ($chr_sizes{'Chr4'} - $pre_tailswap_length + 1);
		
		my $random_coord = int(rand($max_coord));
		
		# is this coordinate on chromosome 1 but more 1/2 a window away from ends?
		if (($random_coord < ($chr_sizes{'Chr1'} - $bp/2)) and ($random_coord > ($bp/2))){
			$breakpoints{$breakpoint}{chr} = 'Chr1';
			$breakpoints{$breakpoint}{pos} = $random_coord;
		} elsif(($random_coord > ($chr_sizes{'Chr1'} + $bp/2)) and ($random_coord < ($max_coord - $bp/2))){
			# is this coordinate on chromosome 4 but more than 1/2 a window away from ends?
	
			$breakpoints{$breakpoint}{chr} = 'Chr4';
			$breakpoints{$breakpoint}{pos} = $random_coord - $chr_sizes{'Chr1'} + $pre_tailswap_length;		
			
		} else{
			# try again
			redo;
		}		
	}

}
