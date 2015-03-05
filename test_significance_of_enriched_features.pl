#!/usr/bin/perl
#
# test_significance_of_enriched_features
#
# A script to look at enrichment of genes and replication origins in 
# 2x and 3x breakpoint regions and test (statistically) whether
# this is significant
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use List::Util qw(shuffle sum);
use Getopt::Long;
use FRAG;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($breakpoint_gff, $feature_gff, $shuffles, $target_feature, $help, $verbose) = (undef, undef, 1000, undef, undef, undef);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <GFF file of gene and replication origin features>
--target_feature <gene or DNA_replication_origin>

Optional arguments:
--shuffles <how many shuffling iterations to perform (default = $shuffles)>
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"breakpoint_gff=s" => \$breakpoint_gff,
	"feature_gff=s"    => \$feature_gff,
	"target_feature=s" => \$target_feature,
	"shuffles=i"       => \$shuffles,
	"help"             => \$help,
	"verbose"          => \$verbose,
);


die $usage if ($help);
die $usage if (not defined $breakpoint_gff);
die $usage if (not defined $feature_gff);
die $usage if (not defined $target_feature);

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

# read features, 
my %features = FRAG::read_feature_data($feature_gff);

# hash that will just store two distances:
# distance to nearest feature before & after tailswap region
my %features_near_tailswap;




####################################
#
# Main loop
#
####################################

my $real_average_distance;
my @shuffled_averages;

# need to do one master run and then a number of shuffles to test significance (as denoted by $shuffles)
for (my $i = 0; $i <= $shuffles; $i++){

	print "\n### SHUFFLE $i ###\n";
	if ($i == 0){
		warn "Main run with unshuffled $target_feature data\n";		

	} else{
		warn "Shuffle $i\n";	
		# need to shuffle location of features
		shuffle_features();
	}

	# work out distance of last Chr1 $target_feature before tailswap and distance of first Chr4 $target_feature after tailswap
	find_features_near_tailswap();

	# can now process breakpoints
	process_breakpoints($i);
	
}

print "\n";
print "REAL: Average distance to nearest $target_feature = $real_average_distance bp\n";
my @sorted_averages = sort {$a <=> $b} @shuffled_averages;
print "SHUFFLED: @sorted_averages\n\n";

my $lower = 0;
foreach my $d (@shuffled_averages){
	$lower++ if ($d < $real_average_distance);
}

print "Shuffled datasets produced lower average distance than real average distance $lower times out of $shuffles\n\n";

exit;


####################################
#
#   S U B R O U T I N E S
#
####################################


sub process_breakpoints{

	my ($shuffle) = @_;
	
	# want to keep track of last feature before each breakpoint to
	# speed up subsequent searches (resume feature loop from this point)
	my $last_feature_before_breakpoint = "1";

	# will store all distances in array
	my @distances;

	BREAKPOINT: foreach my $breakpoint (sort {$a <=> $b} keys %breakpoints){
		
		my $b_pos = $breakpoints{$breakpoint}{pos};
		my $b_chr = $breakpoints{$breakpoint}{chr};
		
		my $distance_to_nearest_feature = 100000000;
		my $nearest_feature = "NA";
		my $nearest_chr = "NA";
		my $nearest_start = "NA";
		my $nearest_end = "NA";
		my $skip_count = 0;

		# find distance to nearest gene or replication origin

		FEATURE: foreach my $feat (sort {$a <=> $b} keys %features){ 
			my $f_chr = $features{$feat}{'chr'};
			my $s   = $features{$feat}{'start'};
			my $e   = $features{$feat}{'end'};

			# can skip to the last feature we looked at for previous breakpoint
			next FEATURE if ($feat < $last_feature_before_breakpoint);

			next FEATURE if ($f_chr ne $b_chr);


			# first check to see whether breakpoint overlaps with the current feature
			# if so we can just skip to the next breakpoint and declare a distance of 0 bp
			if (($f_chr eq $b_chr) and ($b_pos >= $s) and ($b_pos <= $e)){
				print "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\tDistance = 0\t*\n\n";
				push(@distances, 0);
				$last_feature_before_breakpoint = $feat - 1;
				next BREAKPOINT;
			} 

			# breakpoint doesn't overlap feature, so calculate distance
			my $distance;
			if ($f_chr eq $b_chr and ($b_pos < $s)){ # same chromosome, breakpoint is before feature
				$distance = $s - $b_pos;
			} elsif ($f_chr eq $b_chr and ($b_pos > $e)){ # same chromosome, breakpoint is after feature
				$distance = $b_pos - $e;
			} else{
				print "Uh oh 1: shouldn't get here! feat_chr = $f_chr $s-$e\n";

			}


			#print "BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\t$distance\t$distance_to_nearest_feature\tSkip = $skip_count\n" unless ($shuffle == 0);

			# is the current distance closer than previous best?
			if ($distance < $distance_to_nearest_feature){
				$distance_to_nearest_feature = $distance;
				$nearest_feature = $feat;
				$nearest_chr = $f_chr;
				$nearest_start = $s;
				$nearest_end = $e;
				$last_feature_before_breakpoint = $feat - 1;
				$skip_count = 0;
			} else{
				# must be going further away, but could have something strange like gene inside genes (especially with shuffling)
				# will allow looking at up to 20 more genes
				$skip_count++;
				last if ($skip_count > 20);
			}
		}
		print "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $nearest_feature: $nearest_chr $nearest_start-$nearest_end\tDistance = $distance_to_nearest_feature\t**\n";
		push(@distances, $distance_to_nearest_feature);
		# could we be looking at a breakpoint where the nearest feature is actually after the tailswap
		# region? I.e. on the other chromosome of the reference genome?

		if($nearest_chr eq 'Chr1'){
			my $distance_to_end_of_chr = $chr_sizes{'Chr1'} - $nearest_end;
			my $distance_to_post_tailswap_feature = $distance_to_end_of_chr + $features_near_tailswap{'after'};
			print "\tDistance to nearest feature on same chromosome = $distance_to_nearest_feature\n";
			print "\tDistance to nearest feature after tailswap = $distance_to_post_tailswap_feature\n\n";	
			die "Uh oh 2!\n" if ($distance_to_post_tailswap_feature < $distance_to_nearest_feature);
		} else{
			my $distance_to_pre_tailswap_feature = $nearest_start - $pre_tailswap_length + $features_near_tailswap{'before'};
			print "\tDistance to nearest feature on same chromosome = $distance_to_nearest_feature\n";
			print "\tDistance to nearest feature before tailswap = $distance_to_pre_tailswap_feature\n\n";	
			die "Uh oh 3!\n" if ($distance_to_pre_tailswap_feature < $distance_to_nearest_feature);

		}
	}

	print "\n\nDISTANCES: @distances\n";
	my $n = @distances;
	my $average_distance = sprintf("%.1f", sum(@distances) / $n);

	print "Average distance to nearest $target_feature (from $n breakpoints) = $average_distance bp\n\n";
	
	# store results
	if ($shuffle == 0){
		$real_average_distance = $average_distance;	
	} else{
		push(@shuffled_averages, $average_distance);
	}
}

sub find_features_near_tailswap{

	# reset hash details with some default values
	$features_near_tailswap{'before'} = 1000000;
	$features_near_tailswap{'after'}  = 1000000;

	foreach my $feat (sort {$a <=> $b} keys %features){ 

		my $chr = $features{$feat}{'chr'};
		my $s   = $features{$feat}{'start'};
		my $e   = $features{$feat}{'end'};

#		print "$feat $chr $s-$e\t$chr_sizes{'Chr1'}\t$pre_tailswap_length\n";

		# want to see if this feature is closer to the tailswap region than
		# any previous feature, if so store the distance. Do this for 
		# features on chr1 *before* the tailswap, and features on chr4 *after* the tailswap
		if ($chr eq 'Chr1' and $e < $chr_sizes{'Chr1'}){
			my $distance = $chr_sizes{'Chr1'} - $e;
			($features_near_tailswap{'before'} = $distance) if ($distance < $features_near_tailswap{'before'});
		}
		
		if ($chr eq 'Chr4' and $s > $pre_tailswap_length){
			my $distance = $s - $pre_tailswap_length;
			($features_near_tailswap{'after'} = $distance) if ($distance < $features_near_tailswap{'after'});
		}

	}

	print "\n";
	print "Distance of nearest $target_feature *after* tailswap on Chr4: $features_near_tailswap{'after'}\n";
	print "Distance of nearest $target_feature *before* tailswap on Chr1: $features_near_tailswap{'before'}\n";
	print "\n";
	

}

###############################################################
# randomly shuffle feature coordinates
###############################################################

sub shuffle_features{

	# temporarily put shuffled features into two hashes (chr1 and chr4)
	# and then combine afterwards (so we can have them sorted by position)
	my %chr1_feats;
	my %chr4_feats;
	
	foreach my $feat (sort {$a <=> $b} keys %features){
		
		
		#print "$feat) Feature = $features{$feat}{'chr'} $features{$feat}{'start'}-$features{$feat}{'end'}\n";
		my $feat_length = $features{$feat}{'end'} - $features{$feat}{'start'} + 1;
	
		# choose random position anywhere between length of sequence representing
		# chr1 (pre-tailswap) + region on chr4 after tailswap
		# have to factor in size of $bp as will need to extract a window either side 
		# of breakpoint
		my $max_coord = $chr_sizes{'Chr1'} + ($chr_sizes{'Chr4'} - $pre_tailswap_length + 1);
		my $random_coord = int(rand($max_coord));
		
		# is this coordinate on chromosome 1 and less than feature length away from tailswap boundary?
		if ($random_coord < ($chr_sizes{'Chr1'} - $feat_length)){
			my $s = $random_coord;
			my $e = $s + $feat_length - 1;
			
			# add to temp hash as long as we haven't seen this start coordinate before
			if(exists($chr1_feats{$s})){
				redo;
			} else{
				$chr1_feats{$s} = $e;	
			}
			
			
		} elsif(($random_coord > $chr_sizes{'Chr1'}) and ($random_coord < ($max_coord - $feat_length))){
			# coordinate must be on chromosome 4
			my $s = $random_coord - $chr_sizes{'Chr1'} + $pre_tailswap_length;
			my $e = $s + $feat_length - 1;
			
			# add to temp hash as long as we haven't seen this start coordinate before
			if(exists($chr4_feats{$s})){
				redo;
			} else{
				$chr4_feats{$s} = $e;	
			}
		} else{
			# try again
			redo;
		}	
		#print "$feat) Feature = $features{$feat}{'chr'} $features{$feat}{'start'}-$features{$feat}{'end'}\n\n";
	}

	# now empty %features hash, sort %chr1_feats and %chr4_feats and add back to %features
	%features = ();

	my $feat_counter = 0;
	foreach my $start (sort {$a <=> $b} keys %chr1_feats){
		$feat_counter++;
		$features{$feat_counter}{'chr'}   = "Chr1";
		$features{$feat_counter}{'start'} = $start;
		$features{$feat_counter}{'end'}   = $chr1_feats{$start};
	}
	foreach my $start (sort {$a <=> $b} keys %chr4_feats){
		$feat_counter++;
		$features{$feat_counter}{'chr'}   = "Chr4";
		$features{$feat_counter}{'start'} = $start;
		$features{$feat_counter}{'end'}   = $chr4_feats{$start};
	}
	
	
}






