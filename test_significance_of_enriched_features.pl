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

my ($breakpoint_gff, $feature_gff, $shuffles, $help, $verbose) = (undef, undef, 0, undef, undef, undef);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <GFF file that contains a *single* GFF feature>

Optional arguments:
--shuffles <how many shuffling iterations to perform (default = $shuffles)>
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"breakpoint_gff=s" => \$breakpoint_gff,
	"feature_gff=s"    => \$feature_gff,
	"shuffles=i"       => \$shuffles,
	"help"             => \$help,
	"verbose"          => \$verbose,
);


die $usage if ($help);
die $usage if (not defined $breakpoint_gff);
die $usage if (not defined $feature_gff);

# determine copy number
my $copy_number;
$copy_number = "2x" if ($breakpoint_gff =~ m/2x/);
$copy_number = "3x" if ($breakpoint_gff =~ m/3x/);


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

# what are the features…should only be one
my @features = FRAG::get_list_of_GFF_features($feature_gff);
die "ERROR: More than 1 feature in $feature_gff\n\n$usage" unless (@features == 1);
my $target_feature = $features[0];



# hashe that will just store distances and IDs:
# of nearest features and breakpoints before & after tailswap region
my %tailswap_details;

####################################
#
# Main loop
#
####################################

my $real_average_distance;
my $real_standard_deviation;
my @real_distances;
my @shuffled_averages;

my $exception_counter1 = 0;
my $exception_counter2 = 0;


# work out distance of last Chr1 $target_feature before tailswap and distance of first Chr4 $target_feature after tailswap
find_features_near_tailswap();

print "Run\tCopy_number\tFeature\tAverage_distance_to_nearest_feature\tStandard_deviation\t#_of_times_shuffled_average_was_below_real_average\n";

# need to do one master run and then a number of shuffles to test significance (as denoted by $shuffles)
for (my $i = 0; $i <= $shuffles; $i++){
	print "$i\t$copy_number\t$target_feature\t";
	
	# need to shuffle location of breakpoints
	shuffle_breakpoints() if ($i != 0);
		
	# we always need to work out whether breakpoints near the tailswap region might actually be closer
	# to a feature on the other side of the tailswap. More likely to happen when you shuffle breakpoints.
	find_breakpoints_near_tailswap();

	# can now process breakpoints
	process_breakpoints($i);
	
	# can write the real results to a file
	next unless ($i == 0);
	my $outfile = "real_distances_${copy_number}_$target_feature.tsv";
	open(my $out, ">", "$outfile") or die "can't write to $outfile\n";
	foreach my $d (sort {$a <=> $b} @real_distances){
		print $out "$d\n";
	}
	close($out);
}


# can now write the file of average distances from all of the shuffle runs

my @sorted_shuffled_averages = sort {$a <=> $b} @shuffled_averages;
my $outfile = "shuffled_distances_${copy_number}_$target_feature.tsv";

open(my $out, ">", "$outfile") or die "can't write to $outfile\n";
my $lower = 0;
foreach my $d (@sorted_shuffled_averages){
	$lower++ if ($d < $real_average_distance);
	print $out "$d\n";
}

close($out);

# summary
warn "\n";
warn "REAL: Average distance to nearest $target_feature = $real_average_distance bp (std dev = $real_standard_deviation)\n";
warn "Shuffled datasets produced lower average distance than real average distance $lower times out of $shuffles\n\n";
warn "There were $exception_counter1 instances of nearest feature being *after* tailswap\n"  if ($verbose);
warn "There were $exception_counter2 instances of nearest feature being *before* tailswap\n" if ($verbose);
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
			my $s     = $features{$feat}{'start'};
			my $e     = $features{$feat}{'end'};

			# apart from the possible issue of the features next to tailswap region,
			# don't need to look at current feature if on a different chromosome
			next FEATURE unless ($f_chr eq $b_chr);

			# can skip to the last feature we looked at for previous breakpoint
			next FEATURE if ($feat < $last_feature_before_breakpoint);

			
			# first check to see whether breakpoint overlaps with the current feature
			# if so we can just skip to the next breakpoint and declare a distance of 0 bp
			if (($b_pos >= $s) and ($b_pos <= $e)){
				warn "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\tDistance = 0\t*\n" if ($verbose);
				push(@distances, 0);
				$last_feature_before_breakpoint = $feat - 1;
				next BREAKPOINT;
			} 


			# Now we know that current breakpoint doesn't overlap current feature, so calculate distance
			my $distance;

			# Is the current feature before or after the breakpoint
			if ($b_pos < $s){ 
				$distance = $s - $b_pos;
			} else{
				$distance = $b_pos - $e;
			}

			#warn "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\tD=$distance\tDTN=$distance_to_nearest_feature\tSkip = $skip_count\n";

			# is the current distance closer than previous best?
			if ($distance < $distance_to_nearest_feature){
				$distance_to_nearest_feature = $distance;
				
				$nearest_feature = $feat;
				$nearest_chr = $f_chr;
				$nearest_start = $s;
				$nearest_end = $e;
				$last_feature_before_breakpoint = $feat - 1;
				$skip_count = 0;
			} elsif((($s >= $nearest_start) and ($s <= $nearest_end)) or
					(($e >= $nearest_start) and ($e <= $nearest_end))){
				# check for overlapping features (e.g. gene inside gene, or nested repeats)
				# if features overlap, then we don't really mind that the the current distance
				# exceeds the nearest distance. Just move along to next feature on chromosome
				next FEATURE;
			} else{
				# must be going further away, but maybe annotations are not in coordinate order for some reason?
				# will allow looking at up to 5 more features to be sure that the $distance_to_nearest_feature
				# is the best distance possible
				$skip_count++;
				last if ($skip_count > 5);
			}


			# Next check whether this is the last breakpoint before the tailswap		
			# *and* the last feature before the tailswap
			if ($feat       eq $tailswap_details{'feature_before'}{'id'} and
				$breakpoint eq $tailswap_details{'breakpoint_before'}{'id'}){
				
				# set the nearest possible feature to be the first feature following the tailswap
				# but only if we haven't already seen a feature that's closer
				warn "\n*BEFORE TAILSWAP*\n" if ($verbose); 
				my $distance_to_post_tailswap_feature = $chr_sizes{'Chr1'} - $b_pos + $tailswap_details{'feature_after'}{'distance'};
				warn "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\t" if ($verbose);
				warn "Distance = $distance_to_nearest_feature\tDistance to post-tailswap feature = $distance_to_post_tailswap_feature\n" if ($verbose);		
				if ($distance_to_post_tailswap_feature < $distance_to_nearest_feature){
					
					$distance_to_nearest_feature = $distance_to_post_tailswap_feature;
					$nearest_feature = $tailswap_details{'feature_after'}{'id'};
					$nearest_start   = $tailswap_details{'feature_after'}{'s'};
					$nearest_end     = $tailswap_details{'feature_after'}{'e'};
					$nearest_chr     = 'Chr4';
					$skip_count = 0;
					warn "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\tDistance = $distance_to_nearest_feature\n" if ($verbose);
					$exception_counter1++;
				}
				warn "\n" if ($verbose);
			}


			# Now check whether this is the first breakpoint after the tailswap 
			# and the first feature after the tailswap
			if ($feat       eq $tailswap_details{'feature_after'}{'id'} and
				$breakpoint eq $tailswap_details{'breakpoint_after'}{'id'}){
				warn "\n*AFTER TAILSWAP*\n"  if ($verbose);
				my $distance_to_pre_tailswap_feature = $b_pos - $pre_tailswap_length + $tailswap_details{'feature_before'}{'distance'};
				warn "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\t"  if ($verbose);
				warn "Distance = $distance_to_nearest_feature\tDistance to pre-tailswap feature = $distance_to_pre_tailswap_feature\n"  if ($verbose);
		
				if ($distance_to_pre_tailswap_feature < $distance_to_nearest_feature){

					# set the nearest possible feature to be the last feature preceding the tailswap
					$distance_to_nearest_feature = $distance_to_pre_tailswap_feature;
					$nearest_feature = $tailswap_details{'feature_before'}{'id'};				
					$nearest_start = $tailswap_details{'feature_before'}{'s'};
					$nearest_end = $tailswap_details{'feature_before'}{'e'};
					$nearest_chr = 'Chr1';
					$last_feature_before_breakpoint = $feat - 1;
					$skip_count = 0;
					warn "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $feat: $f_chr $s-$e\tDistance = $distance_to_nearest_feature\n"  if ($verbose);
					
					$exception_counter2++;
					#exit unless ($shuffle == 0);
				}
				warn "\n" if ($verbose);
			}
		}

		warn "SHUFFLE $shuffle: BREAKPOINT $breakpoint: $b_chr $b_pos\t$target_feature $nearest_feature: $nearest_chr $nearest_start-$nearest_end\tDistance = $distance_to_nearest_feature\t**\n" if ($verbose);
		push(@distances, $distance_to_nearest_feature);
	}

	# calculate stats
	my $n                = @distances;
	my $average_distance = sprintf("%.1f", sum(@distances) / $n);
	my $stdev            = sprintf("%.1f", sqrt(sum(map {($_ - $average_distance) ** 2} @distances) / ($n-1)));
	
	#print "Average distance to nearest $target_feature (from $n breakpoints) = $average_distance bp\n\n";
	print "$average_distance\t$stdev\t";

	# store results
	if ($shuffle == 0){
		$real_average_distance = $average_distance;	
		$real_standard_deviation = $stdev;
		@real_distances = @distances;
		print ".";
	} else{
		push(@shuffled_averages, $average_distance);

		# how many of the shuffled averages are below the real distance
		my $lower = 0;
		foreach my $d (@shuffled_averages){
			$lower++ if ($d < $real_average_distance);
		}
		print "$lower";
	}
	
	print "\n";

}

sub find_features_near_tailswap{

	# reset hash details with some default values
	$tailswap_details{'feature_before'}{'id'} = "";
	$tailswap_details{'feature_before'}{'distance'} = 10000000;
	$tailswap_details{'feature_after'}{'id'} = "";
	$tailswap_details{'feature_after'}{'distance'} = 10000000;
	

	foreach my $feat (sort {$a <=> $b} keys %features){ 

		my $chr = $features{$feat}{'chr'};
		my $s   = $features{$feat}{'start'};
		my $e   = $features{$feat}{'end'};

		# want to see if this feature is closer to the tailswap region than
		# any previous feature, if so store the distance and gene ID. Do this for 
		# features on chr1 *before* the tailswap, and features on chr4 *after* the tailswap
		if ($chr eq 'Chr1' and $s < $chr_sizes{'Chr1'}){		
			my $distance = $chr_sizes{'Chr1'} - $e;
			
			# special case of feature that might span tailswap region
			$distance = 0 if ($e > $chr_sizes{'Chr1'});
			
			if ($distance < $tailswap_details{'feature_before'}{'distance'}){
				$tailswap_details{'feature_before'}{'distance'} = $distance;
				$tailswap_details{'feature_before'}{'id'} = $feat;
				$tailswap_details{'feature_before'}{'s'} = $s;
				$tailswap_details{'feature_before'}{'e'} = $e;
			}
		}
		
		if ($chr eq 'Chr4' and $e > $pre_tailswap_length){
			my $distance = $s - $pre_tailswap_length;
			
			# special case of feature that might span tailswap region
			$distance = 0 if ($s < $pre_tailswap_length);
			if ($distance < $tailswap_details{'feature_after'}{'distance'}){
				
				$tailswap_details{'feature_after'}{'distance'} = $distance;
				$tailswap_details{'feature_after'}{'id'} = $feat;
				$tailswap_details{'feature_after'}{'s'} = $s;
				$tailswap_details{'feature_after'}{'e'} = $e;
			}
		}

	}

	if ($verbose){
		warn "Nearest $target_feature *before* tailswap on Chr1: ID = $tailswap_details{'feature_before'}{'id'}, ";
		warn "$tailswap_details{'feature_before'}{'s'}-$tailswap_details{'feature_before'}{'e'}, " ;
		warn "Distance = $tailswap_details{'feature_before'}{'distance'} bp\n";


		warn "Nearest $target_feature *after* tailswap on Chr4:  ID = $tailswap_details{'feature_after'}{'id'}, ";
		warn "$tailswap_details{'feature_after'}{'s'}-$tailswap_details{'feature_after'}{'e'}, ";
		warn "Distance = $tailswap_details{'feature_after'}{'distance'} bp\n";
		warn "\n";
	}	
}


sub find_breakpoints_near_tailswap{

	# when we shuffle data, we want to keep track of what is the last breakpoint
	# before the tailswap, and the first breakpoint after the tailswap

	# reset hash details with some default values
	$tailswap_details{'breakpoint_before'}{'id'} = "";
	$tailswap_details{'breakpoint_before'}{'distance'} = 100000000;
	$tailswap_details{'breakpoint_after'}{'id'} = "";
	$tailswap_details{'breakpoint_after'}{'distance'} = 100000000;
	

	foreach my $bp (sort {$a <=> $b} keys %breakpoints){ 

		my $chr = $breakpoints{$bp}{'chr'};
		my $pos = $breakpoints{$bp}{'pos'};
		

		# want to see if this feature is closer to the tailswap region than
		# any previous feature, if so store the distance and gene ID. Do this for 
		# features on chr1 *before* the tailswap, and features on chr4 *after* the tailswap
		if ($chr eq 'Chr1' and $pos <= $chr_sizes{'Chr1'}){
			my $distance = $chr_sizes{'Chr1'} - $pos;
			
			if ($distance < $tailswap_details{'breakpoint_before'}{'distance'}){
				$tailswap_details{'breakpoint_before'}{'distance'} = $distance;
				$tailswap_details{'breakpoint_before'}{'id'} = $bp;
			}
		}
		
		if ($chr eq 'Chr4' and $pos >= $pre_tailswap_length){
			my $distance = $pos - $pre_tailswap_length;

			if ($distance < $tailswap_details{'breakpoint_after'}{'distance'}){
				$tailswap_details{'breakpoint_after'}{'distance'} = $distance;
				$tailswap_details{'breakpoint_after'}{'id'} = $bp;
			}
		}
	}
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


sub shuffle_breakpoints{

	# will use temp arrays to store shuffled breakpoint coordinates before reassigning to hash
	my @tmp_chr1;
	my @tmp_chr4;

	foreach my $breakpoint (keys %breakpoints){
		
		# choose random position anywhere between length of sequence representing
		# chr1 (pre-tailswap) + region on chr4 after tailswap
		# have to factor in size of $bp as will need to extract a window either side 
		# of breakpoint
		
		my $max_coord = $chr_sizes{'Chr1'} + ($chr_sizes{'Chr4'} - $pre_tailswap_length + 1);
		my $random_coord = int(rand($max_coord));
		
		# add feature to chr1 or chr4 temp hash
		if ($random_coord < $chr_sizes{'Chr1'}){
			push(@tmp_chr1, $random_coord);		
		} else{
			# need to change any breakpoint on chr4 to use chr coordinate ranges
			my $adjusted_coord = $random_coord - $chr_sizes{'Chr1'}  + $pre_tailswap_length;
			push(@tmp_chr4, $adjusted_coord);
		}	
	}

	my $breakpoint_counter = 0;
	foreach my $pos (sort {$a <=> $b} @tmp_chr1){
		$breakpoint_counter++;
		$breakpoints{$breakpoint_counter}{chr} = 'Chr1';
		$breakpoints{$breakpoint_counter}{pos} = $pos;
	}
	foreach my $pos (sort {$a <=> $b} @tmp_chr4){
		$breakpoint_counter++;
		$breakpoints{$breakpoint_counter}{chr} = 'Chr4';
		$breakpoints{$breakpoint_counter}{pos} = $pos;
	}
}
