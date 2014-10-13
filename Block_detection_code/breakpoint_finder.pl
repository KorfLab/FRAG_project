#!/usr/bin/perl
#
# breakpoint_finder.pl
#
# A script to find 'peaks' of duplicated/triplicated regions in FRAG project data
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

use strict;
use warnings FATAL => 'all';
use List::Util qw(sum shuffle);
use Getopt::Long;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($file, $win_size, $shuffle_file) = (undef, 20, undef);
my ($min_t, $skip)  = (2.576, undef);
my ($min_mean_diff, $max_var, $help) = (0.65, 0.1, undef);

my $usage = "
$0 --file <TSV file of read count data> --win_size <default = $win_size> 

optional arguments:
		--shuffle_file : use a shuffle file
		--min_t : minimum t-test value to consider significant <default = $min_t>
		--min_mean_diff : minimum difference expected between average copy number within and outside window <default = $min_mean_diff
		--max_var : maximum variance allowed within a window (scaled upwards for small windows) <default = $max_var>
		--skip : skip rest of chromosome when the start chunk coordinate is greater than that specified
        --help
        
";

GetOptions (
	"file=s"            => \$file,
	"win_size=i"        => \$win_size,
	"shuffle_file=s"    => \$shuffle_file,
	"min_t=f"           => \$min_t,
	"min_mean_diff=f"   => \$min_mean_diff,
	"max_var=f"         => \$max_var,
	"skip=i"            => \$skip,
	"help"              => \$help,
);

die $usage if ($help);
die $usage if (not defined $file);

warn "Parameters:
--win_size = $win_size
--min_t = $min_t
--min_mean_diff = $min_mean_diff
--max_var = $max_var\n\n";

####################################
# Global variables
####################################

# keep track of maximum value of t observed at each window size
my %max_t;

# need to store count data in hash of arrays (primary key is chromosome name)
my %read_counts;

# will want to store all significant scoring windows found in each chunk
my %sig_windows;


# read shuffled data from file (if specified)
read_shuffle_data() if ($shuffle_file);

# grab read count data from file
read_counts_file();

# main loop over each chromosome
warn "Looping over chromosome in windows of size $win_size\n";
print "Chr\tWindow\tWin:start-end\tType\tCopy_number\tStart_coord\tEnd_coord\tt\tMean\tVariance\n";
	
# 1st loop, just loop over desired chromosomes worth of data
foreach my $chr (qw(Chr1 Chr4)){
	
	# reset significant windows
	%sig_windows = ();

	my $chr_size = @{$read_counts{$chr}};

	# now process all windows in this chunk
	loop_over_chromosome($chr, 0, $chr_size);

	exit; #DEBUG
	
	# need to keep track of previous significant details in order to generate blocks that span the two edges.
	my @previous_windows;
	my $counter = 0;

	# final output for current chunk, in ascending chromosome order
	WINDOW: foreach my $window (sort {$sig_windows{$a}{min} <=> $sig_windows{$b}{min}} keys %sig_windows){
		$counter++;

		my $t            = sprintf("%.2f", $sig_windows{$window}{t});
		my $mean1        = sprintf("%.2f", $sig_windows{$window}{mean1});
		my $mean2        = sprintf("%.2f", $sig_windows{$window}{mean2});
		my $var1         = $sig_windows{$window}{var1};
		my $var2         = $sig_windows{$window}{var2};
		my $midpoint_A   = $sig_windows{$window}{midA};
		my $midpoint_B   = $sig_windows{$window}{midB};
		my $copy_number1 = $sig_windows{$window}{copy_number1};
		my $copy_number2 = $sig_windows{$window}{copy_number2};

		my ($start, $end)   = split(/-/, $window);

		# now format output for printing
		my $start_pos      = sprintf("%.0f",($start * 1000));
		my $end_pos        = sprintf("%.0f",($end * 1000));
		my $midpoint_A_pos = sprintf("%.0f",($midpoint_A * 1000));
		my $midpoint_B_pos = sprintf("%.0f",($midpoint_B * 1000));

		
		# print first half of window
		print "$chr\t";
		print "$window\t";
		print "$start-$midpoint_A\t";
		print "edge-L\t";
		print "$copy_number1\t";
		print "$start_pos\t$midpoint_A_pos\t";
		print "$t\t";
		print "$mean1\t";
		print "$var1\n";

		# print 2nd half of window
		print "$chr\t";
		print "$window\t";
		print "$midpoint_B-$end\t";
		print "edge-R\t";
		print "$copy_number2\t";
		print "$midpoint_B_pos\t$end_pos\t";
		print "$t\t";
		print "$mean2\t";
		print "$var2\n";

		# create a block for first part of chromosome (might be erroneous assumption?)
		if($counter == 1){
			print "$chr\t";
			print "NA\t";
			print "0-$midpoint_A\t";
			print "block\t";
			print "$copy_number1\t";
			print "0\t$midpoint_A_pos\t";
			print "NA\t";
			print "NA\t";
			print "NA\n";

			# save details of current window and move to next window
			my $hash_ref = {'midA'         => $midpoint_A, 
							'midB'         => $midpoint_B,
							'copy_number1' => $copy_number1, 
							'copy_number2' => $copy_number2, 
							'midB_pos'     => $midpoint_B_pos};
							
			push(@previous_windows, $hash_ref);
			next WINDOW;			
		}				

		# if we are here we are looking at the 2nd or greater window
		# this means can print details of what should be the full block
		# between the end of the previous window (technically mid_point_A) and the start of current window 
		
		# extract details from previous window that we examined
		# will be last element in @previous_windows
		my $previous_midpoint_B      = ${$previous_windows[-1]}{midB};
		my $previous_midpoint_B_pos  = ${$previous_windows[-1]}{midB_pos};
		my $previous_copy_number     = ${$previous_windows[-1]}{copy_number2};
		
		# we can calculate what the variance will be for this block
		# if too high, we want to discount it
		my @block = @{$read_counts{$chr}}[$previous_midpoint_B..$midpoint_A];

		my $n = scalar @block;
		my $sum = sum(@block);		
		my $mean = $sum/$n;
		my $var = sum(map {($_ - $mean) ** 2} @block) / ($n - 1);

		# scale variance based on length of block
		my $scaled_max_var = $max_var * (30/(($n / 2) + 10) + 0.8);			

		print "# Block $previous_midpoint_B-$midpoint_A (L=$n)\tVariance = $var\tScaled max variance = $scaled_max_var";
		if ($var > $max_var and $var < $scaled_max_var){
			print "\t***\n";			
		} else {
			print "\n";
		}

		# only print block if variance is low enough
		if ($var > $scaled_max_var){
			print "# removing block as variance is too high ($var)\n";
			$counter++;
			# need to remove the information from @previous_windows and then skip 
			splice(@previous_windows, -1, 1);
			next WINDOW;				
		}
		
		print "$chr\t";
		print "NA\t";
		print "$previous_midpoint_B-$midpoint_A\t";
		print "block\t";
		print "$previous_copy_number\t";
		print "$previous_midpoint_B_pos\t$midpoint_A_pos\t";
		print "NA\t";
		print "NA\t";
		print "NA\n";
		
		# save details of current window and move to next window
		my $hash_ref = {'midA'         => $midpoint_A, 
						'midB'         => $midpoint_B,
						'copy_number1' => $copy_number1, 
						'copy_number2' => $copy_number2, 
						'midB_pos'     => $midpoint_B_pos};							
		push(@previous_windows, $hash_ref);

	}
}

exit;


####################################
#
# Subroutines
#
####################################

sub loop_over_chromosome{

	my ($chr, $start, $end) = @_;

	# track copy number of last 'half-edge' that we saw. E.g. if an edge spanned a 2x-3x junction,
	# the last half-edge was 3x, also want to store coordinates, so use a hash
	my %last_window;

	# generate initial array of read counts for first window that we inspect
	my @values_A;
	
	# Will compare @values_A to all *other* counts in current chunk
	my @values_B;
		
	# now loop along windows in current chunk (data from @{$read_counts{}} array)
	WINDOW: for(my $i = $start; $i <= $end - $win_size; $i++){

		my ($min, $max) = ($i, $i + $win_size - 1);
		
		# do we want to use our shortcut to end early/
		if($skip and ($min > $skip)){
			last WINDOW;
			#die "window start $min exceeded by --skip $skip option, exiting early\n";
		}

		my $mid_point_A = $min + ($win_size / 2) - 1;		
		my $mid_point_B = $mid_point_A + 1;

		# do things a bit differently for the first window
		if($i == $start){
			# set default contents of @values_A, don't need to do anything else
			@values_A = @{$read_counts{$chr}}[$min..$mid_point_A];
			@values_B = @{$read_counts{$chr}}[$mid_point_B..$max]; 
		} else{
			# get rid of one element from LHS of @values_A, and add a new element on RHS
			my $tmp = shift(@values_A);
			push(@values_A, shift(@values_B));	
						
			# now add one value to @values_B
			push(@values_B, ${$read_counts{$chr}}[$max]);	
		}
		
		# calculate t value
		my ($t, $mean1, $mean2, $var1, $var2) = t_test(\@values_A, \@values_B);			
		 
		# only want to proceed if t-test results are statistically significant (P < 0.001 t = ~3.29)			
		# if we have used a shuffled results file, we can now also ask whether
		# t-value could be expected by chance (ignore current t-value if this is the case)
		next unless ($t >= $min_t);		
		# next WINDOW if ($shuffle_file and $t <= $max_t{"$chr:$start-$end:$win_size"});
		 
		 
#		print "$chr\t$min-$max\t$min-$mid_point_A\t$mid_point_B-$max\tt=$t\tmeans = $mean1 vs $mean2\tvars = $var1 vs $var2\n";
#		warn "# A\t$chr\tCHR\t$min-$max\t$min-$mid_point_A\t$mid_point_B-$max\tblock\t2x\t${min}000\t${max}000\t$t\t$mean1\t$mean2\tmean_diff = $mean_diff\tvars = $var1 vs $var2\n";
		
		# Skip windows if the differences in the two means is too low. We are looking for 
		# cases where a window goes from 1x reads to 2-3x, so difference in means should 
		# really be at least 1.0, but we will be cautious and use a user-defined cutoff 
		# ($min_mean_diff) defaulting to 0.65
		my $mean_diff = sprintf("%.3f", abs($mean1 - $mean2));
		next unless ($mean_diff > $min_mean_diff);

		# scale max variance allowed for block (dependent on window size)
		my $scaled_max_var = $max_var * (30/(($win_size / 2) + 10) + 0.8);			

		# is there too much variance inside this window?
		# e.g. when a block incorrectly spans a 2x block and some of a 3x block, variance
		# will be higher than a smaller block that captures just the 2x region
		next unless ($var1 <= $scaled_max_var and $var2 <= $scaled_max_var);

		# now check if this current window overlaps a pre-existing window
		# which has a lower t-value score. If so, can remove the pre-existing window
		next unless check_for_redundancy($min, $max, $mean1, $t);			
		
		print  "$chr\t$min-$max\t$min-$mid_point_A\t$mid_point_B-$max\tt=$t\tmeans = $mean1 vs $mean2\tmean_diff = $mean_diff\tvars = $var1 vs $var2\tMax_var = $max_var, scaled_max_var = $scaled_max_var\n";
#		exit;

					
		# now calculate whether each half is in a 1x, 2x, or 3x block
		# this is a guesstimate
		my ($copy_number1, $copy_number2);
		
		if ($mean1 < 2){
			$copy_number1 = "1x";
		} elsif ($mean1 < 3){
			$copy_number1 = "2x";
		} else{
			$copy_number1 = "3x";
		}

		if ($mean2 < 2){
			$copy_number2 = "1x";
		} elsif ($mean2 < 3){
			$copy_number2 = "2x";
		} else{
			$copy_number2 = "3x";
		}

		
		# Before we store details of this new 'edge', we should first try to check
		# that we didn't miss an edge earlier that would leave to an errant block detection.
		# E.g. if we previously called a 2x block, and now we have an edge that goes from 3x to 2x
		# then there should have been a *second* 2x-3x edge which has been missed.
		# but don't check this if current window overlaps previous window (we have not yet found the best edge)
		if($last_window{copy_number}){
		
#			warn "MISSING BLOCK: Last half window = $last_window{start}-$last_window{end} $last_window{copy_number}, This half window = $min-$mid_point_A $copy_number1
# C\t$chr\tCHR\t$min-$max\t$min-$mid_point_A\t$mid_point_B-$max\tblock\t2x\t${min}000\t${max}000\t$t\t$mean1\t$mean2\tmean_diff = $mean_diff\tvars = $var1 vs $var2\n";
		
			if(($last_window{copy_number} ne $copy_number1) and	check_for_overlap($last_window{start}, $last_window{end}, $min, $mid_point_A) == 0){
				next WINDOW;
			}	
		}
		
#		warn "# B\t$chr\tCHR\t$min-$max\t$min-$mid_point_A\t$mid_point_B-$max\tblock\t2x\t${min}000\t${max}000\t$t\t$mean1\t$mean2\tmean_diff = $mean_diff\tvars = $var1 vs $var2\n";
#		warn "\n";

		
		# keep track of details of this significant window
		$sig_windows{"$min-$max"}{t}     = $t;
		$sig_windows{"$min-$max"}{min}   = "$min";				
		$sig_windows{"$min-$max"}{max}   = "$max";				
		$sig_windows{"$min-$max"}{mean1} = "$mean1";
		$sig_windows{"$min-$max"}{mean2} = "$mean2";				
		$sig_windows{"$min-$max"}{var1}  = "$var1";
		$sig_windows{"$min-$max"}{var2}  = "$var2";				
		$sig_windows{"$min-$max"}{midA}  = "$mid_point_A";				
		$sig_windows{"$min-$max"}{midB}  = "$mid_point_B";			
		$sig_windows{"$min-$max"}{copy_number1}  = "$copy_number1";
		$sig_windows{"$min-$max"}{copy_number2}  = "$copy_number2";
		
		# store a few details of the last significant window
		$last_window{start} = $mid_point_B;
		$last_window{end} = $max;
		$last_window{copy_number} = $copy_number2;
						
	}
}

# take current significant window and compare to other significant windows we have already 
# seen. If the current window overlaps with another window *and* has a higher t-score, 
# then we will keep it and remove the other window

sub check_for_redundancy{

	my ($start1, $end1, $mean1, $t1) = @_;

	foreach my $window (sort {$a cmp $b} keys %sig_windows){
			my ($start2, $end2) = split(/-/, $window);
			my $t2 = $sig_windows{$window}{t};
			my $mean2 = $sig_windows{$window}{mean1};
			my $size = keys(%sig_windows);
			
			# first ask whether windows overlap (in any way) and skip to next 
			# comparison if windows don't overlap
			next unless check_for_overlap($start1, $end1, $start2, $end2);

			my ($s1, $s2) = ($end1-$start1+1, $end2-$start2+1);
			# print "\tComparing $start1-$end1 (L=$s1, t=$t1) with $start2-$end2 (L=$s2, t=$t2)\n";
			
			# If two significant windows overlap, they might have opposing
			# means (i.e. one is significantly high, and one is significantly low)
			# need to skip these windows
			if (abs($mean1 - $mean2) > $min_mean_diff){			
				# print "\tMeans are too different ($mean1 vs $mean2)\n";
				next;
			}
			
			# Now ask whether the current window is better (higher t-score) than 
			# a pre-existing window? If so, remove pre-existing window from main hash
			if($t1 > $t2){
				delete($sig_windows{$window});
				# print "\tdeleting $start2-$end2\n";
				next;
			} 
			# can now add checks to see whether current window is worse (lower t-score)
			# than pre-existing values, in which case we will ignore it (by returning 0);
			elsif($t1 < $t2){
				# print "\tKeeping $start2-$end2\n";
				return(0);
			}
			# If t-values are the same (extremely unlikely), keep both windows
			elsif($t1 == $t2){
				return(0);
			}
			else{
				# shouldn't ever get here
				die "ERROR: ($start1-$end1 vs $start2-$end2)\n";
			}
	}
	return(1);
}

# compare two pairs of coordinates to detect overlap
sub check_for_overlap{

	my ($start1, $end1, $start2, $end2) = @_;

	# LHS overlap of first window (or same coords)
	if(($start1 <= $start2) and ($end1 >= $start2)){
		return(1);
#		return(0);
	} 
	# RHS overlap of first window (or same coords)
	elsif(($start1 >= $start2) and ($start1 <= $end2)){
		return(1);
#		return(0);
	} 
	# First window completely encloses second window
	elsif(($start1 <= $start2) and ($end1 >= $end2)){
		return(1);
	} 
	# First window is completely enclosed by second window
	elsif(($start1 > $start2) and ($end1 < $end2)){
		return(1);
#		return(0);
	} else{
		return(0); # no overlap
	}
}


sub t_test{
	my ($num1, $num2, $n1, $n2, $sum1, $sum2, $mean1, $mean2, $var1, $var2);
	
	# will receive references to 2 arrays of numbers
	($num1, $num2) = @_;
	 
	# n1, mean1, sum1, and var1 is the same for either method
	$n1 = scalar @$num1;
	$sum1 = sum(@$num1);		
	$mean1 = $sum1/$n1;
	$var1 = sum(map {($_ - $mean1) ** 2} @$num1) / ($n1 - 1);

	# n2, mean2, sum2, and var2 
	$n2 = scalar @$num2;
	$sum2 = sum(@$num2);
	$mean2 = $sum2/$n2;
	$var2 = sum(map {($_ - $mean2) ** 2} @$num2) / ($n2 - 1);	
	
	
	# change $var1 and $var2 to be some small non-zero value if actual value is zero
	# this can happen for some very small window sizes
	$var1 = 0.0001 if ($var1 == 0);
	$var1 = 0.0001 if ($var2 == 0);

	my $std_err_of_diff = sqrt( ( ((($n1-1)*$var1) + (($n2-1)*$var2)) / ($n1+$n2-2) ) * 
	                            ($n1+$n2) / ($n1*$n2)
							  );
		
	# calculate t-test score and return this plus both means
	my $t = ($mean1 - $mean2) / $std_err_of_diff;			
	return(abs($t), sprintf("%.3f", $mean1), $mean2, sprintf("%.4f", $var1), sprintf("%.4f", $var2));
}

# open shuffle file information (if present)
# read data and add to hash
sub read_shuffle_data{

	die "$shuffle_file does not exist" unless (-e $shuffle_file);
	
	warn "Processing shuffle file: $shuffle_file\n";
	open(my $in, "<", $shuffle_file) or die "Can't open $shuffle_file\n";

	while(my $line = <$in>){
		chomp($line);
		my ($window_size, $t) = split(/\t/, $line);
		$max_t{$window_size} = $t;
	}	
	close($in);
	
	# add zero values for large window sizes that may not have been covered in the shuffling
	for(my $i = 2; $i <= $win_size; $i++){
		if (not defined $max_t{$i}){
			$max_t{$i} = 0;
		}
	}
}


# read actual counts data and load to array
sub read_counts_file{
	
	warn "Processing counts data in $file\n";
	open(my $in, "<", $file) or die "Can't read from $file\n";

	while(my $line = <$in>){
		chomp($line);

	
		# skip header line if present
		next if ($line =~ m/^Chrom/);
	
		# skip blank lines
		next if ($line =~ m/^$/);
	
		# skip lines which don't have data (some lines might just have tabs)
		next unless ($line =~ m/^Chr/);

		my @line = split(/\t/, $line);
		my $chr = $line[0];
	
		# only want to process chromosomes 1 and 4
		next unless ($line =~ m/^Chr[14]/);
	
		# what to do with zero or undefined values? Set default to 2
		if (not defined $line[8] or $line[8] == 0){
			push(@{$read_counts{$chr}}, 2);
		} else{
			push(@{$read_counts{$chr}}, $line[8]);
		}
	}

	close($in);
}