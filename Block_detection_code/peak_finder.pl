#!/usr/bin/perl
#
# peak_finder.pl
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

my ($file, $chunk_size, $shuffle, $shuffle_file) = (undef, 1000, undef, undef);
my ($min_t, $min_window_size, $max_window_size, $skip)  = (2.576, 3, undef, undef);
my ($mean_mode, $min_mean_diff, $max_var, $help) = (1, 0.6, 0.04, undef);

my $usage = "
$0 --file <TSV file of read count data> --chunk_size <default = $chunk_size> 

optional arguments:
        --shuffle : produce a shuffle file after <x> shufflings of chunk data
		--shuffle_file : use a shuffle file
		--min_t : minimum t-test value to consider significant <default = $min_t>
		--min_window_size : minimum size of window <default = $min_window_size>
		--max_window_size : maximum size of window <default = undefined>
		--min_mean_diff : minimum difference expected between average copy number within and outside window <default = $min_mean_diff
		--max_var : maximum variance allowed withing a window <default = $max_var>
		--mean_mode : 1 or 2 (default = $mean_mode)
					mode = 1, compares mean of window with mean of all bins outside of window
					mode = 2, compares mean of window with mean of enclosing chunk (quicker for large chunks)
		--skip : skip rest of chromosome when the start chunk coordinate is greater than that specified
        --help
        
";

GetOptions (
	"file=s"            => \$file,
	"chunk_size=i"      => \$chunk_size,
	"shuffle=i"         => \$shuffle,
	"shuffle_file=s"    => \$shuffle_file,
	"min_t=f"           => \$min_t,
	"min_window_size=i" => \$min_window_size,
	"max_window_size=i" => \$max_window_size,
	"min_mean_diff=f"   => \$min_mean_diff,
	"max_var=f"         => \$max_var,
	"mean_mode=i"       => \$mean_mode,
	"skip=i"            => \$skip,
	"help"              => \$help,
);

die $usage if ($help);
die $usage if (not defined $file);
die "$shuffle should be positive integer" if (defined $shuffle and $shuffle == 0);

$max_window_size = $chunk_size / 2 if (not defined $max_window_size);
warn "Parameters:
--chunk_size = $chunk_size
--min_t = $min_t
--min_window_size = $min_window_size
--max_window_size = $max_window_size
--min_mean_diff = $min_mean_diff
--max_var = $max_var
--mean_mode = $mean_mode\n\n";

####################################
# Global variables
####################################

# keep track of maximum value of t observed at each window size
my %max_t;

# need to store count data in hash of arrays (primary key is chromosome name)
my %read_counts;

# will want to store all significant scoring windows found in each chunk
my %sig_windows;

# if $mean_mode == 2, we'll keep track of the mean and variance of each chunk
my ($chunk_mean, $chunk_var);


# read shuffled data from file (if specified)
read_shuffle_data() if ($shuffle_file);

# grab read count data from file
read_counts_file();

# check that $max_window_size doesn't exceed the number of data points available (at least for Chr1)
if ($max_window_size >= scalar @{$read_counts{'Chr1'}}){
	die "ERROR: Max window size ($max_window_size) >= number of data points\n";
}

# what we do next depends on whether we are shuffling data or not
# lots of data processing with shuffled data, or one run with the real (unshuffled) data
if ($shuffle){
	for (my $i = 1; $i <= $shuffle; $i++){
		warn "Shuffling data, run $i";
		process_chromosomes($i);
		
		# print out t-test scores after each shuffling
		my $output_file = "$file.shuffled.$chunk_size.$shuffle";
		open(my $out, ">", $output_file) or die "Can't open $output_file\n";
		foreach my $window (sort {$a cmp $b} keys %max_t){
			print $out "$window\t$max_t{$window}\n";
		}
		close($out);
	}
} else{
		warn "Looping over chromosome in chunks of size $chunk_size\n";
		print "Chr\tChunk\tWin:start-end\tWin_size\tType\tCopy_number\tStart_coord\tEnd_coord\tt\tMean1\tMean2\n";
		process_chromosomes(1);
}


exit;


####################################
#
# Subroutines
#
####################################


####################################
# Main loop, loop over chromosome in chunks
####################################

sub process_chromosomes{

	my ($run) = @_;
	
	# will keep track of final results per chunk, as there may be duplicate
	# significant windows, found from overlapping chunks, and we only want to
	# show one
	my %final_windows;

	# 1st loop, just loop over desired chromosomes worth of data
	foreach my $chr (qw(Chr1 Chr4)){

		# 2nd loop, going over chunks (v. large windows) across each chromosome
		for(my $chunk_start = 0; $chunk_start < @{$read_counts{$chr}} - $max_window_size; $chunk_start += ($max_window_size)){
			# do we want to use our shortcut to end early/
			if($skip){
				die "Chunk start $chunk_start exceeded by --skip $skip option, exiting early\n" if ($chunk_start > $skip);	
			}
			
			my $chunk_end = $chunk_start + $chunk_size - 1;
			$chunk_end = scalar(@{$read_counts{$chr}}) - 1  if ($chunk_end > scalar @{$read_counts{$chr}});

			# shuffle this chunk if we are in shuffle mode
			if($shuffle){
				@{$read_counts{$chr}}[$chunk_start..$chunk_end] = shuffle(@{$read_counts{$chr}}[$chunk_start..$chunk_end]);
			}

			warn "Run $run, looping over CHUNK: $chr:$chunk_start-$chunk_end\n";
			
			# calculate mean and variance of read counts in current chunk (if $mean_mode == 2)
			if ($mean_mode == 2){
				$chunk_mean = sum(@{$read_counts{$chr}}[$chunk_start..$chunk_end]) / $chunk_size;
				$chunk_var = sum(map {($_ - $chunk_mean) ** 2} @{$read_counts{$chr}}[$chunk_start..$chunk_end]) / ($chunk_size - 1);
			}
	
			# reset significant windows
			%sig_windows = ();

			# now process all windows in this chunk
			loop_over_chunk($chunk_start, $chunk_end, $chr);


			# final output for current chunk, in descending t-score order
			foreach my $window (sort {$sig_windows{$b}{t} <=> $sig_windows{$a}{t}} keys %sig_windows){
				my $t        = sprintf("%.2f", $sig_windows{$window}{t});
				my $win_size = $sig_windows{$window}{size};
				my $mean1    = sprintf("%.2f", $sig_windows{$window}{mean1});
				my $mean2    = sprintf("%.2f", $sig_windows{$window}{mean2});
		
				my ($start, $end)   = split(/-/, $window);

				# have we already seen this window from a previous, overlapping, chunk?
				# if so, just keep one window. They may have different t-scores based on 
				# mean of chunk (if using --mean_mode 2)
				next if (exists($final_windows{"$chr:$start-$end"}));

				# now format output for printing
				my $start_pos = sprintf("%.0f",($start * 1000));
				my $end_pos   = sprintf("%.0f",($end * 1000));

				my $win_coords = "$start-$end";

				# some significant windows will correspond to edge of chunks
				# rather than edges of true 2x/3x blocks within chunks. Want to flag these.
				# default is 'block' for whole block, otherwise edge-L, edge-R or possibly
				# edge-LR. Also, if the window is exactly half the chunk size, then this is 
				# not a true block (call this a 'fakeblock', we probably have not found the
				# true ends of these fakeblocks)

				my $type;
				
				if(($end - $start +1) == $max_window_size){
					$type = 'fakeblock';
				} elsif ($start == $chunk_start and $start == 0){
					$type = 'block';
				} elsif ($start == $chunk_start){
					$type = 'edge-R';
				} elsif ($end == $chunk_end){
					$type = 'edge-L';
				} else {
					$type = 'block';
				}
				
				# now calculate whether this is a 1x, 2x, or 3x block
				# this is a guesstimate
				my $copy_number;
				
				if ($mean1 < 2){
					$copy_number = "1x";
				} elsif ($mean1 < 3){
					$copy_number = "2x";
				} else{
					$copy_number = "3x";
				}
				
				print "$chr\t";
				print "$chunk_start-$chunk_end\t";
				print "$win_coords\t";
				print "$win_size\t";
				print "$type\t";
				print "$copy_number\t";
				print "$start_pos\t$end_pos\t";
				print "$t\t";
				print "$mean1\t$mean2\n";
		
				# log this window
				$final_windows{"$chr:$start-$end"} = 1;		
			}
		}
	}
}




sub loop_over_chunk{

	my ($start, $end, $chr) = @_;
		
	# start comparing windows of size 2 up to 'n / 2' (where n is the size of the current chunk)
#	for(my $win_size = $max_window_size; $win_size >= $min_window_size; $win_size--){
	for(my $win_size = $min_window_size; $win_size <= $max_window_size; $win_size++){

#		warn "$chr\tChunk:$start-$end\twin_size = $win_size\n";

		# generate initial array of read counts for first window that we inspect
		my @values_A;
		
		# if $mean_mode == 1, will compare @values_A to all *other* counts in current chunk
		# otherwise, take mean of entire chunk (no real difference for small window sizes,
		# but there will be a difference when window size is big and chunk size is small)
		my @values_B;
		
		# @values_B will be comprised of two regions (to the left and right of @values_A)
	 	my (@values_BL, @values_BR);
	 	
		# now loop along windows in current chunk (data from @{$read_counts{}} array)
		WINDOW: for(my $i = $start; $i <= $end - $win_size + 1; $i++){
			my ($min, $max) = ($i, $i + $win_size - 1);

			warn "$chr\tChunk:$start-$end\twin_size = $win_size\tWindow $min-$max\\n";

			# first window along chunk
			if($i == $start){
				# set default contents of @values_A, don't need to do anything else
				@values_A = @{$read_counts{$chr}}[$min..$max];
				if ($mean_mode == 1){
					@values_B = @{$read_counts{$chr}}[$max+1..$end]; 
					@values_BR = @values_B;			
				}
			} 
			# all other windows
			else{
				# get rid of one element from LHS of @values_A, and add a new element on RHS
				my $tmp = shift(@values_A);
				push(@values_A, ${$read_counts{$chr}}[$max]);	
				if ($mean_mode == 1){
				
					# add one value to left-hand array
					push(@values_BL, $tmp);
				
					# and remove one value from right-hand array
					shift(@values_BR);
			
					# form @values_B from the two separate arrays
					@values_B = (@values_BL, @values_BR);
				}
			}

			# calculate t value, do this slightly differently depending on value of $mean_mode
			my ($t, $mean1, $mean2, $var1);
			if ($mean_mode == 1){
				($t, $mean1, $mean2, $var1) = t_test(\@values_A, \@values_B);			
			} else{
				($t, $mean1, $mean2, $var1) = t_test(\@values_A, $chunk_mean, $chunk_var);			
			}

			# skip if we are not looking at a duplicated block 
			# these will have means of ~1.7, but we'll use a higher cut-off value
			# (based on observations)
#			next if ($mean1 < 2.2);
			
#			print "\t$win_size\t$min-$max\t$t\t$mean1\t$mean2\n" if ($t > $min_t);
			
			# if we are in -shuffle mode, do we have a new highest t_score at current window size, for current chunk?				
			if($shuffle){
				if(not defined $max_t{"$chr:$start-$end:$win_size"}){
					$max_t{"$chr:$start-$end:$win_size"} = $t;			
					my $rounded_t = sprintf("%.3f", $t);
#					warn "\tShuffle) $min-$max ($start-$end:$win_size) - new max t-value $rounded_t\n";				
	
				} elsif($t > $max_t{"$chr:$start-$end:$win_size"}){
					$max_t{"$chr:$start-$end:$win_size"} = $t;
					my $rounded_t = sprintf("%.3f", $t);
#					warn "\tShuffle) $min-$max ($start-$end:$win_size) - new max t-value $rounded_t\n";				
				}
				next WINDOW;
			}

			
			# skip windows if the window mean is too close to the rest-of-the-window mean
			# we are looking for cases where a window goes from 1x reads to 2x or 3x
			# so difference in means should really be at least 1.0, but we will be cautious
			next if (abs($mean1 - $mean2) <= $min_mean_diff);
			
			# only want to record windows which are statistically significant (P < 0.001 t = ~3.29)			
			if ($t >= $min_t){

				# if we have used a shuffled results file, we can now also ask whether
				# t-value could be expected by chance (ignore current t-value if this is the case)
				next WINDOW if ($shuffle_file and $t <= $max_t{"$chr:$start-$end:$win_size"});

				# is there too much variance inside this window?
				# e.g. when a block incorrectly spans a 2x block and some of a 3x block, variance
				# will be higher than a smaller block that captures just the 2x region
				next WINDOW unless ($var1 <= $max_var);

				# now check if this current window encloses a pre-existing window
				# which has a lower t-value score. If so, can remove the smaller window
				next WINDOW unless check_for_redundancy($min, $max, $mean1, $t);

				my $mean_diff = sprintf("%.2f", abs($mean1 - $mean2));
				warn "\t$chr\t$start-$end\t$min-$max\t$win_size\tblock\t2x\t${min}000\t${max}000\t$t\t$mean1\t$mean2\tmean_diff = $mean_diff\tvar = $var1\n";				
	
			
				# keep track of details of this significant window
				$sig_windows{"$min-$max"}{t}     = $t;
				$sig_windows{"$min-$max"}{size}  = $win_size;
				$sig_windows{"$min-$max"}{mean1} = "$mean1";
				$sig_windows{"$min-$max"}{mean2} = "$mean2";				
					
			} 				
		}
	}
}

# take current significant window and compare to other significant windows
# we have already seen. If the current window is bigger than an overlapping window
# *and* has a higher t-score, then we can ignore it

sub check_for_redundancy{

	my ($start1, $end1, $mean1, $t1) = @_;

	foreach my $window (sort {$a cmp $b} keys %sig_windows){
			my ($start2, $end2) = split(/-/, $window);
			my $t2 = $sig_windows{$window}{t};
			my $mean2 = $sig_windows{$window}{mean1};
			my $size = keys(%sig_windows);
			
			# first ask whether windows overlap (in any way) and skip to next 
			# comparison if windows don't overlap
			my $overlap = check_for_overlap($start1, $end1, $start2, $end2);
			next unless ($overlap);

			my ($s1, $s2) = ($end1-$start1+1, $end2-$start2+1);
#			print "Comparing $start1-$end1 (L=$s1, t=$t1) with $start2-$end2 (L=$s2, t=$t2)\n";
			
			# If two significant windows overlap, they might have opposing
			# means (i.e. one is significantly high, and one is significantly low)
			# need to skip these windows
			if (abs($mean1 - $mean2) > $min_mean_diff){			
#				print "\tMeans are too different ($mean1 vs $mean2)\n";
				next;
			}
			
			# Now ask whether the current window better (higher t-score) than 
			# a pre-existing window? If so, remove pre-existing window from main hash
			if($t1 > $t2){
				delete($sig_windows{$window});
#				print "\tdeleting $start2-$end2\n";
				next;
			} 
			# can now add checks to see whether current window is worse (lower t-score)
			# than pre-existing values, in which case we will ignore it (by returning 0);
			elsif($t1 < $t2){
#				print "\tKeeping $start2-$end2\n";
				return(0);
			}
			# If t-values are the same, just keep largest window
			elsif($t1 == $t2){

				# current window is larger
				if(($end1 - $start1) > ($end2 -$start2)){
					die "$start1-$end1 vs $start2-$end2\n";
					delete($sig_windows{$window});
					next;
				} 
				# pre-existing window is larger
				elsif(($end1 - $start1) < ($end2 -$start2)){
					die "$start1-$end1 vs $start2-$end2\n";
					return(0);
				} 
				else {
					# same length, so no way to choose
					# will end up returning 1 (below) so both windows will be kept
				}
			
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
	my ($num1, $num2, $n1, $n2, $sum1, $sum2, $mean1, $mean2, $var1, $var2, $chunk_mean, $chunk_var);
	
	# will receive references to 2 arrays of numbers ($mean_mode == 1)
	# or 1 array of numbers plus mean and variance for entire chunk ($mean_mode == 2)
	if ($mean_mode == 1){
	 	($num1, $num2) = @_;
	} else {
		# using whole chunk for comparison
	 	($num1, $chunk_mean, $chunk_var) = @_;
	}
	 
	# n1, mean1, sum1, and var1 is the same for either method
	$n1 = scalar @$num1;
	$sum1 = sum(@$num1);		
	$mean1 = $sum1/$n1;
	$var1 = sum(map {($_ - $mean1) ** 2} @$num1) / ($n1 - 1);

	# n2, mean2, sum2, and var2 depends on which $mean_mode we are using
	if ($mean_mode == 1){
		$n2 = scalar @$num2;
		$sum2 = sum(@$num2);
		$mean2 = $sum2/$n2;	
		$var2 = sum(map {($_ - $mean2) ** 2} @$num2) / ($n2 - 1);	
		
	} else {
	 	$n2 = $chunk_size;
		$mean2 = $chunk_mean;
		$var2 = $chunk_var;		
	}	

	# change $var1 and $var2 to be some small non-zero value if actual value is zero
	# this can happen for some very small window sizes
	$var1 = 0.0001 if ($var1 == 0);
	$var1 = 0.0001 if ($var2 == 0);

	my $std_err_of_diff = sqrt( ( ((($n1-1)*$var1) + (($n2-1)*$var2)) / ($n1+$n2-2) ) * 
	                            ($n1+$n2) / ($n1*$n2)
							  );
		
	# calculate t-test score and return this plus both means
	my $t = ($mean1 - $mean2) / $std_err_of_diff;			
	return(abs($t), sprintf("%.3f", $mean1), $mean2, sprintf("%.4f", $var1));
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
	for(my $i = 2; $i <= $max_window_size; $i++){
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