#!/usr/bin/perl
#
# breakpoint_shuffler.pl
#
# A script to shuffle data to work out what t-test value could be expected by chance
# to be used with breakpoint_finder.pl script
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

my ($file, $win_size, $shuffle, $help) = (undef, 20, 1000);

my $usage = "
$0 --file <TSV file of read count data> 

optional arguments:
        --win_size <default = $win_size> 
        --shuffle : produce a shuffle file after <n> shufflings of chunk data <default = $shuffle>
        --help
";

GetOptions (
	"file=s"            => \$file,
	"win_size=i"        => \$win_size,
	"shuffle=i"         => \$shuffle,
	"help"              => \$help,
);

die $usage if ($help);
die $usage if (not defined $file);
die "$shuffle should be positive integer" if (defined $shuffle and $shuffle == 0);

####################################
# Global variables
####################################

# keep track of maximum value of t observed at each window size
my %max_t;

# need to store count data in hash of arrays (primary key is chromosome name)
# grab read count data from file
my %read_counts;
read_counts_file();


for (my $i = 1; $i <= $shuffle; $i++){
	warn "Shuffling data, run $i";
	process_chromosomes($i);
	
	# print out t-test scores after each shuffling
	my $output_file = "$file.shuffled.$win_size.$shuffle";
	open(my $out, ">", $output_file) or die "Can't open $output_file\n";
	foreach my $window (sort {$a cmp $b} keys %max_t){
		print $out "$window\t$max_t{$window}\n";
	}
	close($out);
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
	
	# 1st loop, just loop over desired chromosomes worth of data
	foreach my $chr (qw(Chr1 Chr4)){

		# shuffle read counts on this chromosome
		@{$read_counts{$chr}} = shuffle(@{$read_counts{$chr}});

		# generate initial array of read counts for first window that we inspect
		my @values_A;
	
		# Will compare @values_A to all *other* counts in current chunk
		my @values_B;
		
		# now loop along windows in current chromosome (data from @{$read_counts{}} array)
		for(my $i = 0; $i <= @{$read_counts{$chr}} - $win_size; $i++){

			my ($min, $max) = ($i, $i + $win_size - 1);			

			my $mid_point_A = $min + ($win_size / 2) - 1;		
			my $mid_point_B = $mid_point_A + 1;

			# is this the first window along the chromosome?
			if($i == 0){
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
			my ($t, $mean1, $mean2) = t_test(\@values_A, \@values_B);			
		 
	#		print "$chr\t$min-$max\t$min-$mid_point_A\t$mid_point_B-$max\tt=$t\tmeans = $mean1 vs $mean2\tvars = $var1 vs $var2\n";
			if(not defined $max_t{"$chr:$win_size"}){
 				$max_t{"$chr:$win_size"} = $t;			
 				my $rounded_t = sprintf("%.3f", $t);
				warn "\tShuffle $run) $chr:$win_size new max t-value $rounded_t\n";				
 
 			} elsif($t > $max_t{"$chr:$win_size"}){
 				$max_t{"$chr:$win_size"} = $t;
 				my $rounded_t = sprintf("%.3f", $t);
				warn "\tShuffle $run) $chr:$win_size new max t-value $rounded_t\n";				
 			}
		} 				
	}
}




sub t_test{
	my ($num1, $num2, $n1, $n2, $sum1, $sum2, $mean1, $mean2, $var1, $var2,);
	
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