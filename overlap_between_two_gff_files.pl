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

###################################################### 
#              Command-line options                  #
###################################################### 

my ($junction_gff, $feature_gff, $bp, $shuffles, $help) = (undef, undef, 10000, 100, undef);

my $usage = "$0 
Mandatory arguments:
--junction_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>

Optional arguments:
--bp <how many bp to extract from inside/outside region (default = $bp)>
--shuffles <how many shuffling iterations to perform (default = $shuffles)>
--help 
	
";

GetOptions (
	"junction_gff=s" => \$junction_gff,
	"feature_gff=s"  => \$feature_gff,
	"bp=i"           => \$bp,
	"shuffles=i"     => \$shuffles,
	"help"           => \$help,
);


die $usage if ($help);
die $usage if (not defined $junction_gff);
die $usage if (not defined $feature_gff);


# Need a hash to store sizes of chromosomes 
my %chr_sizes;
$chr_sizes{'1'} = 30427671;
$chr_sizes{'2'} = 19698289;
$chr_sizes{'3'} = 23459830;
$chr_sizes{'4'} = 18585056;
$chr_sizes{'5'} = 26975502;



# read junction coordinates and store in hash
my %junctions;
my %junction_counts;

open (my $in, "<", $junction_gff) or die "Can't read $junction_gff\n";

my $tmp = 0;
while(my $line = <$in>){
    my ($chr, undef, undef, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	
	
	# coordinate that we use depends on whether this was the left or right edge
	my $coord;
	$coord = $s if ($comment =~ m/edge=R/);
	$coord = $e if ($comment =~ m/edge=L/);
	
	# now define a range around this coordinate based on value of $bp
	my ($min, $max) = ($coord - $bp/2, $coord + $bp/2);
	$junctions{$chr}{"$min-$max"}{start} = $min;
	$junctions{$chr}{"$min-$max"}{end}   = $max;
	$junction_counts{$chr}++;
	
	my $distance;
	if($comment =~ m/edge=R/){
		$distance = $s - $tmp;	
		$tmp = $s;
	} else{
		$distance = $e - $tmp;
		$tmp = $e;
	}
	print "Distance = $distance\t$line";
}

close($in);


# now loop over master GFF file and count hits inside and outside junctions
 
my %overlaps;

open ($in, "<", $feature_gff) or die "Can't read $feature_gff\n";

OUTER: while(my $line = <$in>){
    my ($chr, undef, $feature, $s, $e, undef, undef, undef, undef) = split(/\t/, $line);	
	
	# skip chromosome features
	next if $feature eq 'chromosome';

	# does current chromosome feature overlap a junction?
	foreach my $key (keys %{$junctions{$chr}}){
		my $start = $junctions{$chr}{$key}{start};
		my $end   = $junctions{$chr}{$key}{end};
		
		# feature either overlaps with edge of junction region (overlap = 1)
		# or is completely self contained inside junction region (overlap = 2)
		my $overlap = 0;
		if ($s >= $start and $e <= $end){
			$overlap = 2;
		} elsif(($s >= $start and $s <= $end) or ($e >= $start and $e <= $end)){
			$overlap = 1;
		}

		# track overlap (or lack of overlaps)
		$overlaps{$chr}{$feature}{outside}++  if ($overlap == 0);
		$overlaps{$chr}{$feature}{spanning}++ if ($overlap == 1);
		$overlaps{$chr}{$feature}{inside}++   if ($overlap == 2);

		if($overlap){
			print "Overlap $overlap\tJunction:$key vs $feature $s-$e\t$line\n";
			next OUTER;
		}
    }
}

close($in);


__END__
# store the final ratio from the real data for comparisons when shuffling
my $real_ratio;

my ($in_region, $out_region);
my ($in_region_ori_bp, $out_region_ori_bp);

# Also want to count absolute number of replication origins that overlap
my ($in_region_ori_n, $out_region_ori_n);



print "\n############################################\n";
print "Feature + Region + bp = $feature $region $bp\n";
print "############################################\n\n";



# represent a chromosome as a string of dashes
my $chr_seq = '-' x $chr_sizes{"$chromosome"};



########################################################################
# 2) process ori data, mark bases in virtual sequence that are replication origins
########################################################################

# count oris across whole chromosome
my ($ori_total_bp, $ori_n) = (0, 0);

open(my $in, "<", "$gff") or die "Can't read from $gff\n";
while(my $line = <$in>){
	my ($chr, undef, undef, $start, $end, undef, undef, undef, undef) = split(/\t/, $line);
	
	# only want to consider matches to target chromosome 
	next unless $chr eq "Chr$chromosome";

	# will avoid the end of Chr1 (where tailswap with Chr4 occurs)
	next if $start > $chr_sizes{"$chromosome"};
	
	my $ori_length = $end - $start + 1;
	$ori_total_bp += $ori_length;
	$ori_n++;
	
	substr($chr_seq, $start, $ori_length) = ("o" x $ori_length);
}
close($in);

# now want to keep the respective ori & non-ori regions in separate arrays
# for when we will shuffle data later on
my @oris     = split(/-+/, $chr_seq);
my @non_oris = split(/o+/, $chr_seq);

# first element of @oris will be empty, so just remove
shift(@oris);


####################################
# Main loop for shuffling
####################################


# keep track of how many times final ratio is exceeded
my $ratio_exceeded = 0;

for (my $i = 0; $i <= $shuffles; $i++){
	
	# do we need to shuffle $chr_seq
	if ($i > 0){
		print "Shuffle $i\n";
		$chr_seq = shuffle_chr();	
	}
	find_overlaps($i);
}


# how often did we beat the real ratio replication origins per Mbp (inside vs outside)
# in the shuffled datasets?
print "Ratio in real data was exceeded $ratio_exceeded in $shuffles shuffles\n";


exit(0);






####################################
#
# Subroutines
#
####################################

# want to make randomized version of chromosome with different locations
# of replication origins, but preserving size

sub shuffle_chr{
	
	# make copies of main arrays with ori & non-ori data and shuffle them
	my @oris2     = shuffle(@oris);
	my @non_oris2 = shuffle(@non_oris);
	my $shuffled_seq;

# alternative way of shuffling
#	my @combined = shuffle(@oris2, @non_oris2);
#	$shuffled_seq = join('', @combined);	
#	return($shuffled_seq);
	
	# first take a non-ori and add to $shuffled_seq
	$shuffled_seq .= shift(@non_oris2);
	
	while(@oris2){

#		print "There are ", scalar @non_oris2, " non-ori elements and ", scalar @oris2, " ori elements\n";

		# first take a non-ori and add to $shuffled_seq
		$shuffled_seq .= shift(@non_oris2);
				
		# now add an ori		
		$shuffled_seq .= shift(@oris2);
	}

	return($shuffled_seq);
}



sub find_overlaps{
	my ($shuffle) = @_;

	########################################################################
	# 3) Parse block data looking for 2x/3x blocks etc
	########################################################################

	# want to count bases inside and outside of 2x/3x blocks  
	# (or whatever is specified by $feature & $region)
	($in_region, $out_region) = (0, 0);

	# likewise want to count bases inside and outside of 2x/3x blocks that overlap 
	# with replication origins. 
	($in_region_ori_bp, $out_region_ori_bp) = (0, 0);

	# Also want to count absolute number of replication origins that overlap
	($in_region_ori_n, $out_region_ori_n) = (0, 0);


	# Main loop to go over Han's block data
	open($in, "<", "$blocks") or die "Can't read from $blocks\n";
	while(my $line = <$in>){
		my ($chr, $block_id, $han_id, $start, $end, $copy_number, undef, undef) = split(/\t/, $line);

		# only considering Chr1 for now
#		next unless $chr eq 'Chr3';
		next unless $chr eq 'Chr1';

		# Only consider regions that are inside/outside a 2x/3x block
		next unless ($copy_number eq $feature);

		my $block_length = $end - $start + 1;

		if ($region eq 'I'){ # INSIDE A REGION
		
				# if block is not twice the length of $bp, just use 1/2 block length
				if(($bp * 2) > $block_length){
					# now count how many bases in this block are also in replication origins
				
					my ($ori_count, $ori_bp) = process_region($start, $block_length);
					print "$han_id\t$block_id ($block_length bp)\tI\t$ori_count\tORI bp = $ori_bp\n" if ($shuffle == 0); 

				} else{
					my ($ori_count1, $ori_bp1) = process_region($start, $bp);
					my ($ori_count2, $ori_bp2) = process_region($end - $bp, $bp);

					print "$han_id\t$block_id ($block_length bp)\tI\t$ori_count1 + $ori_count2\tORI bp = $ori_bp1 + $ori_bp2\n" if ($shuffle == 0);
				
				}
		} elsif ($region eq 'O'){ # OUTSIDE
			my ($ori_count1, $ori_bp1) = process_region($start - $bp, $bp);
			my ($ori_count2, $ori_bp2) = process_region($end, $bp);

			print "$han_id\t$block_id ($block_length bp)\tO\t$ori_count1 + $ori_count2\tORI bp = $ori_bp1 + $ori_bp2\n" if ($shuffle == 0);		
		
		
		
		} elsif ($region eq 'A'){ # ALL of a region
			my ($ori_count, $ori_bp) = process_region($start, $block_length);				
			print "$han_id\t$block_id ($block_length bp)\tA\t$ori_count\tORI bp = $ori_bp\n" if ($shuffle == 0);

	
		} elsif ($region =~ m/[IO]{2}/){ # INSIDE and OUTSIDE
			# if block is not twice the length of $bp, just use 1/2 block length
			if(($bp * 2) > $block_length){
				# now count how many bases in this block are also in replication origins
		
				my ($ori_count, $ori_bp) = process_region($start - $bp, $block_length + $bp + $bp);
				print "$han_id\t$block_id ($block_length bp)\tIO*\t$ori_count\tORI bp = $ori_bp\n" if ($shuffle == 0);

			} else{
				my ($ori_count1, $ori_bp1) = process_region($start - $bp, $bp + $bp);
				my ($ori_count2, $ori_bp2) = process_region($end - $bp, $bp + $bp);

				print "$han_id\t$block_id ($block_length bp)\tIO\t$ori_count1 + $ori_count2\tORI bp = $ori_bp1 + $ori_bp2\n" if ($shuffle == 0);
			}
	
		} elsif ($region =~ m/[AO]{2}/){ # ALL REGION and OUTSIDE

			my ($ori_count, $ori_bp) = process_region($start - $bp, $block_length + $bp + $bp);				
			print "$han_id\t$block_id ($block_length bp)\tAO\t$ori_count\tORI bp = $ori_bp\n" if ($shuffle == 0);
	
		} else{
			die "Should never reach this else statement!\n";
		}
		

	
	}
	close($in);

	# calculate bases outside of 2x/3x regions and number of replication origins
	$out_region = $chr_sizes{$chromosome} - $in_region;
	($out_region_ori_bp) = $chr_seq =~ tr/o/o/;

	# count oris in this block by splitting at runs of Ns (need to subtract 1)
	$out_region_ori_n = split(/o+/, $chr_seq) - 1;

	# work out replication origins per kb
	my $in_ori_per_mb  = sprintf("%.1f", ($in_region_ori_n  / $in_region)  * 1000000);
	my $out_ori_per_mb = sprintf("%.1f", ($out_region_ori_n / $out_region) * 1000000);
	my $chr_ori_per_mb = sprintf("%.1f", ($ori_n / $chr_sizes{$chromosome}) * 1000000);

	my $percent_ori                = sprintf("%.2f", ($ori_total_bp / $chr_sizes{$chromosome}) * 100);
	my $percent_ori_inside_region  = sprintf("%.2f", ($in_region_ori_bp / $in_region) * 100);
	my $percent_ori_outside_region = sprintf("%.2f", ($out_region_ori_bp / $out_region) * 100);

	my $final_ratio = sprintf("%.4f", (($in_region_ori_n  / $in_region)  * 1000000) / (($out_region_ori_n / $out_region) * 1000000));
	
	($real_ratio = $final_ratio) if ($shuffle == 0);

	# have we beaten $real_ratio?
	if(($shuffle > 0) and ($final_ratio > $real_ratio)){
		$ratio_exceeded++;
		print "\tReal ratio exceeded: $final_ratio\n";
	}

	####################################
	# Final output - only for 1st iteration
	####################################

	if($shuffle == 0){
#	if($shuffle <= $shuffles){
		print "\n";
		print "Total number of replication origins on chr1: $ori_n\n";
		print "Total number of bases in replication origins on chr1: $ori_total_bp (%$percent_ori)\n";
		print "Number of replication origins per 1 Mbp of chr1: $chr_ori_per_mb\n\n";

		print "INSIDE $feature$region regions\n";
		print "-------------------\n";
		print "Total length (bp): $in_region\n";
		print "Number of replication origins: $in_region_ori_n\n";
		print "Sum length of regions that overlap replication origins (bp): $in_region_ori_bp (%$percent_ori_inside_region)\n";
		print "Number of replication origins per 1 Mbp of matching region: $in_ori_per_mb\n\n";


		print "OUTSIDE $feature$region regions\n";
		print "-------------------\n";
		print "Total length (bp): $out_region\n";
		print "Number of replication origins: $out_region_ori_n\n";
		print "Sum length of regions that overlap replication origins (bp): $out_region_ori_bp (%$percent_ori_outside_region)\n";
		print "Number of replication origins per 1 Mbp of matching region: $out_ori_per_mb\n\n";

		print "Inside/outside ratio of number of replication origins per Mbp: $final_ratio\n\n\n\n";
	}

}

sub process_region{

	my ($start, $length) = @_;

	# extract region from $chr_seq
	my $region = substr($chr_seq, $start, $length);

		
	# increment size of $in_region by length of current region
	$in_region += $length;	
	
	# count oris in this block by splitting at runs of Ns (need to handle
	# windows where the entire region is replication origin differently).
	my $ori_count = 0;
	if ($region =~ m/^o+$/){
		$ori_count = 1;
	} else {
		$ori_count = split(/o+/, $region) - 1;
	}

	# add to global count of how many replication origins we've seen
	$in_region_ori_n += $ori_count;

	# count replication origin bases in current region
	my ($ori_bp) = $region =~ tr/o/-/;

	# add to global count of how many bp of replication origins we've seen
	$in_region_ori_bp += $ori_bp;

	# at the same time convert 'o' to '-' in sequence representation so that
	# any o's that remain in $chr_seq will all be outside our target blocks
	# just need to use $block and replace that region of $chr
	substr($chr_seq, $start, $length, $region);

	return($ori_count, $ori_bp);
}


sub check_command_line_options{


}

