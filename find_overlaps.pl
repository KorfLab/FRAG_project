#!/usr/bin/perl
#
# find_overlaps.pl
#
# A script to find overlaps between known origins of replication in A. thaliana
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

my ($gff, $blocks, $feature, $region, $bp, $shuffles, $chromosome, $help) = (undef, undef, undef, 'A', 1000, 0, undef, undef);

check_command_line_options();

# Key to $region
# A = extract all bases inside $feature region
# I = extract just $bp bases INSIDE feature (or all of feature if length < $bp)
# O = extract $n bases from OUTSIDE feature
# IO = extract $n bases from inside and outside feature



# Need a hash to store sizes of chromosomes (though only interested in chr1 & chr3 for now)
# can't deal with chr4 translocation on end of chr1, so will make arbitrary cut-off
# the reason why we have the '1' is because we'll be using an array later and this
# is a fudge to deal with working with zero-indexing 
my %chr_sizes;
$chr_sizes{'1'} = 27700001;
$chr_sizes{'3'} = 23500000;


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

	my $usage = "

	$0 
	----------------------

	Mandatory arguments:
	--gff <gff file of replication origins> 
	--blocks <TSV file of Han's 2x/3x block coordinates>
	--feature <2x or 3x>
	--chromosome <which chromosome to match (1 or 3)>

	Optional arguments:
	--region <what type of regions to match: A,I,O,AO,IO (All, Inside, Outside, All+Outside, Inside+Outside), default = $region>
	--bp <how many bp to extract from inside/outside region (default = $bp)>
	--shuffles <how many shuffling iterations to perform (default = $shuffles)>
	--help 
		
	";

	GetOptions (
		"gff=s"        => \$gff,
		"blocks=s"     => \$blocks,
		"feature=s"    => \$feature,
		"region=s"     => \$region,
		"bp=i"         => \$bp,
		"shuffles=i"   => \$shuffles,
		"chromosome=i" => \$chromosome,
		"help"         => \$help,
	);


	die $usage if ($help);
	die $usage if (not defined $gff);
	die $usage if (not defined $blocks);
	die $usage unless ($feature eq "2x" or $feature eq "3x");
	die $usage if (not defined $region);
	die $usage unless ($chromosome =~ m/^[1-5]$/);
	die $usage unless ($region eq "A" or $region eq "I" or $region eq "O" or $region eq "AO" or $region eq "IO");

}

####################################
#              DATA
####################################

# GFF replication origin data looks like:
# Chr1    GutierrezLab    DNA_replication_origin  814096  814844  .       .       .       ID=ori1-0140

# Duplicated/triplicated data looks like this:
# Chr1    BLOCK0001       01a     1       205500  2x      .
# Chr1    BLOCK0002       01b     32500   86500   3x      Contained in BLOCK0001

# Need to look for overlaps, report on significance








# 4) Count total bases in all 3x blocks

# 5) Count total bases *outside* 3x blocks

# 6) Count total ori bases in 3x blocks

# 7) Count total ori bases *outside* 3x blocks

# 8) Calculate % of bases in 3x blocks that are oris, compare to % of bases outside 3x blockst that are in oris

# 9) Repeat for 2x blocks

# 10) Repeat for 'inner edges' of 2x/3x blocks (1 kbp, 10 kbp?) 

# 11) Repeat for flanking regions of 2x/3x blocks (1 kbp either side?)

# 12) Shuffle ORI locations randomly and see if we can easily beat observed results
