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


####################################
# Set up chromosome data
####################################

# Need a hash to store sizes of chromosomes 
my %chr_sizes;
$chr_sizes{'Chr1'} = 30427671;
$chr_sizes{'Chr1'} = 27700001; # a fudge to deal with tail swap
$chr_sizes{'Chr2'} = 19698289;
$chr_sizes{'Chr3'} = 23459830;
$chr_sizes{'Chr4'} = 18585056;
$chr_sizes{'Chr5'} = 26975502;

# where does the tail swap begin on Chr4?
my $pre_tailswap_length = 16500000;

# will end up representing each chromosome as a string of dashes
my %chr_seqs;
# create fake chromosome sequences
foreach my $chr (qw(Chr1 Chr2 Chr3 Chr4 Chr5)){
	$chr_seqs{$chr} = '-' x $chr_sizes{$chr};
	
	# have to do something different for Chromosome 4 as we are only interested in the 
	# tail swap region. So mask out first part of chromosome with a different character
	if($chr eq 'Chr4'){
		substr($chr_seqs{$chr}, 0, $pre_tailswap_length) = ("x" x $pre_tailswap_length);		
	}
}


################################################################
# read junction coordinates and represent in virtual sequences
################################################################

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
	
	# mask where junctions are in chromosome
#	print "Substituting in chr $chr from $min for $bp bp\n";
	substr($chr_seqs{$chr}, $min, $bp) = ("J" x $bp);

}

close($in);

##########################################
# Main loop over each possible feature
##########################################

# track all results by final ratio
my %results_by_ratio;

# not the quickest way to do this, but ensures we don't run out of memory

my @features = qw(CDS DNA_replication_origin exon five_prime_UTR gene mRNA miRNA ncRNA protein pseudogene pseudogenic_exon pseudogenic_transcript satellite snoRNA tRNA three_prime_UTR transposable_element transposable_element_gene transposon_fragment);

foreach my $desired_gff_feature (@features){

	# create copies of virtual chromosome sequences
	my %tmp_chr_seqs = %chr_seqs;
 
	open ($in, "<", $feature_gff) or die "Can't read $feature_gff\n";

	while(my $line = <$in>){
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, undef) = split(/\t/, $line);	

		# only want to look at one feature at a time
		next unless $feature eq $desired_gff_feature;
	
		# temporarily skip if not chr1 or chr4?
		next unless (($chr eq 'Chr1') or ($chr eq 'Chr4'));
	
		# skip tailswap regions of Chr1 and Chr4
		next if ($chr eq 'Chr1' and $s > $chr_sizes{'Chr1'});
		next if ($chr eq 'Chr4' and $s < $pre_tailswap_length);
	
		# mask where feature occurs in virtual chromosome sequence
		my $length = $e - $s + 1;
		substr($tmp_chr_seqs{$chr}, $s, $length) = ("o" x $length);
	}
	close($in);
	
	# now compare patterns in original virtual sequence (just junctions)
	# and tmp virtual sequence (masked with feature)
	my $output_text = "$desired_gff_feature\t";

	my $original_junction_bp      = $chr_seqs{'Chr1'} =~ tr/J/J/;
	$original_junction_bp        += $chr_seqs{'Chr4'} =~ tr/J/J/;

	my $original_non_junction_bp  = $chr_seqs{'Chr1'} =~ tr/-/-/;
	$original_non_junction_bp    += $chr_seqs{'Chr4'} =~ tr/-/-/;

	my $remaining_junction_bp     = $tmp_chr_seqs{'Chr1'} =~ tr/J/J/;
	$remaining_junction_bp       += $tmp_chr_seqs{'Chr4'} =~ tr/J/J/;

	my $remaining_non_junction_bp = $tmp_chr_seqs{'Chr1'} =~ tr/-/-/;
	$remaining_non_junction_bp   += $tmp_chr_seqs{'Chr4'} =~ tr/-/-/;


	my $feature_overlapping_junction     = $original_junction_bp - $remaining_junction_bp;
	my $feature_overlapping_non_junction = $original_non_junction_bp - $remaining_non_junction_bp;

	# likely to have some zero counts for bases overlapping junction regions, so handle accordingly
	my ($percent_overlapping_junction, $percent_overlapping_non_junction);
	
	if ($feature_overlapping_junction == 0){
		$percent_overlapping_junction = "0";
	} else{
		$percent_overlapping_junction = ($feature_overlapping_junction / $original_junction_bp) * 100;
	}
	
	if ($feature_overlapping_non_junction == 0){
		$percent_overlapping_non_junction = "0";
	} else{
		$percent_overlapping_non_junction = ($feature_overlapping_non_junction / $original_non_junction_bp) * 100;
	}
	
	my $ratio = sprintf("%.2f", $percent_overlapping_junction / $percent_overlapping_non_junction);
	$percent_overlapping_junction     = sprintf("%.2f", $percent_overlapping_junction);
	$percent_overlapping_non_junction = sprintf("%.2f", $percent_overlapping_non_junction);
		
	$output_text .= "$feature_overlapping_junction/$original_junction_bp bp (%$percent_overlapping_junction)\t";
	$output_text .= "$feature_overlapping_non_junction/$original_non_junction_bp bp (%$percent_overlapping_non_junction)";

	$results_by_ratio{$output_text} = $ratio;
}

print "Ratio\tGFF_feature\tInside_junction_overlap\tOutside_junction_overlap\n";
foreach my $result (sort {$results_by_ratio{$b} <=> $results_by_ratio{$a}} keys %results_by_ratio){
	print "$results_by_ratio{$result}\t$result\n";
}



