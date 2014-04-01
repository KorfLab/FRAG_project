#!/usr/bin/perl
#
# nearest_feature.pl
#
# A script to find the average nearest distance of various genome features
# to genome breakpoints from our chromosome shattering lines
# will ultimately combine this script with overlap_between_two_gff_files.pl
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

my ($junction_gff, $feature_gff, $help, $verbose) = (undef, undef, undef, undef);

my $usage = "$0 
Mandatory arguments:
--junction_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>

Optional arguments:
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"junction_gff=s" => \$junction_gff,
	"feature_gff=s"  => \$feature_gff,
	"help"           => \$help,
	"verbose"        => \$verbose,
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

my %breakpoints;

# first store details of breakpoints
read_junction_data();

my @features = qw(CDS DNA_replication_origin exon five_prime_UTR gene mRNA miRNA ncRNA protein pseudogene pseudogenic_exon pseudogenic_transcript satellite snoRNA tRNA three_prime_UTR transposable_element transposable_element_gene transposon_fragment);

open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

while(my $line = <$in>){
	my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

	# skip if not chr1 or chr4
	next unless (($chr eq 'Chr1') or ($chr eq 'Chr4'));

	# skip tailswap regions of Chr1 and Chr4
	next if ($chr eq 'Chr1' and $s > $chr_sizes{'Chr1'});
	next if ($chr eq 'Chr4' and $s < $pre_tailswap_length);

	# only want to consider features on our list
	next unless $feature ~~ @features;

	# want to separate protein-coding gene features from non-protein coding gne
	if($feature eq 'gene'){
		if ($comment =~ m/protein_coding_gene/){
			$feature = "protein_coding_gene";
		} else{
			$feature = "non_protein_coding_gene";		
		}
	}
	
	foreach my $breakpoint (keys %breakpoints){
		my $b_pos = $breakpoints{$breakpoint}{pos};
		my $b_chr = $breakpoints{$breakpoint}{chr};
		
		# only want to compare genomic features on same chromosome as breakpoint
		next unless ($chr eq $b_chr);

		# work out far breakpoint positions is from feature
		my $distance_to_breakpoint;
		if (abs($s - $b_pos) < abs($e - $b_pos)){
			$distance_to_breakpoint = abs($s - $b_pos);		
		} else{
			$distance_to_breakpoint = abs($e - $b_pos);				
		}

		
		# have we seen this feature before? If not add new key to hash
		if (not defined ($breakpoints{$breakpoint}{features}{$feature})){
			$breakpoints{$breakpoint}{features}{$feature} = $distance_to_breakpoint;
		} else {
			# is this closer than existing distance?
			if($distance_to_breakpoint < $breakpoints{$breakpoint}{features}{$feature}){
				$breakpoints{$breakpoint}{features}{$feature} = $distance_to_breakpoint;
#				print "$b_chr:$b_pos\tD=$distance_to_breakpoint\t$s-$e\t$feature\t$comment" if ($verbose and $breakpoint == "33" and $distance_to_breakpoint < 10000);

			}
		}
	}
}
close($in);

my %averages;
my $three_count = 0;
my $five_count = 0;

# Now loop through results to work out average distances
foreach my $breakpoint (keys %breakpoints){
	my $b_pos = $breakpoints{$breakpoint}{pos};
	my $b_chr = $breakpoints{$breakpoint}{chr};
	print "$breakpoint $b_chr:$b_pos\n" if ($verbose);
	foreach my $feat (keys %{$breakpoints{$breakpoint}{features}}){
		my $distance = $breakpoints{$breakpoint}{features}{$feat};
		
		if($verbose){
			my $gene_type = "";
			if ($feat eq 'protein_coding_gene' and $distance == $breakpoints{$breakpoint}{features}{'three_prime_UTR'}){
				$gene_type = "=3";
				$three_count++;
			} elsif($feat eq 'protein_coding_gene' and $distance == $breakpoints{$breakpoint}{features}{'five_prime_UTR'}){
				$gene_type = "=5";			
				$five_count++;
			}
			print "\t$feat$gene_type $distance\n";
		} 
		
		# store average distance in hash as a total of all distances for that feature
		$averages{$feat} += $distance;
	}
}

if($verbose){

	print "=3 count = $three_count\n";
	print "=5 count = $five_count\n";
}
my $n = keys %breakpoints;

foreach my $feature (sort {$averages{$a} <=> $averages{$b}} keys %averages){
	my $final_average = sprintf("%.2f", $averages{$feature} / $n);
	print "$feature\t$final_average\n";
}

###############################################################
# read junction coordinates and represent in virtual sequences
###############################################################

sub read_junction_data{
	open (my $in, "<", $junction_gff) or die "Can't read $junction_gff\n";

	my $breakpoint_count = 0;

	while(my $line = <$in>){

		my ($chr, undef, undef, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	
		$breakpoint_count++;
	
		my $mid_point = abs(int(($e - $s) /2 + $s));	
	
		# store details
		$breakpoints{$breakpoint_count}{chr} = $chr;
		$breakpoints{$breakpoint_count}{pos} = $mid_point;
	}
	close($in);
}





__END__ 
##########################################
# Main loop over each possible feature
##########################################

# not the quickest way to do this, but ensures we don't run out of memory

my @features = qw(CDS DNA_replication_origin exon five_prime_UTR gene mRNA miRNA ncRNA protein pseudogene pseudogenic_exon pseudogenic_transcript satellite snoRNA tRNA three_prime_UTR transposable_element transposable_element_gene transposon_fragment);

foreach my $desired_gff_feature (@features){

	# create copies of virtual chromosome sequences
	my %tmp_chr_seqs = %chr_seqs;

	# count how many times we see this feature
	my $feature_count = 0;
 
	open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

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
	
		$feature_count++;
	}
	close($in);

	# possible that we will get here without seeing any of the desire feature
	# (depending on what input GFF is used), so skip if $feature_count == 0
	next if ($feature_count == 0);


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

#		print "$desired_gff_feature\t$original_junction_bp\t$original_non_junction_bp\t$remaining_junction_bp\t$remaining_non_junction_bp\n";

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

	my $ratio = sprintf("%.4f", $percent_overlapping_junction / $percent_overlapping_non_junction);
	$percent_overlapping_junction     = sprintf("%.2f", $percent_overlapping_junction);
	$percent_overlapping_non_junction = sprintf("%.2f", $percent_overlapping_non_junction);
	
	$output_text .= "$feature_overlapping_junction/$original_junction_bp bp (%$percent_overlapping_junction)\t";
	$output_text .= "$feature_overlapping_non_junction/$original_non_junction_bp bp (%$percent_overlapping_non_junction)";

	$results_by_ratio{$output_text} = $ratio;

	# if we are shuffling, want to count whether this ratio beats real ratio in unshuffled data
	next unless ($i > 0);
	
	if($ratio > $final_ratios{$desired_gff_feature}{real_result}){
		$final_ratios{$desired_gff_feature}{over}++;
#			print "Shuffle $i) OVER REAL RATIO: $ratio $output_text\n" if ($verbose);
	} elsif($ratio < $final_ratios{$desired_gff_feature}{real_result}){
		$final_ratios{$desired_gff_feature}{under}++;
#			print "Shuffle $i) UNDER REAL RATIO: $ratio $output_text\n" if ($verbose);		
	} else{
		$final_ratios{$desired_gff_feature}{same}++;
#			print "Shuffle $i) SAME AS REAL RATIO: $ratio $output_text\n" if ($verbose);				
	}
}

# this part only happens in main run with real (unshuffled) data
if ($i == 0){
	print "Ratio\tGFF_feature\tInside_junction_overlap\tOutside_junction_overlap\n";
	foreach my $result (sort {$results_by_ratio{$b} <=> $results_by_ratio{$a}} keys %results_by_ratio){
		print "$results_by_ratio{$result}\t$result\n";

		# want to store result of final ratio in a hash to see if it is beaten in shuffling
		my @results = split(/\t/, $result);
		$final_ratios{$results[0]}{real_result} = $results_by_ratio{$result};
	}
} else {
	foreach my $feature (keys %final_ratios){
		my $real_ratio                              = $final_ratios{$feature}{real_result};
		my $shuffled_results_are_over_real_ratio    = $final_ratios{$feature}{over};
		my $shuffled_results_are_under_real_ratio   = $final_ratios{$feature}{under};
		my $shuffled_results_are_same_as_real_ratio = $final_ratios{$feature}{same};

		$shuffled_results_are_over_real_ratio    = 0 if (not defined $shuffled_results_are_over_real_ratio);
		$shuffled_results_are_under_real_ratio   = 0 if (not defined $shuffled_results_are_under_real_ratio);
		$shuffled_results_are_same_as_real_ratio = 0 if (not defined $shuffled_results_are_same_as_real_ratio);

		if($verbose and $i < $shuffles){
			print "Shuffle $i/$shuffles\t$feature\t$real_ratio\t";
			print "$shuffled_results_are_over_real_ratio\t";
			print "$shuffled_results_are_same_as_real_ratio\t";
			print "$shuffled_results_are_under_real_ratio\n";
		} elsif ($i == $shuffles) {
			print "Shuffle $i/$shuffles\t$feature\t$real_ratio\t";
			print "$shuffled_results_are_over_real_ratio\t";
			print "$shuffled_results_are_same_as_real_ratio\t";
			print "$shuffled_results_are_under_real_ratio\n";			
		}

	}
}



exit;


####################################
#
# S U B R O U T I N E S
#
####################################

# create fake chromosome sequences
sub initalize_virtual_chromosome_sequences{
	foreach my $chr (qw(Chr1 Chr2 Chr3 Chr4 Chr5)){
		$chr_seqs{$chr} = '-' x $chr_sizes{$chr};
	
		# have to do something different for Chromosome 4 as we are only interested in the 
		# tail swap region. So mask out first part of chromosome with a different character
		if($chr eq 'Chr4'){
			substr($chr_seqs{$chr}, 0, $pre_tailswap_length) = ("x" x $pre_tailswap_length);		
		}
	}
}


