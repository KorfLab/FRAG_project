#!/usr/bin/perl
#
# nearest_feature.pl
#
# A script to find the average nearest distance of various genome features
# to genome breakpoints from our chromosome shattering lines
# will ultimately combine this script with overlap_between_two_gff_files.pl and maybe
# check_gene_orientation.pl
# 
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use List::Util qw(shuffle sum);
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




# first store details of breakpoints
my %breakpoints;

read_junction_data();





####################################
#
# Main loop
#
####################################

# main loop to loop over all features in main GFF file and calculate distance to nearest
# breakpoint for each genomic feature



# hash to track counts of each type of genomic feature
my %feature_count;

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
	
	
	# keep track of how many of each feature type we see
	$feature_count{$feature}++;
	
	# now compare current feature to all breakpoints
	foreach my $breakpoint (keys %breakpoints){
		my $b_pos = $breakpoints{$breakpoint}{pos};
		my $b_chr = $breakpoints{$breakpoint}{chr};
		
		# only want to compare genomic features on same chromosome as breakpoint
		next unless ($chr eq $b_chr);

		# work out far breakpoint position is from feature (have to compare to 
		# start & end coordinates of genomic feature)
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
			}
		}
	}
}
close($in);


# For each breakpoint, we now know the distance to the nearest genomic feature (for all
# genomic features). We now want to average distances across all breakpoints. Will
# collect all distances so as to calculate mean and standard deviation later on.
my %stats;

# Main loop is over each breakpoint
foreach my $breakpoint (sort {$a <=> $b} keys %breakpoints){
	my $b_pos = $breakpoints{$breakpoint}{pos};
	my $b_chr = $breakpoints{$breakpoint}{chr};
	print "$breakpoint $b_chr:$b_pos\n" if ($verbose);

	# now loop over all genomic features
	foreach my $feat (keys %{$breakpoints{$breakpoint}{features}}){
		my $distance = $breakpoints{$breakpoint}{features}{$feat};

		print "\t$feat\t$distance\n" if ($verbose);
		
		# push all distances onto hash for processing later
		push(@{$stats{$feat}}, $distance);
	}
}

my %mean_distance;
my %standard_deviation;

foreach my $feature (keys %stats){
	# $n will nearly always equal the number of breakpoints but not always, because
	# some features (e.g. miRNA) might not occur at all on the tailswap region of Chr4
	# so there will be no distance information between these features and breakpoints
	# in tailswap region
	my $n = @{$stats{$feature}};
	my $mean  = sprintf("%.0f",sum(@{$stats{$feature}})/$n);
	my $stdev = sprintf("%.0f",sqrt(sum(map {($_ - $mean) ** 2} @{$stats{$feature}}) / ($n-1)));
	$mean_distance{$feature}  = $mean;
	$standard_deviation{$feature} = $stdev;
}

print "Feature\tAverage_distance_to_nearest_breakpoint\tStandard_deviation\tNumber_of_features\n";
foreach my $feature (sort {$mean_distance{$a} <=> $mean_distance{$b}} keys %mean_distance){
	print "$feature\t$mean_distance{$feature}\t$standard_deviation{$feature}\t$feature_count{$feature}\n";
}



exit;


####################################
#
#  S U B R O U T I N E S
#
####################################


###############################################################
# read junction coordinates and represent in virtual sequences
###############################################################

sub read_junction_data{

	# only want to process some of the data in this file
	# need to extract chromosome_breakpoint features:
	
	# Chr1	t_test.pl	copy_number_gain	1	205575	.	+	.	ID=block0001;Name=01a1_01a2;Note="duplicated block"
	# Chr1	PRICE	chromosome_breakpoint	1	1	.	.	.	ID=breakpoint0001;Parent=block0001;Name=01a1;Note="telomeric end"
	# Chr1	PRICE	chromosome_breakpoint	205575	205575	.	-	.	ID=breakpoint0002;Parent=block0001;Name=01a2;Note="paired with breakpoint0033"
	# Chr1	t_test.pl	copy_number_gain	31632	87466	.	+	.	ID=block0002;Name=01b1_01b2;Note="triplicated block"
	# Chr1	PRICE	chromosome_breakpoint	31632	31632	.	-	.	ID=breakpoint0003;Parent=block0002;Name=01b1;Note="paired with breakpoint0087"

	open (my $in, "<", $junction_gff) or die "Can't read $junction_gff\n";

	my $breakpoint_count = 0;

	while(my $line = <$in>){
	
		# skip comments
		next if ($line =~ m/^#/);
		
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	
		
		# want chromosome breakpoints, but ignore any which are effectively the ends of
		# the chromosomes
		next unless ($feature eq 'chromosome_breakpoint' and $comment !~ m/telomeric end/);
		
		$breakpoint_count++;
	
		# store details
		$breakpoints{$breakpoint_count}{chr} = $chr;
		$breakpoints{$breakpoint_count}{pos} = $s;
	}
	close($in);
}



