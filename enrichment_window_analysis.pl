#!/usr/bin/perl
#
# enrichment_window_analysis.pl
#
# A script to look at orientation of specific features iN TAIR10 GFF file with respect
# to mapped breakpoints in FRAG data
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use List::Util qw(shuffle);
use Getopt::Long;
use FRAG;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($breakpoint_gff, $feature_gff, $help, $verbose, $target_feature, $bp) = (undef, undef, undef, undef, undef, undef);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>
--target_feature <name of GFF feature to assess>
--bp <how many bp the window around breakpoint should be>

Optional arguments:
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"breakpoint_gff=s" => \$breakpoint_gff,
	"feature_gff=s"    => \$feature_gff,
	"help"             => \$help,
	"target_feature=s" => \$target_feature,
	"verbose"          => \$verbose,
	"bp=i"             => \$bp,
);


die $usage if ($help);
die $usage if (not defined $breakpoint_gff);
die $usage if (not defined $feature_gff);
die $usage if (not defined $target_feature);
die $usage if (not defined $bp);

####################################
# Set up chromosome data
####################################

# Use a hash to store sizes of chromosomes 
# keys will be 'chr1', 'chr2' etc.
my %chr_sizes = FRAG::get_chromosome_lengths; 

# for the main run and each shuffle, we will need to generate a fake 
# string corresponding to the length of each chromosome
# will end up representing each chromosome as a string of dashes
my %chr_seqs = FRAG::initalize_virtual_chromosome_sequences();

# get details of breakpoints coordinates
my %breakpoints = FRAG::read_breakpoint_data($breakpoint_gff);



# need to know background percent of target features across all of genome
my $background_ratio;
$background_ratio = 50.47 if ($target_feature eq 'gene');
$background_ratio = 4.17  if ($target_feature eq 'DNA_replication_origin');


# what is the copy number
my $copy_number;
$copy_number = "2x" if ($breakpoint_gff =~ m/2x/);
$copy_number = "3x" if ($breakpoint_gff =~ m/3x/);

##########################################
# Main loop over all genes
##########################################

# want to know whether we are looking at a stranded feature or not
# so far, all features apart from replication origins are stranded
# so set default to be 1
my $stranded = 1;

open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

while(my $line = <$in>){
	my ($chr, undef, $feature, $s, $e, undef, $strand, undef, $comment) = split(/\t/, $line);	

	# skip comment lines if present
	next if ($line =~ m/^#/);
	
	# only want to look at the one desired feature
	next unless ($feature eq $target_feature);

	# skip plastid chromosomes
	next if ($chr eq 'ChrM' or $chr eq 'ChrC');	

	# mask where feature occurs in virtual chromosome sequence
	my $length = $e - $s + 1;
	my $replacement_string = "o";
	substr($chr_seqs{$chr}, $s, $length) = ($replacement_string x $length);
	
}
close($in);


my $breakpoint_counter = 0;


BREAKPOINT: foreach my $breakpoint (sort {$a <=> $b} keys %breakpoints){
	$breakpoint_counter++;
	my $b_pos = $breakpoints{$breakpoint}{pos};
	my $b_chr = $breakpoints{$breakpoint}{chr};

	# define a region around this breakpoint based on value of $bp
	my ($min, $max);

	# slight difference depending on whether $bp is even or odd
	if ($bp % 2 == 0){
		($min, $max) =($b_pos - $bp/2, $b_pos + $bp/2);
	} else{
		($min, $max) =($b_pos - int($bp/2), $b_pos + int($bp/2));
	}

	# extraction start can't be less than 1st bp of chromosome
	($min = 1) if ($min < 1);

	# likewise can't go beyond length of chromosome
	($max = $chr_sizes{$b_chr}) if ($max > $chr_sizes{$b_chr});

	# extract this region 
	my $region = substr($chr_seqs{$b_chr}, $min, $bp);

	my $region_for_display;
	$region_for_display  = substr($chr_seqs{$b_chr}, $min, int($bp/2));
	$region_for_display .= " | ";
	$region_for_display .= substr($chr_seqs{$b_chr}, $b_pos, int($bp/2));

	my ($feature_bp) = $region =~ tr/o/o/;
	my $percent = sprintf("%.2f", $feature_bp / $bp * 100);
	my $enrichment_ratio = sprintf("%.2f", $percent / $background_ratio);
	my $id = "${target_feature}_${copy_number}_${bp}_$breakpoint_counter";
	print "$id\t$feature_bp\t$percent%\t$background_ratio%\t$enrichment_ratio\n";
}

exit;

