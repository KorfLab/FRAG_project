#!/usr/bin/perl
#
# determine_background_levels.pl
#
# A script to work out what fraction of the (Arabidopsis) genome is occupied by different
# genomic features
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use Getopt::Long;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($feature_gff, $help) = (undef, undef);

my $usage = "$0 
Mandatory arguments:
--feature_gff <GFF file of target feature>
--help 
	
";

GetOptions (
	"feature_gff=s"    => \$feature_gff,
	"help"             => \$help,
);


die $usage if ($help);
die $usage if (not defined $feature_gff);


####################################
# Set up chromosome data
####################################

# Need a hash to store sizes of chromosomes 
my %chr_sizes;
$chr_sizes{'Chr1'} = 30427671;
$chr_sizes{'Chr2'} = 19698289;
$chr_sizes{'Chr3'} = 23459830;
$chr_sizes{'Chr4'} = 18585056;
$chr_sizes{'Chr5'} = 26975502;

my $genome_size = $chr_sizes{'Chr1'} + $chr_sizes{'Chr2'} + $chr_sizes{'Chr3'} + $chr_sizes{'Chr4'} + $chr_sizes{'Chr5'};

# will end up representing each chromosome as a string of dashes
my %chr_seqs;

# we will need to generate a fake 
# string corresponding to the length of each chromosome
initalize_virtual_chromosome_sequences();

my @gff_features;
push (@gff_features, qw(CDS DNAseI_hypersensitive_site DNA_replication_origin exon));
push (@gff_features, qw(five_prime_UTR gene mRNA miRNA ncRNA));
push (@gff_features, qw(protein pseudogene pseudogenic_exon pseudogenic_transcript)); 
push (@gff_features, qw(satellite snoRNA tRNA three_prime_UTR transposable_element));
push (@gff_features, qw(transposable_element_gene transposon_fragment));

# this is a bit of a hack, will have 9 different features for open_chromatin_state
# even though they all share the same GFF feature name
push (@gff_features, qw(open_chromatin_state_1 open_chromatin_state_2));
push (@gff_features, qw(open_chromatin_state_3 open_chromatin_state_4));
push (@gff_features, qw(open_chromatin_state_5 open_chromatin_state_6));
push (@gff_features, qw(open_chromatin_state_7 open_chromatin_state_8));
push (@gff_features, qw(open_chromatin_state_9));


foreach my $gff_feature (@gff_features){
	
	# create copies of virtual chromosome sequences
	my %tmp_chr_seqs = %chr_seqs;

	# count how many times we see this feature
	my $feature_count = 0;
 
	open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

	LINE: while(my $line = <$in>){
		chomp($line);
		# skip GFF header lines
		next if ($line =~ m/^#/);
		
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

		# skip non nuclear chromosomes
		next LINE if ($chr eq 'ChrC');
		next LINE if ($chr eq 'ChrM');
		
		# only want to look at one feature at a time
		if ($gff_feature =~ m/open_chromatin_state/){
			my ($number) = $gff_feature =~ m/_(\d)/;
			next LINE unless ($feature =~ m/open_chromatin_state/ and $comment =~ m/state$number/);
		} else{
			next LINE unless $feature eq $gff_feature;
		}	

		# mask where feature occurs in virtual chromosome sequence
		my $length = $e - $s + 1;
		substr($tmp_chr_seqs{$chr}, $s, $length) = ("o" x $length);
	}
	close($in);

	my $feature_bp = 0;
	
	foreach my $chr (keys %tmp_chr_seqs){
		$feature_bp += $tmp_chr_seqs{$chr} =~ tr/o/o/;	
	}
	my $percent_feature = sprintf("%.2f", ($feature_bp / $genome_size) * 100);
	
	print "$gff_feature\t$genome_size\t$feature_bp\t$percent_feature\n";
}




exit;


####################################
#
# S U B R O U T I N E S
#
####################################

# create fake chromosome sequences as strings of '-'
sub initalize_virtual_chromosome_sequences{
	foreach my $chr (qw(Chr1 Chr2 Chr3 Chr4 Chr5)){
		$chr_seqs{$chr} = '-' x $chr_sizes{$chr};	
	}
}

