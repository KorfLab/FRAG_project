#!/usr/bin/perl
#
# determine_background_feature_enrichment.pl
#
# A script to work out what fraction of the (Arabidopsis) genome is occupied by different
# genomic features, and look at the levels of these features when the genome is 
# separated into duplicated blocks, triplicated blocks, and everything else.
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use FRAG;

###################################################### 
#              Command-line options                  #
###################################################### 

my ($feature_gff, $breakpoint_gff, $help) = (undef, undef, undef);

my $usage = "$0 
Mandatory arguments:
--feature_gff <GFF file of all target features>
--breakpoint_gff <GFF file of breakpoint junction coordinates> 
--help <this help>
	
";

GetOptions (
	"feature_gff=s"    => \$feature_gff,
	"breakpoint_gff=s" => \$breakpoint_gff,
	"help"             => \$help,
);


die $usage if ($help);
die $usage if (not defined $feature_gff);
die $usage if (not defined $breakpoint_gff);


####################################
# Set up chromosome data
####################################

# Use a hash to store sizes of chromosomes 
# keys will be 'chr1', 'chr2' etc.
# want the shorter 'tailswap' length for Chr1
my %chr_sizes = FRAG::get_chromosome_lengths(); 


# where does the tail swap begin on Chr4?
my $pre_tailswap_length = FRAG::get_chr4_pre_tailswap_length();


# get list of all target GFF features we are interested in
my @target_gff_features = sort (FRAG::get_list_of_GFF_features($feature_gff));



# want to perform analysis at four levels: whole genome, 1x, 2x, and 3x
# for the latter 3 categories we'll only look at the Chr1/Chr4 regions involved
# in genome shattering

foreach my $type (qw(genome 1x 2x 3x)){
	

	foreach my $target_feature (@target_gff_features){
		# will create virtual  chromosome sequences that are first masked for 
		# the parent 1x, 2x, 3x blocks and then masked for specific genome features
		my %chr_seqs;
	
		print "$type\t$target_feature\t";

		my $bp_before = 0;
		
		# we will mask sequences slightly differently for 1x-3x cases
		if ($type eq 'genome'){
			%chr_seqs = FRAG::initalize_virtual_chromosome_sequences();
			# bp_before in this case is just equal to the genome size
			$bp_before = $chr_sizes{'Chr1'} + $chr_sizes{'Chr2'} + $chr_sizes{'Chr3'} + $chr_sizes{'Chr4'} + $chr_sizes{'Chr5'};
		} else{
			%chr_seqs = FRAG::initalize_virtual_chromosome_sequences('tailswap');
			$bp_before = read_block_data($type, \%chr_seqs);
		}
		
		open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

		LINE: while(my $line = <$in>){
			chomp($line);

			# skip GFF header lines
			next if ($line =~ m/^#/);
		
			my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

			# skip non-nuclear chromosomes
			next LINE if ($chr eq 'ChrC');
			next LINE if ($chr eq 'ChrM');
		
			#  want to separate out chromatin state data (if it exists)
			if($feature eq 'open_chromatin_state'){
				my ($state) = $comment =~ m/Note=\"state(\d+)\"/;
				$feature = "open_chromatin_state_${state}";
			}

			# skip lines that are not the same as our current target feature
			next unless ($feature eq $target_feature);
			

			# for the entire genome, we want to mask every feature but for 1x-3x types
			# only want to focus on regions of Chr1 and Chr4
			if ($type =~ m/\dx/){
		
				# For most calculations only want to mask features that occur on chr1 (before tailswap) 
				# and on chr4 (after tailswap) 
				next LINE unless (($chr eq 'Chr1') or ($chr eq 'Chr4'));

				# skip tailswap regions of Chr1 and Chr4
				next LINE if ($chr eq 'Chr1' and $s > $chr_sizes{'Chr1'});
				next LINE if ($chr eq 'Chr4' and $s < $pre_tailswap_length);
			}

			# can now mask where feature occurs in virtual chromosome sequences
			my $length = $e - $s + 1;
			substr($chr_seqs{$chr}, $s, $length) = ("o" x $length);
		}
		close($in);

		# how many block or feature bp remain masked
		my $bp_after = 0;
		my $genome_wide_bp = 0;
		my $percent;
		
		# will do things a bit differently for the genome-level stats
		if ($type eq 'genome'){
			foreach my $chr (keys %chr_seqs){
				$genome_wide_bp += $chr_seqs{$chr}  =~ tr/o/o/;	
			}
			$percent = sprintf("%.2f", ($genome_wide_bp / $bp_before) * 100);
			print "\t$genome_wide_bp\t$percent%\n";

		} else{
			foreach my $chr (qw(Chr1 Chr4)){
				$bp_after += $chr_seqs{$chr} =~ tr/B/B/;	
			}
			my $bp_occupied_by_feature = $bp_before - $bp_after;
			$percent = sprintf("%.2f", ($bp_occupied_by_feature / $bp_before) * 100);
			print "\t$bp_occupied_by_feature\t$percent%\n";
		}
	}
}

exit;


####################################
#
# S U B R O U T I N E S
#
####################################


sub read_block_data{

	my ($type, $hash_ref) = @_;
	
	# how many bp inside blocks of type $type 	
	my $bp_before = 0;
	
	open (my $in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

	while(my $line = <$in>){
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

		# skip comments
		next if ($line =~ m/^#/);

		# only interested in 1x, 2x or 3x blocks in FRAG GFF file
		next unless ($feature eq 'copy_number_gain' or $comment =~ m/no copy number gain detected/);

		
		my $copy_number;
		if ($comment =~ m/no copy number gain detected/){
			$copy_number = "1x";	
		} elsif ($comment =~ m/duplicated block/){
			$copy_number = "2x";	
		} else{
			$copy_number = "3x";			
		}

		# only mask duplicated/triplicated blocks in certain hash keys
		next unless ($copy_number eq $type);

		my $block_length = $e - $s + 1;

		substr(${$hash_ref}{$chr}, $s, $block_length) = ("B" x $block_length);
		
		# keep running tally of total block length
		$bp_before += $block_length;
	}
	close($in);

	return($bp_before);
}