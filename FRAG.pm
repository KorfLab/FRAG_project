#!/usr/bin/perl
#
# FRAG.pm
#
# Perl module with common utilities needed by FRAG project code
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.
use strict;
use warnings FATAL => 'all';


package FRAG;

# will store chromosome lengths in hash that will potentially be used by other subroutines
my %chr_sizes;


# will store details of the location of all breakpoints in another hash	
my %breakpoints;


# where does the tail swap begin on Chr4?
# this is needed for only analyzing first 4/5ths (ish) of Chr1 and just tailswap region
# of Chr4
my $pre_tailswap_length = 16541500; # coordinate from Han

####################################
# Set up chromosome data
####################################

sub get_chromosome_lengths{

	my ($tailswap) = @_;
	 
	# Need a hash to store sizes of chromosomes 
	$chr_sizes{'Chr1'} = 30427671;
	$chr_sizes{'Chr2'} = 19698289;
	$chr_sizes{'Chr3'} = 23459830;
	$chr_sizes{'Chr4'} = 18585056;
	$chr_sizes{'Chr5'} = 26975502;

	# do we want tailswap coordinates for chr1 (truncates the length for analysis
	# purposes)
	if ($tailswap){
		$chr_sizes{'Chr1'} = 28315915; # coordinate from Han
	}
	return %chr_sizes;
}

####################################
# Tailswap length for Chr4
####################################

sub get_chr4_pre_tailswap_length{

	return $pre_tailswap_length;
}




####################################
# Create virtual chromosome sequences
####################################

sub initalize_virtual_chromosome_sequences{

	my ($tailswap) = @_;

	my %chr_seqs;
	foreach my $chr (qw(Chr1 Chr2 Chr3 Chr4 Chr5)){
		$chr_seqs{$chr} = '-' x $chr_sizes{$chr};
	}
	
	# have to do something different for Chromosome 4 as we are only interested in the
	# tail swap region. So mask out first part of chromosome with a different character
  	if($tailswap){
    	substr($chr_seqs{'Chr4'}, 0, $pre_tailswap_length) = ("x" x $pre_tailswap_length);
    }
	return %chr_seqs;
}




###############################################################
# read breakpoint coordinates and represent in virtual sequences
###############################################################

sub read_breakpoint_data{

	my ($breakpoint_gff) = @_;
	# only want to process some of the data in this file
	# need to extract chromosome_breakpoint features:
	
	# Chr1	t_test.pl	copy_number_gain	1	205575	.	+	.	ID=block0001;Name=01a1_01a2;Note="duplicated block"
	# Chr1	PRICE	chromosome_breakpoint	1	1	.	.	.	ID=breakpoint0001;Parent=block0001;Name=01a1;Note="telomeric end"
	# Chr1	PRICE	chromosome_breakpoint	205575	205575	.	-	.	ID=breakpoint0002;Parent=block0001;Name=01a2;Note="paired with breakpoint0033"
	# Chr1	t_test.pl	copy_number_gain	31632	87466	.	+	.	ID=block0002;Name=01b1_01b2;Note="triplicated block"
	# Chr1	PRICE	chromosome_breakpoint	31632	31632	.	-	.	ID=breakpoint0003;Parent=block0002;Name=01b1;Note="paired with breakpoint0087"

	open (my $in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

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

	return(%breakpoints);
}



####################################
# Mask breakpoint regions
####################################

sub mask_breakpoint_regions{

	my ($bp, $seqs_ref, $breakpoints_ref) = @_;
	
	my %seqs = %{$seqs_ref};
	my %breakpoints = %{$breakpoints_ref};
	
	foreach my $breakpoint (keys %breakpoints){
	
		# grab breakpoint details for each breakpoint
		my $chr = $breakpoints{$breakpoint}{chr};
		my $pos = $breakpoints{$breakpoint}{pos};

		# define a region around this breakpoint based on value of $bp
		my ($min, $max) = ($pos - $bp/2, $pos + $bp/2);

		# mask where breakpoint regions are in chromosome
		substr($seqs{$chr}, $min, $bp) = ("B" x $bp);
	}
	return(%seqs);
}


###############################################################
# read junction coordinates and represent in virtual sequences
###############################################################

sub read_block_data{
	
	my ($breakpoint_gff) = @_;
	
	my %blocks;
	
	open (my $in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

	my $block_counter = 0;
	while(my $line = <$in>){
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

		# skip comments
		next if ($line =~ m/^#/);

		# only interested in 2x or 3x blocks
		next unless ($feature eq 'copy_number_gain');

		$block_counter++;
		
		# 3 things to store for each block
		$blocks{$block_counter}{'chr'}   = $chr;
		$blocks{$block_counter}{'left'}  = $s;
		$blocks{$block_counter}{'right'} = $e;
	}

	close($in);
	return(%blocks);
}


####################################
# Get list of all features in current GFF file
####################################

sub get_list_of_GFF_features{
	my ($GFF_file) = @_;
	
	my %features;
	
	open (my $in, "<", $GFF_file) or die "Can't read $GFF_file\n";

	while(my $line = <$in>){
		# skip comments
		next if ($line =~ m/^#/);
		
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	
		

		# skip unworkable chromosome feature (only 1 per chromosome)
		next if ($feature eq 'chromosome');

		# need a fix for chromatin state data as 9 different states use same GFF feature
		# will make 9 versions of the GFF feature
		if ($feature =~ m/open_chromatin_state/){
			my ($number) = $comment =~ m/state(\d)/;
			$feature .= "_$number";
		}

		$features{$feature} = 1;	
	
	}
	close($in);

	return(keys %features);
	
}

###############################################################
# read feature data from single-feature GFF file 
###############################################################

sub read_feature_data{

	my ($feature_gff) = @_;
	
	my %features;
	
	my $feature_count = 0;
	
	open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

	while(my $line = <$in>){
		chomp($line);
		# skip GFF header lines
		next if ($line =~ m/^#/);
		
		my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

		#  skip if not chr1 or chr4?
		next unless (($chr eq 'Chr1') or ($chr eq 'Chr4'));

		# skip tailswap regions of Chr1 and Chr4
		next if ($chr eq 'Chr1' and $s > $chr_sizes{'Chr1'});
		next if ($chr eq 'Chr4' and $s < $pre_tailswap_length);

		$feature_count++;
		
		# add data to hash
		$features{$feature_count}{'chr'}   = $chr;
		$features{$feature_count}{'start'} = $s;
		$features{$feature_count}{'end'}   = $e;

	}
	close($in);
	return(%features);
}

1;