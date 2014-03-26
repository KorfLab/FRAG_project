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

my ($junction_gff, $feature_gff, $bp, $help, $verbose) = (undef, undef, 10000, undef, undef);

my $usage = "$0 
Mandatory arguments:
--junction_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>

Optional arguments:
--bp <how many bp to extract from inside/outside region (default = $bp)>
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"junction_gff=s" => \$junction_gff,
	"feature_gff=s"  => \$feature_gff,
	"bp=i"           => \$bp,
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
#$chr_sizes{'Chr1'} = 28315915; # a fudge to deal with tail swap, coordinate from Han
$chr_sizes{'Chr2'} = 19698289;
$chr_sizes{'Chr3'} = 23459830;
$chr_sizes{'Chr4'} = 18585056;
$chr_sizes{'Chr5'} = 26975502;

# where does the tail swap begin on Chr4?
#my $pre_tailswap_length = 16541500; # coordinate from Han

# will end up representing each chromosome as a string of dashes
# do this for all data and just for junction regions
my %chr_seqs;
my %chr_junction_seqs;

# for the main run and each shuffle, we will need to generate a fake 
# string corresponding to the length of each chromosome
initalize_virtual_chromosome_sequences();

# now read junction data
#read_junction_data();


##########################################
# Main loop over each possible feature
##########################################


my $desired_gff_feature = "gene";


open (my $in, "<", $feature_gff) or die "Can't read $feature_gff\n";

while(my $line = <$in>){
	my ($chr, undef, $feature, $s, $e, undef, $strand, undef, $comment) = split(/\t/, $line);	

	# only want protein coding genes
	next unless ($feature eq $desired_gff_feature);
	next unless ($comment =~ m/protein_coding_gene/);

	# skip plastid chromosomes
	next if ($chr eq 'ChrM' or $chr eq 'ChrC');	

	# temporarly skip if not chr1 or chr4?
#	next unless (($chr eq 'Chr1') or ($chr eq 'Chr4'));

	# skip tailswap regions of Chr1 and Chr4
#	next if ($chr eq 'Chr1' and $s > $chr_sizes{'Chr1'});
#	next if ($chr eq 'Chr4' and $s < $pre_tailswap_length);

	# mask where feature occurs in virtual chromosome sequence
	my $length = $e - $s + 1;
	my $replacement_string;
	$replacement_string = ">" if ($strand eq '+');
	$replacement_string = "<" if ($strand eq '-');
	substr($chr_seqs{$chr}, $s, $length) = ($replacement_string x $length);
	my $tmp_seq = substr($chr_seqs{$chr}, $s-100, $length+200);
	
}
close($in);


my %orientation;
my $total = 0;
foreach my $chr (sort keys %chr_seqs){
	print "$chr\n";
	my $seq = $chr_seqs{$chr};

	# shorten  sequences
	$seq =~ s/\>{2,}/\>/g;
	$seq =~ s/\<{2,}/\</g;
	$seq =~ s/\-//g;
	print "$seq\n";

	for (my $i = 0; $i < length($seq) - 1; $i++){
		my $two_genes = substr($seq, $i, 2);
		$orientation{$two_genes}++;
		$total++;
		#		print "$two_genes ";
	}
	print "\n";
}

foreach my $key (keys %orientation){
	my $percent = sprintf("%.2f", $orientation{$key} / $total * 100);
	print "$key\t$orientation{$key}\t%$percent\n";
}


# now loop over junctions
my %junction_details;
$total = 0;
open ($in, "<", $junction_gff) or die "Can't read $junction_gff\n";

while(my $line = <$in>){
	my ($chr, undef, undef, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

	# coordinate that we use depends on whether this was the left or right edge
	my $coord = $s;
	$coord = $s if ($comment =~ m/edge=R/);
	$coord = $e if ($comment =~ m/edge=L/);

	my $junction_status = 0;
	my $length = 25;
	print "$line";
	$total++;
	while($junction_status == 0){
		my $region = substr($chr_seqs{$chr}, $coord-$length, $length);
		$region .= " | ";
		$region .= substr($chr_seqs{$chr}, $coord + 1, $length);
	
		if ($region =~ m/> \| >/){
			$junction_details{'>>>|>>>'}++;
			$junction_status = 1;
			print "1) $region\n\n";
		} elsif ($region =~ m/< \| </){
			$junction_details{'<<<|<<<'}++;
			$junction_status = 1;
			print "2) $region\n\n";
		} elsif ($region =~ m/< \| >/){
			$junction_details{'<<<|>>>'}++;
			$junction_status = 1;
			print "3) $region\n\n";
		} elsif ($region =~ m/> \| </){
			$junction_details{'>>>|<<<'}++;
			$junction_status = 1;
			print "4) $region\n\n";
		} elsif ($region =~ m/<\-+ \| \-+</ or $region =~ m/<\-+ \| </){
			$junction_details{'<<<---|---<<<'}++;
			$junction_status = 1;
			print "5) $region\n\n";		
		} elsif ($region =~ m/>\-+ \| \-+>/){
			$junction_details{'>>>---|--->>>'}++;
			$junction_status = 1;
			print "6) $region\n\n";		
		} elsif ($region =~ m/<\-+ \| \-+>/){
			$junction_details{'<<<---|--->>>'}++;
			$junction_status = 1;
			print "7) $region\n\n";		
		} elsif ($region =~ m/>\-+ \| \-+</){
			$junction_details{'>>>---|---<<<'}++;
			$junction_status = 1;
			print "8) $region\n\n";		
		} else{
			$length *= 2;
#			print "Increasing length to $length\n$region\n\n";
		}

	}
}

close($in);

foreach my $key (keys %junction_details){
	my $percent = sprintf("%.2f", $junction_details{$key} / $total * 100);
	print "$key\t$junction_details{$key}\t%$percent\n";
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
		$chr_junction_seqs{$chr} = '-' x $chr_sizes{$chr};
	
		# have to do something different for Chromosome 4 as we are only interested in the 
		# tail swap region. So mask out first part of chromosome with a different character
#		if($chr eq 'Chr4'){
#			substr($chr_seqs{$chr}, 0, $pre_tailswap_length) = ("x" x $pre_tailswap_length);		
#		}
	}
}


###############################################################
# read junction coordinates and represent in virtual sequences
###############################################################

sub read_junction_data{
	open (my $in, "<", $junction_gff) or die "Can't read $junction_gff\n";

	while(my $line = <$in>){
		my ($chr, undef, undef, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

		# coordinate that we use depends on whether this was the left or right edge
		my $coord = $s;
		$coord = $s if ($comment =~ m/edge=R/);
		$coord = $e if ($comment =~ m/edge=L/);
	
		# now define a range around this coordinate based on value of $bp
		my ($min, $max) = ($coord - $bp/2, $coord + $bp/2);
	
		# mask where junctions are in chromosome
		substr($chr_junction_seqs{$chr}, $min, $bp) = ("J" x $bp);
	}

	close($in);

}


