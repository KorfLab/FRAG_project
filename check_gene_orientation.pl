#!/usr/bin/perl
#
# check_gene_orientation.pl
#
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

my ($breakpoint_gff, $feature_gff, $bp, $help, $verbose, $target_feature) = (undef, undef, 10000, undef, undef, undef);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>
--target_feature <name of GFF feature to assess>

Optional arguments:
--bp <how many bp to extract from inside/outside region (default = $bp)>
--verbose - turn on extra output
--help 
	
";

GetOptions (
	"breakpoint_gff=s" => \$breakpoint_gff,
	"feature_gff=s"    => \$feature_gff,
	"bp=i"             => \$bp,
	"help"             => \$help,
	"target_feature=s" => \$target_feature,
	"verbose"          => \$verbose,
);


die $usage if ($help);
die $usage if (not defined $breakpoint_gff);
die $usage if (not defined $feature_gff);
die $usage if (not defined $target_feature);


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
	
	# only want protein coding genes
#	next unless ($feature eq 'gene' and $comment =~ m/protein_coding_gene/);
	next unless ($feature eq $target_feature);

	# skip plastid chromosomes
	next if ($chr eq 'ChrM' or $chr eq 'ChrC');	

	# mask where feature occurs in virtual chromosome sequence
	# type of masking will depend on whether feature is stranded or not
	my $length = $e - $s + 1;
	my $replacement_string;
	if ($strand eq '+'){
		$replacement_string = ">";
	} elsif ($strand eq '-'){
		$replacement_string = "<";	
	} elsif ($strand eq '.'){
		$replacement_string = "o";
		$stranded = 0;
	} else {
		die "Can't find suitable strand character\n";
	}

	substr($chr_seqs{$chr}, $s, $length) = ($replacement_string x $length);
	
}
close($in);

# some extra details if we are looking at gene features
if($target_feature eq 'gene'){
	my $coding_bases = 0;
	my $total_bases = 0;
	foreach my $chr (keys %chr_seqs){

		my ($forward) = $chr_seqs{$chr} =~ tr/>/>/;
		my ($reverse) = $chr_seqs{$chr} =~ tr/</</;

		$coding_bases += $forward;
		$coding_bases += $reverse;
		$total_bases += length($chr_seqs{$chr});
	}

	my $percent = sprintf("%.2f", $coding_bases/$total_bases * 100);
	print "$coding_bases coding bp in genome out of $total_bases bp ($percent%)\n";
}

# now loop over breakpoints
my %breakpoint_details;
my $total = 0;
open ($in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

while(my $line = <$in>){
	my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

	# skip comments
	next if ($line =~ m/^#/);
	
	# only want breakpoint features
	next unless ($feature eq 'chromosome_breakpoint' and $comment !~ m/telomeric end/);

	my $junction_status = 0;
	my $length = 25;
	print "$line" if ($verbose);
	$total++;
	while($junction_status == 0){
		my $region = substr($chr_seqs{$chr}, $s-$length, $length);
		$region .= " | ";
		$region .= substr($chr_seqs{$chr}, $s + 1, $length);
	
		if ($region =~ m/> \| >/){
			$breakpoint_details{'>>>|>>>'}++;
			$junction_status = 1;
			print "1) $region\n\n" if ($verbose); 
		} elsif ($region =~ m/< \| </){
			$breakpoint_details{'<<<|<<<'}++;
			$junction_status = 1;
			print "2) $region\n\n" if ($verbose);
		} elsif ($region =~ m/< \| >/){
			$breakpoint_details{'<<<|>>>'}++;
			$junction_status = 1;
			print "3) $region\n\n" if ($verbose);
		} elsif ($region =~ m/> \| </){
			$breakpoint_details{'>>>|<<<'}++;
			$junction_status = 1;
			print "4) $region\n\n" if ($verbose);
		} elsif ($region =~ m/<\-+ \| \-+</ or $region =~ m/<\-+ \| </){
			$breakpoint_details{'<<<---|---<<<'}++;
			$junction_status = 1;
			print "5) $region\n\n" if ($verbose);		
		} elsif ($region =~ m/>\-+ \| \-+>/){
			$breakpoint_details{'>>>---|--->>>'}++;
			$junction_status = 1;
			print "6) $region\n\n" if ($verbose);		
		} elsif ($region =~ m/<\-+ \| \-+>/){
			$breakpoint_details{'<<<---|--->>>'}++;
			$junction_status = 1;
			print "7) $region\n\n" if ($verbose);		
		} elsif ($region =~ m/>\-+ \| \-+</){
			$breakpoint_details{'>>>---|---<<<'}++;
			$junction_status = 1;
			print "8) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/ooo \| ooo/){
			$breakpoint_details{'ooo|ooo'}++;
			$junction_status = 1;
			print "9) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/\-+ \| \-+/){
			$breakpoint_details{'---|---'}++;
			$junction_status = 1;
			print "10) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/o+ \| \-+/){
			$breakpoint_details{'ooo|---'}++;
			$junction_status = 1;
			print "11) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/\-+ \| o+/){
			$breakpoint_details{'---|ooo'}++;
			$junction_status = 1;
			print "12) $region\n\n" if ($verbose);	
		} else{
			$length *= 2;
#			print "Increasing length to $length\n$region\n\n";
		}
	}
}

close($in);

foreach my $key (keys %breakpoint_details){
	my $percent = sprintf("%.2f", $breakpoint_details{$key} / $total * 100);
	print "$key\t$breakpoint_details{$key}\t%$percent\n";
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

