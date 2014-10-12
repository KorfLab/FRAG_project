#!/usr/bin/perl
#
# check_feature_orientation.pl
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

my ($breakpoint_gff, $feature_gff, $help, $verbose, $target_feature) = (undef, undef, undef, undef, undef);

my $usage = "$0 
Mandatory arguments:
--breakpoint_gff <gff file of breakpoint junction coordinates> 
--feature_gff <master GFF file of all TAIR10 features>
--target_feature <name of GFF feature to assess>

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
);


die $usage if ($help);
die $usage if (not defined $breakpoint_gff);
die $usage if (not defined $feature_gff);
die $usage if (not defined $target_feature);


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


# now loop over breakpoints from our FRAG GFF file
my %breakpoint_details;
my $total = 0;

open ($in, "<", $breakpoint_gff) or die "Can't read $breakpoint_gff\n";

while(my $line = <$in>){
	my ($chr, undef, $feature, $s, $e, undef, undef, undef, $comment) = split(/\t/, $line);	

	# skip comments
	next if ($line =~ m/^#/);
	
	# only want breakpoint features
	next unless ($feature eq 'chromosome_breakpoint' and $comment !~ m/telomeric end/);
	
	# default length to explore either side of breakpoint, will increase this 
	# if we don't find features
	my $length = 25;

	# when we know orientation of features that occur either side of breakpoint
	# we will set $breakpoint_status to 1
	my $breakpoint_status = 0;

	print "$line" if ($verbose);
	$total++;

	
	while($breakpoint_status == 0){

		# extract region either side of breakpoint
		my $region;

		
		my $extraction_start;
		my $extraction_length;
		
		# extraction start can't be less than 1st bp of chromosome
		if(($s - $length) < 1){
			$extraction_start = 1;
		} else{
			$extraction_start = $s - $length;
		}
		
		# likewise can't go beyond length of chromosome
		if ($s + 1 + $length > $chr_sizes{$chr}){
			$extraction_length = $chr_sizes{$chr} - $s;
		} else{
			$extraction_length = $length;
		}
		
		# 5' region		
		$region = substr($chr_seqs{$chr}, $extraction_start, $extraction_length);
		# add symbol for breakpoint
		$region .= " | ";
		# 3' region
		$region .= substr($chr_seqs{$chr}, $s + 1, $extraction_length);
	
		# first look for forward strand orientation either side of breakpoint
		if ($region =~ m/> \| >/){
			$breakpoint_details{'>>>|>>>'}++;
			$breakpoint_status = 1;
			print "1) $region\n\n" if ($verbose); 
		} elsif ($region =~ m/< \| </){
			# now look for reverse strand orientation either side of breakpoint
			$breakpoint_details{'<<<|<<<'}++;
			$breakpoint_status = 1;
			print "2) $region\n\n" if ($verbose);
		} elsif ($region =~ m/< \| >/){
			# now look for reverse + forward strand orientation either side of breakpoint
			$breakpoint_details{'<<<|>>>'}++;
			$breakpoint_status = 1;
			print "3) $region\n\n" if ($verbose);
		} elsif ($region =~ m/> \| </){
			# now look for forward + reverse strand orientation either side of breakpoint (very rare!)
			$breakpoint_details{'>>>|<<<'}++;
			$breakpoint_status = 1;
			print "4) $region\n\n" if ($verbose);
		} elsif ($region =~ m/<\-+ \| \-+</ or $region =~ m/<\-+ \| </){
			# now look for reverse strand orientation either side of breakpoint but with non-feature sequence
			$breakpoint_details{'<<<---|---<<<'}++;
			$breakpoint_status = 1;
			print "5) $region\n\n" if ($verbose);		
		} elsif ($region =~ m/>\-+ \| \-+>/){
			# now look for forward strand orientation either side of breakpoint but with non-feature sequence
			$breakpoint_details{'>>>---|--->>>'}++;
			$breakpoint_status = 1;
			print "6) $region\n\n" if ($verbose);		
		} elsif ($region =~ m/<\-+ \| \-+>/){
			# now look for reverse/forward strand orientation either side of breakpoint but with non-feature sequence
			$breakpoint_details{'<<<---|--->>>'}++;
			$breakpoint_status = 1;
			print "7) $region\n\n" if ($verbose);		
		} elsif ($region =~ m/>\-+ \| \-+</){
			# now look for forward/reverse strand orientation either side of breakpoint but with non-feature sequence
			$breakpoint_details{'>>>---|---<<<'}++;
			$breakpoint_status = 1;
			print "8) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/ooo \| ooo/){
			# now look for non-stranded orientation either side of breakpoint
			$breakpoint_details{'ooo|ooo'}++;
			$breakpoint_status = 1;
			print "9) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/\-+ \| \-+/){
			# non-target-features only either side
			$breakpoint_details{'---|---'}++;
			$breakpoint_status = 1;
			print "10) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/o+ \| \-+/){
			# non-strand & non-target-feature
			$breakpoint_details{'ooo|---'}++;
			$breakpoint_status = 1;
			print "11) $region\n\n" if ($verbose);		
		} elsif ($stranded == 0 and $region =~ m/\-+ \| o+/){
			# non-target-feature & non-stranded feature
			$breakpoint_details{'---|ooo'}++;
			$breakpoint_status = 1;
			print "12) $region\n\n" if ($verbose);	
		} else{

			# deal with edge cases this is a bit of a fudge.
			# $length shouldn't be more than chr length
			
			if ($length >= $chr_sizes{$chr}){
				$breakpoint_details{'???|???'}++;
				$breakpoint_status = 1;
				print "13) EDGE CASE\n\n" if ($verbose);
				
			} 
			
			# if we haven't seen neighboring features double the size of the region
			$length *= 2;
#			print "Increasing length to $length\n";
#			print "$region\n\n";
		}
	}
}

close($in);

foreach my $key (keys %breakpoint_details){
	my $percent = sprintf("%.2f", $breakpoint_details{$key} / $total * 100);
	print "$key\t$breakpoint_details{$key}\t%$percent\n";
}
exit;

