FRAG_project
============

To contain code and data files for the Fragmentation project, investigating 
shattered chromosome phenotypes in _Arabidopsis thaliana_.
	
	
## Background ## 

Raw data from Han's mapping of reads back to reference genome is in the file:

	CGRs_1kb_plots_forkeith.tsv

This file contains raw and normalized counts of reads mapped to 1 Kbp bins in three lines 
(FRAG62, FRAG133, and FRAG80 - a control line). The normalized counts don't just use 
FRAG80 but are based on a concatenated file of all the F1 diploids that were sequenced. 
Han calls this the 'SuperF1diploid'.

The normalized counts come out to approximately 2.0 for duplications, and 3.0 for 
triplications. Actual values are slightly lower than this.



## Approach 1: Peak finder ##

Each row in the above TSV file corresponds to a 1 Kbp bin. Want to take all windows of 
different sized bins and perform t-tests to compare average copy number inside the bin 
with either:

a) all bins outside that window 
b) some maximum sized window either side of the target window.

Rather than treat the whole chromosome as one thing, we can instead split it into 
overlapping 'chunks' and pretty much treat each chunk as its own chromosome. Chunks should 
be big enough to contain ends of large 2x/3x blocks, or the entirety of a small duplicated 
block. For any chunk, we can take a target window in that chunk (minimum size of 2 data 
points) and run a a t-test in one of two ways:

1. Calculate mean1 from target window and mean2 from every data point outside target window 
(still within same chunk)
2. Calculate mean1 from target window and mean2 from all data points inside chunk

When window sizes are small and chunk sizes are large, method 2 is much quicker and 
produces a mean read count for the entire chunk (mean2) that is not very different from 
the mean of all data points outside the target window (mean1). As the window size gets 
larger though, the use of this 'chunk mean' gets a little more misleading.

Method 1 is slower but produces higher t-test values. Can compare results from both 
methods as follows: 

	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 200 --mean_mode 1 > run_200_mode1.tsv
	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 200 --mean_mode 2 > run_200_mode2.tsv
	
--mean_mode 1 = method 1 above, --mean_mode 2 = method 2.

These runs use a chunk size of 200 (data points) which equals 200 Kbp. This means the max 
window size will be 100 Kbp.



## Find maximum possible t-values by shuffling real data##

With so many t-tests, there will be significant results by chance, so we can first find 
out what is the maximum t-test score we would see if we shuffled the order of the bins n 
times. When doing this we record the maximum t-test value for windows of different sizes, 
within each chunk, within each chromosome. E.g. 

	Chr1:0-199:10   2.83146179714909 

This says that the largest t-test value from a window of 10 data points (in chunk 0-199 on 
chromosome 1) is 2.81. Ideally want to do 100 or more shufflings for any given chunk size.

This test makes 100 shuffles using window sizes at 2-100 (chunk size 200). The script 
updates the output file (CGRs_1kb_plots_forkeith.tsv.shuffled.2-100.100) after each 
shuffling so you can start using it for testing. 

	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 200 --mean_mode 1 --shuffle 100


## Final testing using shuffled data ##

Using above shuffled file, we can now look for significant t-test results in the real data 
but use the shuffled results as a control. We still want to use a minimum t-test score for 
standard significance (but we can potentially go higher than what is statistically 
significant at P < 0.01 just to be safe). The default 'min_t' value in the script is 5.

	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 200 --mean_mode 1 --shuffle_file CGRs_1kb_plots_forkeith.tsv.shuffled.2-100.100 > run_200s3_mode_1.tsv
 
's3' in the output file name refers to the fact that at runtime, the shuffled output file 
had undergone 3 shufflings.


## Results part 1: shuffling ##

At different chunk sizes (100, 200, 300, 500, and 1000) I have produced background shuffle 
results from 100 (or sometimes 1000) shufflings:

	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 100 --mean_mode 1 --shuffle 100 
	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 200 --mean_mode 1 --shuffle 100 
	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 200 --mean_mode 1 --shuffle 1000
	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 300 --mean_mode 1 --shuffle 100 
	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 500 --mean_mode 1 --shuffle 100 
	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 1000 --mean_mode 1 --shuffle 100 

This produces the following output files:

	CGRs_1kb_plots_forkeith.tsv.shuffled.1000.100
	CGRs_1kb_plots_forkeith.tsv.shuffled.100.100
	CGRs_1kb_plots_forkeith.tsv.shuffled.200.100
	CGRs_1kb_plots_forkeith.tsv.shuffled.200.1000
	CGRs_1kb_plots_forkeith.tsv.shuffled.300.100
	CGRs_1kb_plots_forkeith.tsv.shuffled.500.100
	
The	difference between 100 and 1000 shufflings means that higher t-values can be produced, 
which will ultimately reduce the number of 2x and 3x blocks that are called.


## Results part 2: final block calling ##

Can always do this 3 ways, using no background file (unshuffled = u), or using background 
file (shuffled = s) with either combination of --mean_mode 1 or 2. E.g. 

	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 100 --mean_mode 1  --shuffle_file CGRs_1kb_plots_forkeith.tsv.shuffled.100.100 > run_100s100_mode1.tsv
	./peak_finder.pl -file CGRs_1kb_plots_forkeith.tsv --chunk_size 100 --mean_mode 2  --shuffle_file CGRs_1kb_plots_forkeith.tsv.shuffled.100.100 > run_100s100_mode2.tsv

Final output files look like this:

	Chr     Chunk   Win:start-end   Win_size        Type    Copy_number     Start_coord     End_coord       t       Mean1   Mean2
	Chr1    0-99    33-82   50      fakeblock       3x      33000   82000   15.06   3.38    2.63
	Chr1    0-99    0-31    32      edge-R  2x      0       31000   9.93    2.54    3.23
	Chr1    50-149  50-87   38      edge-R  3x      50000   87000   19.46   3.35    2.51
	Chr1    50-149  98-147  50      fakeblock       2x      98000   147000  9.93    2.50    3.15
	Chr1    100-199 142-145 4       block   2x      142000  145000  6.90    2.00    2.53
	Chr1    150-249 206-249 44      edge-L  1x      206000  249000  32.34   1.71    2.52
	Chr1    150-249 155-204 50      fakeblock       2x      155000  204000  17.05   2.53    1.80
	Chr1    200-299 200-205 6       edge-R  2x      200000  205000  18.52   2.55    1.69
	Chr1    400-499 459-499 41      edge-L  3x      459000  499000  45.35   3.38    1.71
	Chr1    400-499 408-457 50      fakeblock       1x      408000  457000  14.68   1.69    3.09
	Chr1    450-549 450-458 9       edge-R  1x      450000  458000  18.74   1.76    3.34
	Chr1    450-549 450-463 14      edge-R  2x      450000  463000  8.53    2.37    3.33
	Chr1    500-599 572-599 28      edge-L  2x      572000  599000  15.14   2.53    3.37
	Chr1    500-599 522-571 50      fakeblock       3x      522000  571000  7.74    3.41    2.86
	Chr1    550-649 550-571 22      edge-R  3x      550000  571000  19.38   3.48    2.55
	
'Fakeblocks' are where a potential block has a size which is the maximum window size 
(half the chunk size by default). So this block may have a significantly higher copy 
number, and may lie in the middle of a chunk, but it could possibly be extended even 
further in 5' and/or 3' directions. Hits one edge of the current chunk, so we don't know 
if this block should really be extended

Tried creating new command-line options to:

1. control the maximum level of variance within a window (default = 0.04). This is to 
stop calling a block which mostly spans a 2x region (for instance) but overlaps into a 
1x or 3x region by a small amount. 
2. control the minimum difference in means (within window vs outside window). Default = 0.6

Will now include this info in file names:

	./peak_finder2.pl --shuffle_file CGRs_1kb_plots_forkeith.tsv.shuffled.100.100 --file CGRs_1kb_plots_forkeith.tsv --chunk_size 100 --mean_mode 1  > run_c100s100m1v0.04d0.6.tsv 2>errors	
	
---

## Approach 2: quicker and bigger ##

The above approach takes too long to run exhaustively across many different chunk/window
sizes. Tried new approach to simply look for the breaks at a block boundary rather than
finding the entire block.


	breakpoint_finder.pl --shuffle_file CGRs_1kb_plots_forkeith.tsv.shuffled.1000.100 --file CGRs_1kb_plots_forkeith.tsv --win_size 20 > test20.tsv 2>errors20 &



## Overlap analysis ##

Want to compare the locations of our junctions with all known genomic features from TAIRs
GFF files. 

>*Null hypothesis* - Regions surrounding breakpoints are not enriched for any particular
genomic feature compared to non-breakpoint regions

### Plan of action ###

1. Take a list of current junction coordinates (2 for each breakpoint) in GFF format
2. Assemble a master GFF feature file for A. thaliana
3. For each feature in file 2, compare junctions in file 1 to determine:
	1. Are there more features (per Kbp) inside junctions than outside junctions?
	2. May need to shuffle chromosomes to assess significance

Will need to choose threshold for what constitutes a 'junction region'. Will start by 
using 10 Kbp around junction coordinate (5 Kbp each side).


### Making Master GFF file ###

TAIR has separate GFF files for various features on their FTP (not all in one master file).
So:

	cat TAIR10_GFF3_genes_transposons.gff TAIR_GFF3_ssrs.gff Gutierrez_DNA_replication_origin_TAIR10_GBROWSE.gff > all_TAIR10_features.gff

How many unique features do we have?

	cut -f 2,3 all_TAIR10_features.gff | sort -u
	GutierrezLab	DNA_replication_origin
	TAIR10	CDS
	TAIR10	chromosome
	TAIR10	exon
	TAIR10	five_prime_UTR
	TAIR10	gene
	TAIR10	mRNA
	TAIR10	miRNA
	TAIR10	ncRNA
	TAIR10	protein
	TAIR10	pseudogene
	TAIR10	pseudogenic_exon
	TAIR10	pseudogenic_transcript
	TAIR10	rRNA
	TAIR10	snRNA
	TAIR10	snoRNA
	TAIR10	tRNA
	TAIR10	three_prime_UTR
	TAIR10	transposable_element
	TAIR10	transposable_element_gene
	TAIR10	transposon_fragment
	TandemRepeatsFinder_v4.04	satellite

Can reorder master GFF file to sort by chromosome and then coordinate:

	(sort -k1,1 -k4n,4n all_TAIR10_features.gff  > tmp) && mv tmp all_TAIR10_features.gff
	
	
### Running analysis ### 	

Certain GFF features (chromosome, rRNA, snRNA) are ignored from analysis.

If running in default mode it will just calculate ratios of bp-of-GFF-feature that overlaps
junction region vs bp-of-GFF-feature that occur outside junction region.

	overlap_between_two_gff_files.pl --junction_gff junctions_FRAG00062.gff --feature_gff all_TAIR10_features.gff
	
The result will look something like this:

	Ratio	GFF_feature	Inside_junction_overlap	Outside_junction_overlap
	1.45	DNA_replication_origin	32483/648271 bp (%5.01)	1027907/29711200 bp (%3.46)
	1.20	three_prime_UTR	33578/648271 bp (%5.18)	1284166/29711200 bp (%4.32)
	1.20	protein	347914/648271 bp (%53.67)	13251865/29711200 bp (%44.60)
	1.19	gene	405999/648271 bp (%62.63)	15597770/29711200 bp (%52.50)
	1.18	CDS	222609/648271 bp (%34.34)	8653725/29711200 bp (%29.13)
	1.13	mRNA	421893/648271 bp (%65.08)	17156532/29711200 bp (%57.74)
	1.07	exon	287654/648271 bp (%44.37)	12317229/29711200 bp (%41.46)
	1.04	miRNA	207/648271 bp (%0.03)	9129/29711200 bp (%0.03)
	0.99	five_prime_UTR	18016/648271 bp (%2.78)	831111/29711200 bp (%2.80)
	0.68	satellite	27885/648271 bp (%4.30)	1876037/29711200 bp (%6.31)
	0.63	ncRNA	2507/648271 bp (%0.39)	181189/29711200 bp (%0.61)
	0.55	transposon_fragment	59012/648271 bp (%9.10)	4883365/29711200 bp (%16.44)
	0.55	transposable_element	60168/648271 bp (%9.28)	5011077/29711200 bp (%16.87)
	0.47	transposable_element_gene	17310/648271 bp (%2.67)	1681100/29711200 bp (%5.66)
	0.22	tRNA	72/648271 bp (%0.01)	15268/29711200 bp (%0.05)
	0.22	pseudogenic_exon	1009/648271 bp (%0.16)	206478/29711200 bp (%0.69)
	0.21	pseudogene	1009/648271 bp (%0.16)	224494/29711200 bp (%0.76)
	0.21	pseudogenic_transcript	1009/648271 bp (%0.16)	224494/29711200 bp (%0.76)
	0.00	snoRNA	0/648271 bp (%0.00)	2174/29711200 bp (%0.01)

This suggests that DNA replication origins are 1.45 x more likely to occur inside junction
regions that outside junction regions (in FRAG00062). Other GFF features are also over/under
represented in junction regions.

Now to include shuffling to see whether we see similar ratios when we randomize the location
of all of the breakpoints (for the tailswap region we allow the possibility of all junctions 
occurring in Chr1, or Chr4, or any combination of both).

	overlap_between_two_gff_files.pl --junction_gff junctions_FRAG00062.gff --feature_gff all_TAIR10_features.gff --shuffles 100 --verbose > FRAG00062_junction_output.tsv

Shuffling results suggest that enrichment of gene features in junction regions is 
significant, but enrichment of DNA replication origins is not. For each feature, the 
table below shows how many times (out of a 1000 shuffles of the junction coordinates),
the observed ratio in the real data was exceeded, equalled, or not exceeded.

	DNA_replication_origin	1.4483	289	0	711
	protein	1.2033	0	0	1000
	three_prime_UTR	1.1984	26	1	973
	gene	1.1930	0	0	1000
	CDS	1.1790	0	0	1000
	mRNA	1.1270	7	0	993
	exon	1.0703	120	1	879
	miRNA	1.0392	226	1	773
	five_prime_UTR	0.9935	529	1	470
	satellite	0.6812	998	0	2
	ncRNA	0.6341	666	0	334
	transposon_fragment	0.5538	1000	0	0
	transposable_element	0.5503	1000	0	0
	transposable_element_gene	0.4719	996	0	4
	pseudogenic_exon	0.2240	924	0	76
	tRNA	0.2161	819	0	181
	pseudogenic_transcript	0.2060	931	0	69
	pseudogene	0.2060	931	0	69
	snoRNA	0.0000	358	642	0

### Checking gene orientation ###

At this point we realized that we would like to know whether the enrichment of genes 
inside junction regions followed any pattern. I.e. are there more likely to be convergently
transcribed genes (hence more 3' UTRs) than divergently or tandemly transcribed genes.

I put only gene features from all_TAIR10_feature.gff into a new file (genes.gff) and made
a new script to test this:

	./check_gene_orientation.pl --junction_gff junctions_FRAG00062.gff --feature_gff genes.gff
	
This revealed that across the entire genome, (protein-coding) gene orientation is 
effectively random (as expected):

	>>      6987    %25.82
	<>      6508    %24.05
	<<      7058    %26.08
	><      6508    %24.05
	
I then looked at genes inside junction regions. The breakpoint itself can either be in a 
gene or between genes:

	>>>|>>> 18      %24.32
	<<<|<<< 26      %35.14	
	>>>---|--->>>   6       %8.11
	<<<---|---<<<   9       %12.16
	>>>---|---<<<   4       %5.41
	<<<---|--->>>   11      %14.86
	
So 44 out of 74 breakpoints (59%) are inside genes. The remaining 30 show a slight increase
towards being between divergently transcribed genes (37% (11/30) of all intergenic 
breakpoints), but can't do much with small data size.

Re-ran this script with newer breakpoint data for FRAG00062 (now with 90 breakpoints):

	>>>|>>>	25	%27.78
	<<<|<<<	31	%34.44
	>>>---|---<<<	4	%4.44
	<<<---|---<<<	8	%8.89
	<<<---|--->>>	14	%15.56
	>>>---|--->>>	8	%8.89

Breakpoints inside genes now up to 62% of all breakpoints (56/90) and those that are inside
divergently transcribed genes now up to 41% (14/34).


### New nomenclature and data formats needed ###

Decided to work on reorganizing our data to make it easier going forward. To clarify, we
have contigs/blocks representing duplicated (or triplicated regions). We can represent
them using the Sequence Ontology term 'copy_number_gain' (SO:0001742):

http://www.sequenceontology.org/browser/current_svn/term/SO:0001742

Each block has two ends which define breakpoints. These can be represented with the
Sequence Ontology term 'chromosome_breakpoint' (SO:0001021):

http://www.sequenceontology.org/browser/current_svn/term/SO:0001021 

Each breakpoint has a connected breakpoint (the end of another duplicated block). In most
cases these pairs of breakpoints will contain intervening sequence that is not present
in the reference. Collectively, a pair of breakpoints (and inserted sequence) defines a 
junction. These could potentially also be represented in Sequence Ontology terms with
'insertion_site' (SO:0000366):

http://www.sequenceontology.org/browser/current_svn/term/SO:0000366

We now need anonymous (and unique) identifiers for blocks, breakpoints and junctions. 
E.g. block0001, breakpoint0001 etc. Each block (copy_number_gain) should connect to two
breakpoint objects, and each breakpoint should have a parent block ID, and a paired 
breakpoint ID. Both blocks and breakpoints will be stored in a single GFF file. 

Junctions will point to the flanking breakpoint IDs but will probably not be stored in a 
GFF file (we could list the pair of insertion sites...which will be the same 
coordinates as the breakpoints, but it gets confusing when we are dealing with sequences
that differ from the reference genome). We'll use a spreadsheet instead so that we can 
store the sequence that occurs in junction regions.


	##gff-version 3
	##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=3702
	##genome-build TAIR TAIR10
	#Data from A. thaliana line FRAG00062
	#Generated 20140402
	Chr1	t_test.pl	copy_number_gain	1	205575	.	+	.	ID=block0001;Name=01a1_01a2;Note="duplicated block"
	Chr1	PRICE	chromosome_breakpoint	1	1	.	.	.	ID=breakpoint0001;Parent=block0001;Name=01a1;Note="telomeric end"
	Chr1	PRICE	chromosome_breakpoint	205575	205575	.	-	.	ID=breakpoint0002;Parent=block0001;Name=01a2;Note="paired with breakpoint0033"
	Chr1	t_test.pl	copy_number_gain	31632	87466	.	+	.	ID=block0002;Name=01b1_01b2;Note="triplicated block"
	Chr1	PRICE	chromosome_breakpoint	31632	31632	.	-	.	ID=breakpoint0003;Parent=block0002;Name=01b1;Note="paired with breakpoint0087"


## Data on nearest feature to each breakpoint ##

Using new GFF file, I wrote a script to calculate the average distance of any genomic feature
to each breakpoint. I.e. for every breakpoint find nearest gene/UTR/satellite etc. Then 
average nearest distances across all breakpoints. 

	./nearest_feature.pl --junction_gff FRAG00062.gff --feature_gff all_TAIR10_features.gff

	Feature	Average_distance_to_nearest_breakpoint	Standard_deviation	Number_of_features
	satellite	380	478	73112
	exon	604	1193	57589
	CDS	911	2246	53113
	mRNA	952	1214	9953
	protein	1195	2215	9271
	protein_coding_gene	1218	2195	7092
	three_prime_UTR	2104	2796	8185
	five_prime_UTR	2356	3058	9224
	transposable_element	6257	8084	7107
	transposon_fragment	6257	8084	7858
	non_protein_coding_gene	44014	33069	432
	DNA_replication_origin	44476	58616	376
	ncRNA	108457	98168	149
	pseudogene	136192	279742	233
	pseudogenic_exon	136192	279742	327
	pseudogenic_transcript	136192	279742	233
	transposable_element_gene	142647	166882	681
	tRNA	157066	153789	236
	miRNA	254926	247939	49
	snoRNA	2082570	2450550	22

Results are probably biased towards higher density of certain features. I.e. breakpoints
are most likely to be nearest a satellite feature, but there are more satellite features
than anything else.


## Updated overlap analysis with new data ##

Now updating earlier analysis with more data (90 breakpoints) and using new GFF file. 
Note new names of command-line options. Two runs with different sizes of breakpoint regions 
(1 & 10 Kbp):

	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062.gff --feature_gff all_TAIR10_features.gff --shuffles 1000 --verbose --bp 1000 > FRAG00062_breakpoints_s1000_L1000.tsv
	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062.gff --feature_gff all_TAIR10_features.gff --shuffles 1000 --verbose > FRAG00062_breakpoints_s1000_L10000.tsv

Results from region of 1,000 bp:

	Ratio   GFF_feature     Inside_breakpoint_overlap       Outside_breakpoint_overlap
	1.8220  three_prime_UTR 6888/87307 bp (%7.89)   1310856/30272664 bp (%4.33)
	1.6351  tRNA    72/87307 bp (%0.08)     15268/30272664 bp (%0.05)
	1.1862  five_prime_UTR  2895/87307 bp (%3.32)   846232/30272664 bp (%2.80)
	1.1476  gene    52793/87307 bp (%60.47) 15950976/30272664 bp (%52.69)
	1.1142  protein 43561/87307 bp (%49.89) 13556218/30272664 bp (%44.78)
	1.1004  mRNA    55610/87307 bp (%63.69) 17522815/30272664 bp (%57.88)
	1.0909  exon    39534/87307 bp (%45.28) 12565349/30272664 bp (%41.51)
	1.0765  CDS     27473/87307 bp (%31.47) 8848861/30272664 bp (%29.23)
	0.9328  DNA_replication_origin  2845/87307 bp (%3.26)   1057545/30272664 bp (%3.49)
	0.7439  satellite       4076/87307 bp (%4.67)   1899846/30272664 bp (%6.28)
	0.6135  transposable_element_gene       3000/87307 bp (%3.44)   1695410/30272664 bp (%5.60)
	0.5809  transposable_element    8482/87307 bp (%9.72)   5062763/30272664 bp (%16.72)
	0.5551  transposon_fragment     7899/87307 bp (%9.05)   4934478/30272664 bp (%16.30)
	0.2096  ncRNA   111/87307 bp (%0.13)    183585/30272664 bp (%0.61)
	0.0000  pseudogenic_transcript  0/87307 bp (%0.00)      225503/30272664 bp (%0.74)
	0.0000  pseudogene      0/87307 bp (%0.00)      225503/30272664 bp (%0.74)
	0.0000  snoRNA  0/87307 bp (%0.00)      2174/30272664 bp (%0.01)
	0.0000  pseudogenic_exon        0/87307 bp (%0.00)      207487/30272664 bp (%0.69)
	0.0000  miRNA   0/87307 bp (%0.00)      9336/30272664 bp (%0.03)

And now results after shuffling files 250 times (script still running):

        three_prime_UTR 1.8220  1       0       249
		gene    1.1476  6       0       244
        snoRNA  0.0000  17      233     0
		miRNA   0.0000  23      227     0
        tRNA    1.6351  35      0       215
        five_prime_UTR  1.1862  62      0       188
        ncRNA   0.2096  136     0       114
		pseudogene      0.0000  183     67      0
        DNA_replication_origin  0.9328  171     0       79
    	pseudogenic_exon        0.0000  183     67      0
        pseudogenic_transcript  0.0000  183     67      0
        exon    1.0909  39      0       211
        CDS     1.0765  38      0       212
        transposable_element_gene       0.6135  219     0       31
		mRNA    1.1004  18      0       232
        protein 1.1142  16      0       234
        satellite       0.7439  242     0       8
        transposable_element    0.5809  249     0       1
        transposon_fragment     0.5551  250     0       0

Only 3' UTR and gene features seem to be significant at one end, and transposons at the 
other end. Results seem more significant when you switch to using a larger region around
the breakpoint (10 Kbp):

	less FRAG00062_breakpoints_s1000_L1000.tsv | grep '250/100'  | cut -f 2-6 | sort -nk 3

	Ratio   GFF_feature     Inside_breakpoint_overlap       Outside_breakpoint_overlap
	1.2774  DNA_replication_origin  33668/759948 bp (%4.43) 1026722/29604523 bp (%3.47)
	1.2578  miRNA   292/759948 bp (%0.04)   9044/29604523 bp (%0.03)
	1.1566  protein 392120/759948 bp (%51.60)       13207659/29604523 bp (%44.61)
	1.1495  three_prime_UTR 37770/759948 bp (%4.97) 1279974/29604523 bp (%4.32)
	1.1447  gene    456841/759948 bp (%60.11)       15546928/29604523 bp (%52.52)
	1.1276  CDS     249705/759948 bp (%32.86)       8626629/29604523 bp (%29.14)
	1.0783  mRNA    473476/759948 bp (%62.30)       17104949/29604523 bp (%57.78)
	1.0341  five_prime_UTR  21958/759948 bp (%2.89) 827169/29604523 bp (%2.79)
	1.0223  exon    322310/759948 bp (%42.41)       12282573/29604523 bp (%41.49)
	0.9621  pseudogenic_exon        5001/759948 bp (%0.66)  202486/29604523 bp (%0.68)
	0.8835  pseudogene      5001/759948 bp (%0.66)  220502/29604523 bp (%0.74)
	0.8835  pseudogenic_transcript  5001/759948 bp (%0.66)  220502/29604523 bp (%0.74)
	0.7181  satellite       34461/759948 bp (%4.53) 1869461/29604523 bp (%6.31)
	0.6098  transposon_fragment     76172/759948 bp (%10.02)        4866205/29604523 bp (%16.44)
	0.6048  transposable_element    77527/759948 bp (%10.20)        4993718/29604523 bp (%16.87)
	0.5390  ncRNA   2507/759948 bp (%0.33)  181189/29604523 bp (%0.61)
	0.4205  transposable_element_gene       18136/759948 bp (%2.39) 1680274/29604523 bp (%5.68)
	0.1837  tRNA    72/759948 bp (%0.01)    15268/29604523 bp (%0.05)
	0.0000  snoRNA  0/759948 bp (%0.00)     2174/29604523 bp (%0.01)

The ratios are lower, and many features have different ratios (3' UTR feature drops down 
the list and DNA replication origin moves up). And after 250 shuffles:

	less FRAG00062_breakpoints_s1000_L10000.tsv | grep '250/100'  | cut -f 2-6 | sort -nk 3
	gene	1.1447	0	0	250
	protein	1.1566	0	0	250
	CDS	1.1276	1	0	249
	three_prime_UTR	1.1495	3	0	247
	mRNA	1.0783	9	0	241
	miRNA	1.2578	35	0	215
	five_prime_UTR	1.0341	92	0	158
	pseudogenic_exon	0.9621	93	0	157
	DNA_replication_origin	1.2774	94	0	156
	exon	1.0223	95	0	155
	pseudogene	0.8835	105	0	145
	pseudogenic_transcript	0.8835	105	0	145
	snoRNA	0.0000	125	125	0
	ncRNA	0.5390	180	0	70
	tRNA	0.1837	225	0	25
	satellite	0.7181	250	0	0
	transposable_element	0.6048	250	0	0
	transposable_element_gene	0.4205	250	0	0
	transposon_fragment	0.6098	250	0	0

Now we see that gene, protein, CDS, 3' UTR, and mRNA all see significantly enriched in the
(larger) breakpoint regions. Unclear as to why there is such a difference between 5' and
3' UTR features, unless breakpoints are more likely to occur nearer to 3' end of genes?

Also want to repeat this analysis but separately on breakpoints that flank either 
duplicated or triplicated blocks. Need to extract these GFF breakpoint features into 
separate files (I guess I could make the script do the work):

	grep duplicated FRAG00062.gff   | sed 's/.*ID=\(block[0-9]*\);N.*/\1/' > duplicated_blocks.txt
	grep triplicated FRAG00062.gff  | sed 's/.*ID=\(block[0-9]*\);N.*/\1/' > triplicated_blocks.txt
	grep -f duplicated_blocks.txt FRAG00062.gff > FRAG00062_2x.gff
	grep -f triplicated_blocks.txt FRAG00062.gff > FRAG00062_3x.gff
	
Now re-run overlap script (but maybe with just 1,000 shuffles). Han expects there to be a
difference with more notable enrichment of genomic features around 3x blocks rather than 2x.
Will use default breakpoint region size of 10,000 bp:

	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff all_TAIR10_features.gff --shuffles 100 --verbose > FRAG00062_2x_breakpoints_s100_L10000.tsv
	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff all_TAIR10_features.gff --shuffles 100 --verbose > FRAG00062_3x_breakpoints_s100_L10000.tsv


### Try to look for overlap between breakpoint regions and DNase I hypersensitive sites (DHS) ###

This paper by Zhang et al. (2012) describes a set of hypersensitive sites in A. thaliana:

[Genome-Wide Identification of Regulatory DNA Elements and Protein-Binding Footprints 
Using Signatures of Open Chromatin in Arabidopsis](http://www.plantcell.org/content/24/7/2719.full)

In this paper, they generate DHS maps for seed and flower tissues. These data were submitted
to the GEO database and are available under accession [GSE34318](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34318).
One of the data files in in BED format (GSE34318_dhsites_region.bed) and this includes the 
location of DHS regions (presumably after some threshold value has been exceeded) for both 
leaf and seed tissues.

I extracted just the DHS regions that were identified in wildtype leaf libraries, and 
converted to GFF.

First quick replacement of spaces with tabs in downloaded bed file and then convert to 
desired GFF format with simple Perl script:

	tr ' ' '\t' < GSE34318_dhsites_region.bed  | grep wtleaf| sort -k 1,1 -k 2n,2n   > GSE34318_dhsites_region_sorted.bed
	./bed2gff.pl GSE34318_dhsites_region_sorted.bed  > DHS.gff






