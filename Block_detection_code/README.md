FRAG_project
============

To contain code and data files for the peak calling part of the Fragmentation project, 
investigating shattered chromosome phenotypes in _Arabidopsis thaliana_.
	
	
## Background ## 

Raw data from Han's mapping of reads back to reference genome is in the file:

	CGRs_1kb_plots_forkeith.tsv

This file contains raw and normalized counts of reads mapped to 1 Kbp bins in three lines 
(FRAG00062, FRAG000133, and FRAG00080 - a control line). The normalized counts don't just use 
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

This project was paused as this point. The code had some issues that still needed to be
worked on.
