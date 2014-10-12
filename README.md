FRAG_project
============

# Overview

To contain code and data files for the Fragmentation project, investigating shattered 
chromosome phenotypes in _Arabidopsis thaliana_. The main focus of this repository
is to develop code to help determine breakpoints between duplicated and triplicated
regions of chromosomes (see 'Block_detection_code/' directory) and also to see whether
certain genomic features are over- or under-represented at either the breakpoint regions
or in the blocks themselves.
	

# Generating Data: Making Master GFF file of all genomic features #

## 1) Get GFF features from TAIR ##

TAIR has separate GFF files for various features on their FTP site (not all in one 
master file). So:

	cat TAIR10_GFF3_genes_transposons.gff TAIR_GFF3_ssrs.gff > all_TAIR10_features.gff

	


## 2) Making replication origin GFF file ##

I wrote a simple script `ori2gff.pl` that takes the raw data from the Gutierrez et al.
paper and converts it to a GFF format. E.g. 	

	./ori2gff.pl Gutierrez.txt > Gutierrez_DNA_replication_origin_TAIR10_GBROWSE.gff
	
This file will be appended to the main TAIR10 GFF file of all genomic features.

	cat Gutierrez_DNA_replication_origin_TAIR10_GBROWSE.gff  >> all_TAIR10_features.gff




## 3) Making DNAse I hypersensitive site GFF file ##

This paper by Zhang et al. (2012) describes a set of hypersensitive sites in A. thaliana:

[Genome-Wide Identification of Regulatory DNA Elements and Protein-Binding Footprints 
Using Signatures of Open Chromatin in Arabidopsis](http://www.plantcell.org/content/24/7/2719.full)

In this paper, they generate DHS maps for seed and flower tissues. These data were 
submitted to the GEO database and are available under accession 
[GSE34318](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34318). One of the data 
files in in BED format (GSE34318_dhsites_region.bed) and this includes the location of 
DHS regions (presumably after some threshold value has been exceeded) for both leaf and 
seed tissues.

I extracted just the DHS regions that were identified in wildtype leaf libraries, and 
converted to GFF.

First quick replacement of spaces with tabs in downloaded bed file and then convert to 
desired GFF format with simple Perl script:

	tr ' ' '\t' < GSE34318_dhsites_region.bed  | grep wtleaf | sort -k 1,1 -k 2n,2n  > GSE34318_dhsites_region_sorted.bed
	./bed2gff.pl GSE34318_dhsites_region_sorted.bed  > DHS.gff
	
Now combine with main file:

	cat DHS.gff >> all_TAIR10_features.gff
	grep -vE "^#" all_TAIR10_features.gff | sort -k 1,1 -k 4n,4n  all_TAIR10_features.gff > tmp.gff; mv tmp.gff all_TAIR10_features.gff	




### 4) Making Chromatin state information GFF file ###

This new paper by Sequeira-Mendes et al. (2014) describes a set of distinct chromatin 
states in *A. thaliana*: 

[The Functional Topography of the Arabidopsis Genome Is Organized in a Reduced Number of 
Linear Motifs of Chromatin States](http://www.plantcell.org/content/early/2014/06/11/tpc.114.124578.long)

They mostly use existing histone and other epigenetic modification data, but have some of
their own. Using principal components analysis, they end up defining 9 different states
of chromatin. These are available in supplemental data file 2 (Excel spreadsheet).

I converted this to a GFF file. Sequence Ontology only has one term for chromatin
(SO:0001747 open_chromatin_state) so I will use the 9th column of GFF to distinguish the
9 different states.

First export each tab in Excel spreadsheet to a text file (Windows formatted text), e.g
state1.txt, state2.txt etc. Then simple Perl script to convert them to GFF.

	./chromatin_state2gff.pl > tmp.gff
	cat tmp.gff >> all_TAIR10_features.gff
	sort -k 1,1 -k 4n,4n  all_TAIR10_features.gff > tmp.gff; mv tmp.gff all_TAIR10_features.gff	


## Summary of TAIR10 GFF data ##

How many unique features do we have?

	cut -f 2,3 all_TAIR10_features.gff | sort -u
	GEO     DNAseI_hypersensitive_site
	GutierrezLab    DNA_replication_origin
	Sequeira-Mendes_2014_paper      open_chromatin_state
	TAIR10  CDS
	TAIR10  chromosome
	TAIR10  exon
	TAIR10  five_prime_UTR
	TAIR10  gene
	TAIR10  miRNA
	TAIR10  mRNA
	TAIR10  ncRNA
	TAIR10  protein
	TAIR10  pseudogene
	TAIR10  pseudogenic_exon
	TAIR10  pseudogenic_transcript
	TAIR10  rRNA
	TAIR10  snoRNA
	TAIR10  snRNA
	TAIR10  three_prime_UTR
	TAIR10  transposable_element
	TAIR10  transposable_element_gene
	TAIR10  transposon_fragment
	TAIR10  tRNA
	TandemRepeatsFinder_v4.04       satellite

Can reorder master GFF file to sort by chromosome and then coordinate:

	(sort -k1,1 -k4n,4n all_TAIR10_features.gff  > tmp) && mv tmp all_TAIR10_features.gff




# Generating a GFf file to represent experimental data



## Nomenclature and GFF data formats used to represent breakpoint/block data ##

To summarize, we have regions of the genome that have been identified as representing 
duplicated (or triplicated) 'blocks'. We can represent these regions by using the Sequence 
Ontology term 'copy_number_gain' (SO:0001742):

http://www.sequenceontology.org/browser/current_svn/term/SO:0001742

Each duplicated or triplicated block Ñ henceforth referred to as either 2x or 3x for
simplicity ÑÊhas two ends which define breakpoints. These can be represented with the
Sequence Ontology term 'chromosome_breakpoint' (SO:0001021):

http://www.sequenceontology.org/browser/current_svn/term/SO:0001021 

Each breakpoint has a connected breakpoint (the end of another 2x/3x block). In most
cases these pairs of breakpoints will contain intervening sequence that is not present
in the reference. Collectively, a pair of breakpoints (and inserted sequence) defines a 
junction. These could potentially also be represented in Sequence Ontology terms with
'insertion_site' (SO:0000366):

http://www.sequenceontology.org/browser/current_svn/term/SO:0000366

Finally, there are also the regions of the genome that punctuate 2x and 3x blocks and which
can be thought of as unduplicated regions, or '1x' blocks. These can be represented with
the Sequence Ontology term 'region' (SO:0000001):

http://www.sequenceontology.org/browser/current_svn/term/SO:0000001


## Making a new GFF file ##

We now need anonymous (and unique) identifiers for blocks and breakpoints. E.g. block0001, 
breakpoint0001 etc. Each block (copy_number_gain) should connect to two breakpoint objects, 
and each breakpoint should have a parent block ID, and a paired breakpoint ID. Both blocks 
and breakpoints will be stored in a single GFF file. 

Example GFF file:

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
	Chr1    .       region  205575  458539  .       +       .       ID=block0091;Note="no copy number gain detected"



## Treating FRAG00062 data by copy-number level ##

For FRAG00062 we can also separate out the breakpoint data into two subsets, those 
breakpoints that flank either 2x or 3x blocks:

	grep duplicated FRAG00062.gff   | sed 's/.*ID=\(block[0-9]*\);N.*/\1/' > duplicated_blocks.txt
	grep triplicated FRAG00062.gff  | sed 's/.*ID=\(block[0-9]*\);N.*/\1/' > triplicated_blocks.txt
	grep -f duplicated_blocks.txt FRAG00062.gff > FRAG00062_2x.gff
	grep -f triplicated_blocks.txt FRAG00062.gff > FRAG00062_3x.gff








# Analysis: part 1



### Checking feature orientation ###

At this point we realized that we would like to know whether the enrichment of genes 
inside junction regions followed any pattern. I.e. are there more likely to be convergently
transcribed genes (hence more 3' UTRs) than divergently or tandemly transcribed genes.

I put only gene features from all_TAIR10_feature.gff into a new file (genes.gff) and made
a new script to test this. Running this script with newer breakpoint data for FRAG00062 
(90 breakpoints, excluding the 2 that meet telomeres):

	./check_feature_orientation.pl --breakpoint_gff GFF_files/FRAG00062.gff --feature_gff GFF_files/genes.gff --target gene

	>>>|>>>	25	%27.78
	<<<|<<<	31	%34.44
	>>>---|---<<<	4	%4.44
	<<<---|---<<<	8	%8.89
	<<<---|--->>>	14	%15.56
	>>>---|--->>>	8	%8.89

Breakpoints inside genes account for 62% of all breakpoints (56/90) and those that are inside
divergently transcribed genes account for 41% of all intergenic breakpoints (14/34).

And for FRAG00045 (30 breakpoints):

	./check_feature_orientation.pl --breakpoint_gff FRAG00045.gff --feature_gff genes.gff --target gene

	>>>|>>>	12	%40.00
	<<<|<<<	5	%16.67
	>>>---|---<<<	3	%10.00
	<<<---|---<<<	3	%10.00
	>>>---|--->>>	5	%16.67
	<<<---|--->>>	2	%6.67

This time only 57% of breakpoints were inside genes (17/30) and only 15% of intergenic 
breakpoints are between divergently transcribed genes. So no strong pattern.





### Data on nearest feature to each breakpoint ###

I wrote a script to calculate the average distance of any genomic feature to each 
breakpoint. I.e. for every breakpoint find nearest gene/UTR/satellite etc. Then averaged
the nearest distances across all breakpoints. 

	./nearest_feature.pl --breakpoint_gff GFF_files/FRAG00062.gff --feature_gff GFF_files/all_TAIR10_features.gff

	Feature Average_distance_to_nearest_breakpoint  Standard_deviation      Number_of_features
	CDS     879     2251    53113
	DNA_replication_origin  44837   59143   376
	DNAseI_hypersensitive_site      1507    2403    10187
	chromosome      7407070 3791549 1
	exon    567     1162    57589
	five_prime_UTR  2325    3085    9224
	mRNA    922     1194    9953
	miRNA   257654  248174  49
	ncRNA   108994  98485   149
	non_protein_coding_gene 44618   33137   432
	open_chromatin_state_1  3126    3856    3155
	open_chromatin_state_2  2581    3245    3833
	open_chromatin_state_3  3749    4999    3056
	open_chromatin_state_4  2834    3494    4050
	open_chromatin_state_5  9277    11084   1897
	open_chromatin_state_6  5153    6944    2303
	open_chromatin_state_7  12548   15714   1248
	open_chromatin_state_8  27134   31662   1275
	open_chromatin_state_9  343955  388532  468
	protein 1169    2223    9271
	protein_coding_gene     1194    2205    7092
	pseudogene      116394  230257  233
	pseudogenic_exon        116394  230257  327
	pseudogenic_transcript  116394  230257  233
	satellite       384     482     73112
	snRNA   4342430 2745429 2
	snoRNA  2116014 2467006 22
	tRNA    158832  154954  236
	three_prime_UTR 2039    2791    8185
	transposable_element    6236    8143    7107
	transposable_element_gene       131352  146540  681
	transposon_fragment     6236    8143    7858

Results are probably biased towards higher density of certain features. I.e. breakpoints
are most likely to be nearest a satellite feature, but there are more satellite features
than anything else.






# Analysis: part 2




## Analysis of enriched features in breakpoint regions ##

If we define a 'breakpoint region' as a window of sequence around each breakpoint location
(mapped to the reference genome), we can ask whether any genomic features are enriched
in these breakpoint regions. E.g. take 1,000 bp around all breakpoints (potentially
overlapping other breakpoints) and ask whether the total bp of a feature such as 'coding
exons' is higher (as a percentage) *inside* those regions vs all DNA *outside* those
regions.

Can try this for many different sizes of breakpoint region (100 bp up to 50 Kbp). In 
this analysis certain GFF features (e.g. chromosome) are ignored. The final result is 
calculated as a ratio of %breakpoint-region-occupied-by-feature compared to
%non-breakpoint-region-occupied-by-feature. E.g.

 ./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062.gff  --feature_gff GFF_files/all_TAIR10_features.gff --bp 100

	Run     Real_ratio      Bp      Feature Breakpoint_region_bp    Non_breakpoint_region_bp        Feature_bp_inside       %Inside Feature_bp_outside      %Outside        Shuffled_ratio  Above   Same    Below
	0       1.1372  100     CDS     8957    30350564        2978    33.25   8873356 29.24   1.1372  0       0       0
	0       0.9589  100     DNA_replication_origin  8957    30350564        300     3.35    1060090 3.49    0.9589  0       0       0
	0       1.4907  100     DNAseI_hypersensitive_site      8957    30350564        1749    19.53   3975577 13.10   1.4907  0       0       0
	0       1.1628  100     exon    8957    30350564        4324    48.28   12600559        41.52   1.1628  0       0       0
	0       1.4412  100     five_prime_UTR  8957    30350564        361     4.03    848766  2.80    1.4412  0       0       0
	0       1.1497  100     gene    8957    30350564        5428    60.60   15998341        52.71   1.1497  0       0       0
	0       1.0995  100     mRNA    8957    30350564        5702    63.66   17572723        57.90   1.0995  0       0       0
	0       0.0000  100     miRNA   8957    30350564        0       0.00    9336    0.03    0.0000  0       0       0
	0       0.0000  100     ncRNA   8957    30350564        0       0.00    183696  0.61    0.0000  0       0       0
	0       1.2536  100     open_chromatin_state_1  8957    30350564        1275    14.23   3446325 11.36   1.2536  0       0       0
	0       2.0929  100     open_chromatin_state_2  8957    30350564        1910    21.32   3092404 10.19   2.0929  0       0       0
	0       1.1439  100     open_chromatin_state_3  8957    30350564        1053    11.76   3119097 10.28   1.1439  0       0       0
	0       0.9309  100     open_chromatin_state_4  8957    30350564        1309    14.61   4764641 15.70   0.9309  0       0       0
	0       0.6904  100     open_chromatin_state_5  8957    30350564        900     10.05   4417500 14.55   0.6904  0       0       0
	0       0.9582  100     open_chromatin_state_6  8957    30350564        935     10.44   3306415 10.89   0.9582  0       0       0
	0       0.7249  100     open_chromatin_state_7  8957    30350564        669     7.47    3127281 10.30   0.7249  0       0       0
	0       0.9969  100     open_chromatin_state_8  8957    30350564        706     7.88    2399744 7.91    0.9969  0       0       0
	0       0.1266  100     open_chromatin_state_9  8957    30350564        100     1.12    2677100 8.82    0.1266  0       0       0
	0       1.0635  100     protein 8957    30350564        4267    47.64   13595512        44.79   1.0635  0       0       0
	0       0.0000  100     pseudogene      8957    30350564        0       0.00    225503  0.74    0.0000  0       0       0
	0       0.0000  100     pseudogenic_exon        8957    30350564        0       0.00    207487  0.68    0.0000  0       0       0
	0       0.0000  100     pseudogenic_transcript  8957    30350564        0       0.00    225503  0.74    0.0000  0       0       0
	0       0.7886  100     satellite       8957    30350564        443     4.95    1903479 6.27    0.7886  0       0       0
	0       0.0000  100     snRNA   8957    30350564        0       0.00    337     0.00    0.0000  0       0       0
	0       0.0000  100     snoRNA  8957    30350564        0       0.00    2174    0.01    0.0000  0       0       0
	0       5.7529  100     tRNA    8957    30350564        26      0.29    15314   0.05    5.7529  0       0       0
	0       1.8473  100     three_prime_UTR 8957    30350564        718     8.02    1317026 4.34    1.8473  0       0       0
	0       0.5152  100     transposable_element    8957    30350564        771     8.61    5070474 16.71   0.5152  0       0       0
	0       0.5986  100     transposable_element_gene       8957    30350564        300     3.35    1698110 5.59    0.5986  0       0       0
	0       0.4601  100     transposon_fragment     8957    30350564        671     7.49    4941706 16.28   0.4601  0       0       0

The first set of columns of output are as follows:

1. Run number (starts at 0 for unshuffled results)
2. Real ratio from unshuffled data (the ratio of columns 6 & 7)
3. bp of individual breakpoint regions
4. Feature 
5. Breakpoint_region_bp (for small values of bp, this equals number of breakpoints * bp)
6. Non_breakpoint_region_bp      
7. Feature_bp_inside      
8. %Inside 
9. Feature_bp_outside      
10. %Outside        

To assess the significance of these ratios, I perform shuffling experiments to see whether 
we see similar ratios when we randomize the location of all of the breakpoints 
(for the tailswap region we allow the possibility of all junctions occurring in Chr1, 
or Chr4, or any combination of both, but it is proportional to the length of each region
on Chr1 and Chr4).

E.g. if you want to see how significant the above observed enrichment ratio for genes is
(1.1497), you could run:

	./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062.gff  --feature_gff GFF_files/genes.gff --bp 100 --shuffles 1000


You should now note the final four columns of output:

11. Shuffled_ratio  
12. Above (number of times column 9 exceeds column 2)
13. Same (number of times column 9 equals column 2)   
14. Below (number of times column 9 is below column 2)

There will also be rows of output generated for each shuffle run. Just want to really look
at first run (run = 0) and last run in file:

	Run     Real_ratio      Bp      Feature Breakpoint_region_bp    Non_breakpoint_region_bp        Feature_bp_inside       %Inside Feature_bp_outside      %Outside        Shuffled_ratio  Above   Same    Below
	0       1.1497  100     gene    8957    30350564        5428    60.60   15998341        52.71   1.1497  0       0       0
	1000    1.1497  100     gene    9000    30350471        4064    45.16   15999705        52.72   0.8566  62      0       938

This suggests that the observed enrichment ratio (1.1497) was exceeded in 62 out of 1000 
shuffles. Ideally, we want to do this for all genome features, for both 2x and 3x regions,
and for different sizes of breakpoint regions. E.g. set up lots of runs like this:

	./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062.gff    --feature_gff GFF_files/all_TAIR10_features.gff --shuffles 1000 --verbose --bp 10000 > Results/FRAG00062_breakpoints_s1000_L10000.tsv
	./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062_2x.gff --feature_gff GFF_files/all_TAIR10_features.gff --shuffles 1000 --verbose --bp 10000 > Results/FRAG00062_2x_breakpoints_s1000_L10000.tsv
	./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062_3x.gff --feature_gff GFF_files/all_TAIR10_features.gff --shuffles 1000 --verbose --bp 10000 > Results/FRAG00062_3x_breakpoints_s1000_L10000.tsv
    

Here are the principle results for a breakpoint region size of 1 Kbp with 1000 shuffles. The 
last 3 columns of output count how many times the observed ratio from the real data was
exceeded, equalled, or not exceeded/equalled in the 1000 shufflings:

	tail -n 19 FRAG00062_breakpoints_1Kbp.tsv | cut -f 2,3,11-13
	1.0765  CDS     121     0       879
	0.9328  DNA_replication_origin  669     0       331
	1.0909  exon    156     0       844
	1.1862  five_prime_UTR  241     0       759
	1.1476  gene    26      0       974
	1.1004  mRNA    89      0       911
	0.0000  miRNA   129     871     0
	0.2096  ncRNA   574     0       426
	1.1142  protein 61      0       939
	0.0000  pseudogene      672     328     0
	0.0000  pseudogenic_exon        672     328     0
	0.0000  pseudogenic_transcript  672     328     0
	0.7439  satellite       967     0       33
	0.0000  snoRNA  60      940     0
	1.6351  tRNA    149     0       851
	1.8220  three_prime_UTR 3       0       997
	0.5809  transposable_element    997     0       3
	0.6135  transposable_element_gene       935     0       65
	0.5551  transposon_fragment     996     0       4

So in this result file, genes are the only feature that are enriched (P < 0.05) and
various repeat elements are significantly underrepresented (e.g. P < 0.01 for 
transposable_element features). However, this might be a logical consequence of genes
being enriched (two features can't always occupy the same space).

	grep -w gene FRAG00062_breakpoints_1Kbp.tsv  | grep -E "^0" | cut -f 3-9
	gene    87307   30272664        52793   60.47   15950976        52.69

Gene features in this dataset account for ~60% of breakpoint regions compared to ~53%
of non-breakpoint regions. When we extend to a larger size of breakpoint region (10 Kbp),
we see similar results, but more enrichment of replication origins:

	tail -n 19 FRAG00062_breakpoints_10Kbp.tsv | cut -f 2,3,11-13
	1.1276  CDS     5       0       995
	1.2774  DNA_replication_origin  423     0       577
	1.0223  exon    434     1       565
	1.0341  five_prime_UTR  382     0       618
	1.1447  gene    2       0       998
	1.0783  mRNA    38      0       962
	1.2578  miRNA   160     0       840
	0.5390  ncRNA   770     0       230
	1.1566  protein 0       0       1000
	0.8835  pseudogene      379     0       621
	0.9621  pseudogenic_exon        328     0       672
	0.8835  pseudogenic_transcript  379     0       621
	0.7181  satellite       1000    0       0
	0.0000  snoRNA  432     568     0
	0.1837  tRNA    894     0       106
	1.1495  three_prime_UTR 47      0       953
	0.6048  transposable_element    1000    0       0
	0.4205  transposable_element_gene       998     0       2
	0.6098  transposon_fragment     999     0       1

With a larger region size, many gene-related features become significantly enriched
(P < 0.01 for CDS, gene, and protein). Enrichment of replication origins is not significant.
At the largest breakpoint region size that I chose (50 Kbp), gene-related features
remain significantly enriched, even if the actual ratios become much lower:

	tail -n 19 FRAG00062_breakpoints_50Kbp.tsv | cut -f 2,3,11-13
	1.0764  CDS     9       0       991
	0.9004  DNA_replication_origin  879     0       121
	0.9801  exon    878     1       121
	1.1176  five_prime_UTR  54      0       946
	1.0800  gene    8       0       992
	1.0115  mRNA    352     4       644
	0.8165  miRNA   312     0       688
	0.5308  ncRNA   938     0       62
	1.0847  protein 6       0       994
	1.2193  pseudogene      17      0       983
	0.9997  pseudogenic_exon        178     0       822
	1.2193  pseudogenic_transcript  17      0       983
	0.7842  satellite       988     0       12
	0.0000  snoRNA  944     56      0
	1.3975  tRNA    174     0       826
	1.0690  three_prime_UTR 57      3       940
	0.6428  transposable_element    1000    0       0
	0.3878  transposable_element_gene       1000    0       0
	0.6452  transposon_fragment     1000    0       0

The percentage of bp of gene features inside and outside breakpoint regions now narrows
to just 56.4% and 52.2%. At the other extreme, a breakpoint region of just 100 bp still
sees significant enrichment of genes (P < 0.05):

	tail -n 19 FRAG00062_breakpoints_100bp.tsv | cut -f 2,3,11-13
	1.1372  CDS     112     0       888
	0.9589  DNA_replication_origin  573     0       427
	1.1628  exon    67      0       933
	1.4412  five_prime_UTR  185     0       815
	1.1497  gene    25      0       975
	1.0995  mRNA    112     0       888
	0.0000  miRNA   32      968     0
	0.0000  ncRNA   424     576     0
	1.0635  protein 173     0       827
	0.0000  pseudogene      441     559     0
	0.0000  pseudogenic_exon        435     565     0
	0.0000  pseudogenic_transcript  441     559     0
	0.7886  satellite       873     0       127
	0.0000  snoRNA  16      984     0
	5.7529  tRNA    62      0       938
	1.8473  three_prime_UTR 57      0       943
	0.5152  transposable_element    993     0       7
	0.5986  transposable_element_gene       898     0       102
	0.4601  transposon_fragment     998     0       2


#### Conclusion 1: genes are significantly enriched in breakpoint regions of all sizes

Can now look to see whether the same pattern holds true when considering 2x or 3x regions
separately. Again, starting at 1 Kbp breakpoint regions (2x results on top):

	tail -n 19 FRAG00062_2x_breakpoints_1Kbp.tsv | cut -f 2,3,11-13; echo; tail -n 19 FRAG00062_3x_breakpoints_1Kbp.tsv | cut -f 2,3,11-13
	1.1842  CDS     63      0       937
	0.0000  DNA_replication_origin  920     80      0
	1.1536  exon    103     0       897
	1.5011  five_prime_UTR  124     0       876
	1.2893  gene    4       0       996
	1.1736  mRNA    35      0       965
	0.0000  miRNA   52      948     0
	0.0000  ncRNA   396     604     0
	1.2076  protein 37      0       963
	0.0000  pseudogene      434     566     0
	0.0000  pseudogenic_exon        434     566     0
	0.0000  pseudogenic_transcript  434     566     0
	0.5646  satellite       989     0       11
	0.0000  snoRNA  36      964     0
	0.0000  tRNA    221     779     0
	2.1291  three_prime_UTR 5       0       995
	0.3110  transposable_element    997     0       3
	0.0000  transposable_element_gene       991     9       0
	0.3191  transposon_fragment     997     0       3

	0.9882  CDS     446     3       551
	1.9026  DNA_replication_origin  211     0       789
	1.0435  exon    408     0       592
	0.8264  five_prime_UTR  623     0       377
	1.0136  gene    374     0       626
	1.0364  mRNA    385     0       615
	0.0000  miRNA   62      938     0
	0.4276  ncRNA   320     0       680
	1.0320  protein 333     0       667
	0.0000  pseudogene      397     603     0
	0.0000  pseudogenic_exon        397     603     0
	0.0000  pseudogenic_transcript  397     603     0
	0.9180  satellite       678     0       322
	0.0000  snoRNA  23      977     0
	3.3352  tRNA    70      0       930
	1.5318  three_prime_UTR 72      0       928
	0.8548  transposable_element    814     0       186
	1.2515  transposable_element_gene       489     0       511
	0.7935  transposon_fragment     867     0       133

There is a striking difference in the ratios of replication origins, and to a lesser-extent
to gene-related features. Breakpoints that flank duplicated (2x) blocks are more likely to
be enriched for gene-related features, but triplicated blocks (3x) have more replication
origins (though this may not be significant).

At 10 Kbp (though not at other sizes), the increase of replication origins in 3x blocks
appears significant (P < 0.05):

	tail -n 19 FRAG00062_2x_breakpoints_10Kbp.tsv | cut -f 2,3,11-13; echo; tail -n 19 FRAG00062_3x_breakpoints_10Kbp.tsv | cut -f 2,3,11-13
	1.1017  CDS     65      0       935
	0.2077  DNA_replication_origin  986     0       14
	0.9516  exon    836     0       164
	0.9697  five_prime_UTR  572     0       428
	1.1664  gene    6       0       994
	1.0600  mRNA    175     0       825
	0.6745  miRNA   380     0       620
	0.1337  ncRNA   831     0       169
	1.1786  protein 4       0       996
	1.4328  pseudogene      146     0       854
	1.5598  pseudogenic_exon        126     0       874
	1.4328  pseudogenic_transcript  146     0       854
	0.6646  satellite       996     0       4
	0.0000  snoRNA  225     775     0
	0.0000  tRNA    861     139     0
	1.1905  three_prime_UTR 52      0       948
	0.5516  transposable_element    995     0       5
	0.0000  transposable_element_gene       1000    0       0
	0.5631  transposon_fragment     994     0       6

	1.1737  CDS     11      0       989
	2.4161  DNA_replication_origin  41      0       959
	1.1120  exon    55      0       945
	1.0853  five_prime_UTR  291     0       709
	1.1366  gene    30      0       970
	1.1127  mRNA    26      0       974
	1.8390  miRNA   112     0       888
	0.9709  ncRNA   411     0       589
	1.1525  protein 21      0       979
	0.2467  pseudogene      822     0       178
	0.2682  pseudogenic_exon        807     0       193
	0.2467  pseudogenic_transcript  822     0       178
	0.7702  satellite       935     0       65
	0.0000  snoRNA  247     753     0
	0.3825  tRNA    555     0       445
	1.1063  three_prime_UTR 171     0       829
	0.6545  transposable_element    975     0       25
	0.8754  transposable_element_gene       808     0       192
	0.6524  transposon_fragment     973     0       27
	
However, at this region size, 3x regions also appear significantly enriched for gene features
but to a lesser extent in 3x regions. Here are the results for just genes and replication
origin features across all 2x/3x output files:

	tail -n 19 FRAG00062_*x_*.tsv | grep -Ew "==>|gene|DNA_replication_origin" | cut -f 2,3,11-13
	==> FRAG00062_2x_breakpoints_100bp.tsv <==
	0.0000  DNA_replication_origin  865     135     0
	1.3193  gene    4       0       996
	==> FRAG00062_2x_breakpoints_10Kbp.tsv <==
	0.2077  DNA_replication_origin  986     0       14
	1.1664  gene    6       0       994
	==> FRAG00062_2x_breakpoints_1Kbp.tsv <==
	0.0000  DNA_replication_origin  920     80      0
	1.2893  gene    4       0       996
	==> FRAG00062_2x_breakpoints_25Kbp.tsv <==
	0.4650  DNA_replication_origin  968     0       32
	1.0901  gene    46      1       953
	==> FRAG00062_2x_breakpoints_2Kbp.tsv <==
	0.1613  DNA_replication_origin  929     0       71
	1.2647  gene    2       0       998
	==> FRAG00062_2x_breakpoints_500bp.tsv <==
	0.0000  DNA_replication_origin  900     100     0
	1.3100  gene    2       0       998
	==> FRAG00062_2x_breakpoints_50Kbp.tsv <==
	0.4673  DNA_replication_origin  994     0       6
	1.0762  gene    56      0       944
	==> FRAG00062_2x_breakpoints_5Kbp.tsv <==
	0.1925  DNA_replication_origin  963     1       36
	1.1705  gene    15      0       985

	==> FRAG00062_3x_breakpoints_100bp.tsv <==
	1.9523  DNA_replication_origin  125     146     729
	0.9795  gene    483     0       517
	==> FRAG00062_3x_breakpoints_10Kbp.tsv <==
	2.4161  DNA_replication_origin  41      0       959
	1.1366  gene    30      0       970
	==> FRAG00062_3x_breakpoints_1Kbp.tsv <==
	1.9026  DNA_replication_origin  211     0       789
	1.0136  gene    374     0       626
	==> FRAG00062_3x_breakpoints_25Kbp.tsv <==
	1.7445  DNA_replication_origin  165     0       835
	1.1510  gene    3       0       997
	==> FRAG00062_3x_breakpoints_2Kbp.tsv <==
	2.0706  DNA_replication_origin  150     0       850
	1.0800  gene    148     0       852
	==> FRAG00062_3x_breakpoints_500bp.tsv <==
	1.8882  DNA_replication_origin  250     0       750
	0.9880  gene    421     1       578
	==> FRAG00062_3x_breakpoints_50Kbp.tsv <==
	1.3704  DNA_replication_origin  318     0       682
	1.0852  gene    56      0       944
	==> FRAG00062_3x_breakpoints_5Kbp.tsv <==
	2.2055  DNA_replication_origin  104     0       896
	1.1063  gene    82      0       918

Genes are most enriched in 2x breakpoint regions of size 100 bp (1.32 fold increase) and
significantly enriched (P < 0.01) at all sizes apart from 25 Kbp (P < 0.05) and 50 Kbp
(NS). Genes are moderately enriched in all 3x breakpoint regions apart from 100 bp and
only highly significantly enriched (P < 0.01) in 25 Kbp regions (weakly significant, 
P < 0.05 in 10 Kbp regions).

In contrast, replication origins were most enriched in 10 Kbp 3x regions (2.42 fold
increase) and this was the only size where they were significantly enriched (P < 0.05).
In 2x breakpoint regions, replication origins are never enriched.
	
#### Conclusion 2: significant enrichment of genes mostly occurs in 2x regions 

#### Conclusion 3: enrichment of replication origins only occurs in 3x regions


### FRAG00045 results

The FRAG00045 line only consists of duplicated blocks, so no 3x regions are present. In
the 10 Kbp region result file, no significantly enriched features are detected apart 
from 5' UTRs (P < 0.05):

	tail -n 19 FRAG00045_breakpoints_10Kbp.tsv | cut -f 2,3,11-13
	1.0584  CDS     170     0       830
	1.4244  DNA_replication_origin  340     0       660
	1.0134  exon    486     0       514
	1.3360  five_prime_UTR  48      0       952
	1.0818  gene    108     0       892
	1.0445  mRNA    285     0       715
	0.0000  miRNA   282     718     0
	2.1837  ncRNA   108     0       892
	1.0712  protein 137     1       862
	0.0000  pseudogene      846     154     0
	0.0000  pseudogenic_exon        846     154     0
	0.0000  pseudogenic_transcript  846     154     0
	0.9193  satellite       635     0       365
	0.0000  snoRNA  155     845     0
	0.0000  tRNA    729     271     0
	1.1155  three_prime_UTR 199     0       801
	0.6484  transposable_element    952     0       48
	0.6349  transposable_element_gene       900     0       100
	0.6654  transposon_fragment     948     0       52

Here are the results that focus only on genes and replication origins across all breakpoint
region sizes:

	tail -n 19 FRAG00045_*.tsv | grep -Ew "==>|gene|DNA_replication_origin" | cut -f 2,3,11-13
	==> FRAG00045_breakpoints_100bp.tsv <==
	0.9543  DNA_replication_origin  390     328     282
	1.0655  gene    318     1       681
	==> FRAG00045_breakpoints_10Kbp.tsv <==
	1.4244  DNA_replication_origin  340     0       660
	1.0818  gene    108     0       892
	==> FRAG00045_breakpoints_1Kbp.tsv <==
	0.9543  DNA_replication_origin  503     172     325
	1.0242  gene    355     0       645
	==> FRAG00045_breakpoints_25Kbp.tsv <==
	1.2519  DNA_replication_origin  437     0       563
	1.0914  gene    74      0       926
	==> FRAG00045_breakpoints_2Kbp.tsv <==
	0.9543  DNA_replication_origin  492     106     402
	1.0247  gene    365     0       635
	==> FRAG00045_breakpoints_500bp.tsv <==
	0.9543  DNA_replication_origin  451     235     314
	1.0350  gene    330     0       670
	==> FRAG00045_breakpoints_50Kbp.tsv <==
	0.9370  DNA_replication_origin  675     0       325
	1.0927  gene    63      0       937
	==> FRAG00045_breakpoints_5Kbp.tsv <==
	1.4082  DNA_replication_origin  382     0       618
	1.0977  gene    130     0       870

Genes are only enriched (in terms of their ratio being positive) in 5 Kbp, 10 Kbp, 25 Kbp 
region sizes, but none of these are significant. Strangely, 5' UTRs (but no 3' UTRs) are
significantly enriched in some large breakpoint sizes (10 Kbp, 25 Kbp, 50 Kbp) with
the latter two sizes being highly significant (P < 0.01):

	tail -n 19 FRAG00045_*.tsv | grep -E "==>|UTR" | cut -f 2,3,11-13
	==> FRAG00045_breakpoints_100bp.tsv <==
	2.5508  five_prime_UTR  59      2       939
	0.5452  three_prime_UTR 698     1       301
	==> FRAG00045_breakpoints_10Kbp.tsv <==
	1.3360  five_prime_UTR  48      0       952
	1.1155  three_prime_UTR 199     0       801
	==> FRAG00045_breakpoints_1Kbp.tsv <==
	1.7759  five_prime_UTR  75      0       925
	1.0760  three_prime_UTR 373     1       626
	==> FRAG00045_breakpoints_25Kbp.tsv <==
	1.3489  five_prime_UTR  11      0       989
	1.1645  three_prime_UTR 59      0       941
	==> FRAG00045_breakpoints_2Kbp.tsv <==
	1.6390  five_prime_UTR  62      0       938
	1.0842  three_prime_UTR 320     0       680
	==> FRAG00045_breakpoints_500bp.tsv <==
	2.0724  five_prime_UTR  73      1       926
	0.4192  three_prime_UTR 828     1       171
	==> FRAG00045_breakpoints_50Kbp.tsv <==
	1.3177  five_prime_UTR  7       0       993
	1.1445  three_prime_UTR 64      0       936
	==> FRAG00045_breakpoints_5Kbp.tsv <==
	1.2854  five_prime_UTR  126     0       874
	1.2761  three_prime_UTR 92      0       908

#### Conclusion 4: not all FRAG lines may exhibit the same patterns of feature enrichment





# Analysis: part 3


## Count breakpoints that contain features of interest ##

Now that we have a good idea about the pattern of genes and replication origins in 
breakpoint regions, I wanted to more simply ask 'how many breakpoint regions' contain
at least 1 gene or replication origin (simply overlapping by a single bp counts as 
'contained' in this case). Previously we are looking at the total amount of bp that a 
feature occupies in breakpoint regions, but maybe just knowing that 'at least one' 
feature is in a breakpoint region is enough?

My new script can do this for all features at once at various sizes of breakpoint region.
E.g.

	./count_breakpoints_with_features.pl --breakpoint_gff GFF_files/FRAG00062_2x.gff --feature_gff GFF_files/all_TAIR10_features.gff --bp 1000

	FINAL: 36/46 breakpoint regions(78.3%) overlap with CDS
	FINAL: 0/46 breakpoint regions(0.0%) overlap with DNA_replication_origin
	FINAL: 19/46 breakpoint regions(41.3%) overlap with DNAseI_hypersensitive_site
	FINAL: 38/46 breakpoint regions(82.6%) overlap with exon
	FINAL: 13/46 breakpoint regions(28.3%) overlap with five_prime_UTR
	FINAL: 39/46 breakpoint regions(84.8%) overlap with gene
	FINAL: 39/46 breakpoint regions(84.8%) overlap with mRNA
	FINAL: 0/46 breakpoint regions(0.0%) overlap with miRNA
	FINAL: 0/46 breakpoint regions(0.0%) overlap with ncRNA
	FINAL: 16/46 breakpoint regions(34.8%) overlap with open_chromatin_state_1
	FINAL: 15/46 breakpoint regions(32.6%) overlap with open_chromatin_state_2
	FINAL: 12/46 breakpoint regions(26.1%) overlap with open_chromatin_state_3
	FINAL: 9/46 breakpoint regions(19.6%) overlap with open_chromatin_state_4
	FINAL: 6/46 breakpoint regions(13.0%) overlap with open_chromatin_state_5
	FINAL: 14/46 breakpoint regions(30.4%) overlap with open_chromatin_state_6
	FINAL: 6/46 breakpoint regions(13.0%) overlap with open_chromatin_state_7
	FINAL: 6/46 breakpoint regions(13.0%) overlap with open_chromatin_state_8
	FINAL: 0/46 breakpoint regions(0.0%) overlap with open_chromatin_state_9
	FINAL: 37/46 breakpoint regions(80.4%) overlap with protein
	FINAL: 0/46 breakpoint regions(0.0%) overlap with pseudogene
	FINAL: 0/46 breakpoint regions(0.0%) overlap with pseudogenic_exon
	FINAL: 0/46 breakpoint regions(0.0%) overlap with pseudogenic_transcript
	FINAL: 0/46 breakpoint regions(0.0%) overlap with rRNA
	FINAL: 32/46 breakpoint regions(69.6%) overlap with satellite
	FINAL: 0/46 breakpoint regions(0.0%) overlap with snRNA
	FINAL: 0/46 breakpoint regions(0.0%) overlap with snoRNA
	FINAL: 0/46 breakpoint regions(0.0%) overlap with tRNA
	FINAL: 15/46 breakpoint regions(32.6%) overlap with three_prime_UTR
	FINAL: 5/46 breakpoint regions(10.9%) overlap with transposable_element
	FINAL: 0/46 breakpoint regions(0.0%) overlap with transposable_element_gene
	FINAL: 5/46 breakpoint regions(10.9%) overlap with transposon_fragment

I also combined just the gene and replication origins into a new GFF file (to make things 
fast). We are really interested in comparing 2x and 3x breakpoint regions, and testing 
different sizes of breakpoint regions:

	# 1000 bp
	./count_breakpoints_with_features.pl --breakpoint_gff GFF_files/FRAG00062_2x.gff --feature_gff GFF_files/genes_and_origins.gff --bp 1000
	FINAL: 0/46 breakpoint regions (0.0%) overlap with DNA_replication_origin
	FINAL: 39/46 breakpoint regions (84.8%) overlap with gene

	./count_breakpoints_with_features.pl --breakpoint_gff GFF_files/FRAG00062_3x.gff --feature_gff GFF_files/genes_and_origins.gff --bp 1000
	FINAL: 5/44 breakpoint regions (11.4%) overlap with DNA_replication_origin
	FINAL: 33/44 breakpoint regions (75.0%) overlap with gene

	# 10,000 bp
	./count_breakpoints_with_features.pl --breakpoint_gff GFF_files/FRAG00062_2x.gff --feature_gff GFF_files/genes_and_origins.gff --bp 10000
	FINAL: 4/46 breakpoint regions (8.7%) overlap with DNA_replication_origin
	FINAL: 45/46 breakpoint regions (97.8%) overlap with gene

	./count_breakpoints_with_features.pl --breakpoint_gff GFF_files/FRAG00062_3x.gff --feature_gff GFF_files/genes_and_origins.gff --bp 10000
	FINAL: 17/44 breakpoint regions (38.6%) overlap with DNA_replication_origin
	FINAL: 42/44 breakpoint regions (95.5%) overlap with gene

	# 25,000 bp
	./count_breakpoints_with_features.pl --breakpoint_gff GFF_files/FRAG00062_2x.gff --feature_gff GFF_files/genes_and_origins.gff --bp 25000
	FINAL: 10/46 breakpoint regions (21.7%) overlap with DNA_replication_origin
	FINAL: 46/46 breakpoint regions (100.0%) overlap with gene

	./count_breakpoints_with_features.pl --breakpoint_gff GFF_files/FRAG00062_3x.gff --feature_gff GFF_files/genes_and_origins.gff --bp 25000
	FINAL: 20/44 breakpoint regions (45.5%) overlap with DNA_replication_origin
	FINAL: 44/44 breakpoint regions (100.0%) overlap with gene	

#### Conclusion 5: Breakpoint regions of 10 Kbp or larger are almost certain to contain at least one gene (in either 2x or 3x data sets)	

#### Conclusion 6: The majority of 3x breakpoint regions do *not* contain replication origins, however they are more likely to contain them than 2x breakpoint regions


####################################
#
# CHECKING
#
####################################


### Updated results for paper ###

Now want to reproduce results with more sizes of breakpoint regions. Only interested in 
FRAG00062 at first, and will just use percentage of overlapping bp inside/outside breakpoint
regions in the final figure.

	#!/bin/bash

	for i in 100 1000 10000;
	do
        echo "./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062.gff    --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00062_${i}bp_S1000.tsv \& ";
        ./overlap_between_two_gff_files.pl       --breakpoint_gff GFF_files/FRAG00062.gff    --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00062_${i}bp_S1000.tsv \& 

        echo "./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00045.gff    --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00045_2x_${i}bp_S1000.tsv \& ";
        ./overlap_between_two_gff_files.pl       --breakpoint_gff GFF_files/FRAG00045.gff    --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00045_2x_${i}bp_S1000.tsv \& 

        echo "./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062_2x.gff --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00062_2x_${i}bp_S1000.tsv \& ";
        ./overlap_between_two_gff_files.pl       --breakpoint_gff GFF_files/FRAG00062_2x.gff --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00062_2x_${i}bp_S1000.tsv \& 

        echo "./overlap_between_two_gff_files.pl --breakpoint_gff GFF_files/FRAG00062_3x.gff --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00062_3x_${i}bp_S1000.tsv \& ";
        ./overlap_between_two_gff_files.pl       --breakpoint_gff GFF_files/FRAG00062_3x.gff --feature_gff GFF_files/all_TAIR10_features.gff --verbose --bp $i --shuffles 1000 > Results/FRAG00062_3x_${i}bp_S1000.tsv \& 


	done

Combine results into two files (not properly sorted):

	cat FRAG00062_2x_*for_paper* | sort -ru > 2x_for_paper_results_log.tsv
	cat FRAG00062_3x_*for_paper* | sort -ru > 3x_for_paper_results_log.tsv

Then try with modified loop (100 to 5000, in 100 bp increments. Just change first part of loop:

	for i in {100..5000..100}
	
And then combine results as before (but with better sorting):

	cat FRAG00062_2x_*for_paper* | sort -u | sort -k4,4 -k3n,3n > 2x_for_paper_results_100-5000.tsv
	cat FRAG00062_3x_*for_paper* | sort -u | sort -k4,4 -k3n,3n > 3x_for_paper_results_100-5000.tsv



### New script ###

Finally attempted what I should have done all along. A script that uses a sliding window
around the left or right edges of 2x or 3x blocks. Have to run separately for left vs right
vs both, 2x vs 3x, and genes vs replication origins. Can use different ranges (how far out to go
either side of breakpoint), window size, and step sizes. Defaults to +- 50,000 bp (range),
2,500 bp (window size), and 500 bp (step size). 

Needs 18 runs in total for all data:

	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes.gff --v --mode left > figure_2x_genes_L.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes.gff --v --mode right > figure_2x_genes_R.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes.gff --v --mode both > figure_2x_genes_B.tsv &

	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff replication_origins.gff --v --mode left > figure_2x_origins_L.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff replication_origins.gff --v --mode right > figure_2x_origins_R.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff replication_origins.gff --v --mode both > figure_2x_origins_B.tsv &

	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff DHS.gff --v --mode left > figure_2x_DHS_L.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff DHS.gff --v --mode right > figure_2x_DHS_R.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff DHS.gff --v --mode both > figure_2x_DHS_B.tsv &
	
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes.gff --v --mode left > figure_3x_genes_L.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes.gff --v --mode right > figure_3x_genes_R.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes.gff --v --mode both > figure_3x_genes_B.tsv &

	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff replication_origins.gff --v --mode left > figure_3x_origins_L.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff replication_origins.gff --v --mode right > figure_3x_origins_R.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff replication_origins.gff --v --mode both > figure_3x_origins_B.tsv &

	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff DHS.gff --v --mode left > figure_3x_DHS_L.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff DHS.gff --v --mode right > figure_3x_DHS_R.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff DHS.gff --v --mode both > figure_3x_DHS_B.tsv &


Decided to then focus on combining data from both ends of break point and also started to 
look at chromatin state 2 data. The problem with combining L & R sides, is that when
my region of interest is something like -1000 to -900 bp, this is *outside* the block
on the left edge, but *inside* the block on the right edge. 

A solution is to use a new --flip option which inverts the coordinates of the right side.
This means that we can focus on comparing regions that are inside vs outside of each 
block edge. First test to see if this makes a difference:

	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes.gff         > figure_2x_genes_B_3000_100_25.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes.gff  --flip > figure_2x_genes_B_3000_100_25_flipped.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes.gff         > figure_3x_genes_B_3000_100_25.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes.gff  --flip > figure_3x_genes_B_3000_100_25_flipped.tsv &

	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff replication_origins.gff         > figure_2x_origins_B_25000_2000_200.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff replication_origins.gff  --flip > figure_2x_origins_B_25000_2000_200_flipped.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff replication_origins.gff         > figure_3x_origins_B_25000_2000_200.tsv &
	./find_bias_around_breakpoints.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff replication_origins.gff  --flip > figure_3x_origins_B_25000_2000_200_flipped.tsv &


This approach of flipping seems to work and gives a better picture of what is happening
inside and outside duplicated and triplicated blocks. Now to set up some final runs
using suitable ranges, bin sizes, and step factors for each particular genomic feature:


	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat genes.gff --flip --range 25000 --bin 100 --step 25 > figure_2x_genes_25000_100_25.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat genes.gff --flip --range 25000 --bin 100 --step 25 > figure_3x_genes_25000_100_25.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat replication_origins.gff --flip --range 25000 --bin 2000 --step 200 > figure_2x_origins_25000_2000_200.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat replication_origins.gff --flip --range 25000 --bin 2000 --step 200 > figure_3x_origins_25000_2000_200.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat DHS.gff --flip --range 25000 --bin 100 --step 25 > figure_2x_DHS_25000_100_25.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat DHS.gff --flip --range 25000 --bin 100 --step 25 > figure_3x_DHS_25000_100_25.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat state2.gff --flip --range 25000 --bin 100 --step 25 > figure_2x_state2_25000_100_25.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat state2.gff --flip --range 25000 --bin 100 --step 25 > figure_3x_state2_25000_100_25.tsv &

	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat genes.gff --flip  > figure_2x_genes_25000_500_50.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat genes.gff --flip  > figure_3x_genes_25000_500_50.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat replication_origins.gff --flip > figure_2x_origins_25000_500_50.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat replication_origins.gff --flip  > figure_3x_origins_25000_500_50.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat DHS.gff --flip  > figure_2x_DHS_25000_500_50.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat DHS.gff --flip  > figure_3x_DHS_25000_500_50.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat state2.gff --flip > figure_2x_state2_25000_500_50.tsv &
	./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat state2.gff --flip  > figure_3x_state2_25000_500_50.tsv &



We have 5 shorter 2x regions which may (somehow) be more like 3x regions. To test 
this I made a version of FRAG00062_2x.gff which excludes these duplicated blocks:

	cat FRAG00062_2x.gff | grep -vE "block0010|block0017|block0032|block0033|block0035|block0042" > FRAG00062_2x_lite.gff
	
Can try seeing whether there is a difference to the 2x genes plot when you exclude these
5 data points:

	./find_bias_around_breakpoints.pl --break FRAG00062_2x_lite.gff --feat genes.gff --flip --range 3000 --bin 100 --step 25 > figure_2x_lite_genes_25000_100_25.tsv &

This doesn't seem to make too much difference. I.e. short 2x regions are possibly not that
different from the longer regions (at least in terms of gene density around edge of 
blocks).

Now we seem to have settled on some window and step sizes that work, maybe time to finalize
on some final output files. First need to make some new GFF files:

	grep -w transposable_element all_TAIR10_features.gff  > transposable_element.gff
	grep -w pseudogene all_TAIR10_features.gff > pseudogene.gff
	grep "Note=\"state1" all_TAIR10_features.gff > state1.gff
	grep "Note=\"state2" all_TAIR10_features.gff > state2.gff
	grep "Note=\"state3" all_TAIR10_features.gff > state3.gff
	grep "Note=\"state4" all_TAIR10_features.gff > state4.gff
	grep "Note=\"state5" all_TAIR10_features.gff > state5.gff
	grep "Note=\"state6" all_TAIR10_features.gff > state6.gff
	grep "Note=\"state7" all_TAIR10_features.gff > state7.gff
	grep "Note=\"state8" all_TAIR10_features.gff > state8.gff
	grep "Note=\"state9" all_TAIR10_features.gff > state9.gff
	grep -w satellite all_TAIR10_features.gff > satellite.gff	
	
Will plot at 3 different scales for various features (many of these might only 
appear as supplemental results):

+ 2,000 bp (window) and 200 bp (step)
+ 500 bp (window) and 50 bp (step)
+ 100 bp (window) and 10 bp (step)

All step sizes are now one tenth of window size. Will run all three combinations for
each feature (using range of -50,000 to +50,000 bp). Time to use a bash looping script
again.

	#!/bin/bash

	for feature in genes pseudogene transposable_element satellite replication_origins DHS state1 state2 state3 state4 state5 state6 state7 state8 state9
	do

			for window in 2000 500 100
			do
			let "step=$window/10"
			echo "Running: $feature $window $step"
			echo "./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat ${feature}.gff --flip  --bin $window --step $step > figure_2x_${feature}_${window}_${step}.tsv &"
			./find_bias_around_breakpoints.pl --break FRAG00062_2x.gff --feat ${feature}.gff --flip  --bin $window --step $step > figure_2x_${feature}_${window}_${step}.tsv &

			echo "./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat ${feature}.gff --flip  --bin $window --step $step > figure_3x_${feature}_${window}_${step}.tsv &"
			./find_bias_around_breakpoints.pl --break FRAG00062_3x.gff --feat ${feature}.gff --flip  --bin $window --step $step > figure_3x_${feature}_${window}_${step}.tsv &
			echo
			done
	done