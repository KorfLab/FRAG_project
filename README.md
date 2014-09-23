FRAG_project
============

To contain code and data files for the Fragmentation project, investigating shattered 
chromosome phenotypes in _Arabidopsis thaliana_.
	
	
## Background

Raw data from Han's mapping of reads back to reference genome is in the file:

	CGRs_1kb_plots_forkeith.tsv

This file contains raw and normalized counts of reads mapped to 1 Kbp bins in three lines 
(FRAG62, FRAG133, and FRAG80 - a control line). The normalized counts don't just use 
FRAG80 but are based on a concatenated file of all the F1 diploids that were sequenced. 
Han calls this the 'SuperF1diploid'.

The normalized counts come out to approximately 2.0 for duplications, and 3.0 for 
triplications. Actual values are slightly lower than this.



### Making Master GFF file of all genomic features

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




### Checking gene orientation ###

At this point we realized that we would like to know whether the enrichment of genes 
inside junction regions followed any pattern. I.e. are there more likely to be convergently
transcribed genes (hence more 3' UTRs) than divergently or tandemly transcribed genes.

I put only gene features from all_TAIR10_feature.gff into a new file (genes.gff) and made
a new script to test this:

	./check_gene_orientation.pl --breakpoint_gff FRAG00062.gff --feature_gff genes.gff --target gene
	
This revealed that across the entire genome, (protein-coding) gene orientation is 
effectively random (as expected):

	>>      6987    %25.82
	<>      6508    %24.05
	<<      7058    %26.08
	><      6508    %24.05

Running this script with newer breakpoint data for FRAG00062 (now with 90 breakpoints):

	./check_gene_orientation.pl --breakpoint_gff FRAG00062.gff --feature_gff genes.gff --target gene

	60130895 coding bp in genome out of 119146348 bp (50.47%)
	>>>|>>>	25	%27.78
	<<<|<<<	31	%34.44
	>>>---|---<<<	4	%4.44
	<<<---|---<<<	8	%8.89
	<<<---|--->>>	14	%15.56
	>>>---|--->>>	8	%8.89

Breakpoints inside genes account for 62% of all breakpoints (56/90) and those that are inside
divergently transcribed genes account for 41% of all intergenic breakpoints (14/34).

And for FRAG00045 (30 breakpoints):

	./check_gene_orientation.pl --breakpoint_gff FRAG00045.gff --feature_gff genes.gff --target gene

	60130895 coding bp in genome out of 119146348 bp (50.47%)
	>>>|>>>	12	%40.00
	<<<|<<<	5	%16.67
	>>>---|---<<<	3	%10.00
	<<<---|---<<<	3	%10.00
	>>>---|--->>>	5	%16.67
	<<<---|--->>>	2	%6.67

This time only 57% of breakpoints were inside genes (17/30) and only 15% of intergenic 
breakpoints are between divergently transcribed genes. So no strong pattern.




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




## Analysis of enriched features in breakpoint regions ##

If we define a 'breakpoint region' as a window of sequence around each breakpoint location
(mapped to the reference genome), we can ask whether any genomic features are enriched
in these breakpoint regions. E.g. take 1,000 bp around all breakpoints (potentially
overlapping other breakpoints) and ask whether the total bp of a feature such as 'coding
exons' is higher (as a percentage) *inside* those regions vs all DNA *outside* those
regions.

Can try this for many different sizes of breakpoint (100 bp up to 50 Kbp). In this analysis
certain GFF features (chromosome, rRNA, snRNA) are ignored.

The final result is calculated as a ratio of %breakpoint-region-occupied-by-feature compared
to %non-breakpoint-region-occupied-by-feature. To assess the significance of these ratios,
I perform shuffling experiments to see whether we see similar ratios when we randomize the 
location of all of the breakpoints (for the tailswap region we allow the possibility of 
all junctions  occurring in Chr1, or Chr4, or any combination of both).


#### FRAG00062 ####

For this dataset we can also separate out the breakpoint data into two subsets,
those breakpoints that flank either 2x or 3x blocks:

	grep duplicated FRAG00062.gff   | sed 's/.*ID=\(block[0-9]*\);N.*/\1/' > duplicated_blocks.txt
	grep triplicated FRAG00062.gff  | sed 's/.*ID=\(block[0-9]*\);N.*/\1/' > triplicated_blocks.txt
	grep -f duplicated_blocks.txt FRAG00062.gff > FRAG00062_2x.gff
	grep -f triplicated_blocks.txt FRAG00062.gff > FRAG00062_3x.gff

<<<<<<< HEAD
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




=======
This means that each 'run' of  analysis for any given size of breakpoint region, occurs 
3 times:

	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062.gff    --feature_gff all_TAIR10_features.gff --shuffles 1000 --verbose --bp 10000 > FRAG00062_breakpoints_s1000_L10000.tsv
	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff all_TAIR10_features.gff --shuffles 1000 --verbose --bp 10000 > FRAG00062_2x_breakpoints_s1000_L10000.tsv
	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff all_TAIR10_features.gff --shuffles 1000 --verbose --bp 10000 > FRAG00062_3x_breakpoints_s1000_L10000.tsv

Note that only the first row of the output file (row 0) contains data from the unshuffled
results. All other results rows reflect shuffled data. The full set of columns of output are:

1. Run number (starts at 0 for unshuffled results)
2. Real ratio from unshuffled data (the ratio of columns 6 & 7)
3. Feature Breakpoint_region_bp    
4. Non_breakpoint_region_bp        
5. Feature_bp_inside      
6. %Inside 
7. Feature_bp_outside      
8. %Outside        
8. Shuffled_ratio  
9. Above (number of times column 8 exceeds column 2)
10. Same (number of times column 8 equals column 2)   
11. Below (number of times column 8 is below column 2)


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


### New script to just count breakpoints that contain features of interest ###

Now that we have a good idea about the pattern of genes and replication origins in 
breakpoint regions, I wanted to more simply ask 'how many breakpoint regions' contain
at least 1 gene or replication origin (simply overlapping by a single bp counts as 
'contained' in this case). Previously we are looking at the total amount of bp that a 
feature occupies in breakpoint regions, but maybe just knowing that 'at least one' 
feature is in a breakpoint region is enough?

I wrote a script to just do this. First I combined just the gene and replication origins
into a new GFF file (to make things fast). We are really interested in comparing 2x and 3x
breakpoint regions, and testing different sizes of breakpoint regions:

	# 1000 bp
	./count_breakpoints_with_features.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes_and_origins.gff --bp 1000
	FINAL: Processed 46 breakpoint regions, found 39 (84.8%) overlap with gene
	FINAL: Processed 46 breakpoint regions, found 0 (0.0%) overlap with DNA_replication_origin

	./count_breakpoints_with_features.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes_and_origins.gff --bp 1000
	FINAL: Processed 44 breakpoint regions, found 33 (75.0%) overlap with gene
	FINAL: Processed 44 breakpoint regions, found 5 (11.4%) overlap with DNA_replication_origin

	# 10,000 bp
	./count_breakpoints_with_features.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes_and_origins.gff --bp 10000
	FINAL: Processed 46 breakpoint regions, found 45 (97.8%) overlap with gene
	FINAL: Processed 46 breakpoint regions, found 4 (8.7%) overlap with DNA_replication_origin

	./count_breakpoints_with_features.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes_and_origins.gff --bp 10000
	FINAL: Processed 44 breakpoint regions, found 42 (95.5%) overlap with gene
	FINAL: Processed 44 breakpoint regions, found 17 (38.6%) overlap with DNA_replication_origin

	# 25,000 bp
	./count_breakpoints_with_features.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes_and_origins.gff --bp 25000
	FINAL: Processed 46 breakpoint regions, found 46 (100.0%) overlap with gene
	FINAL: Processed 46 breakpoint regions, found 10 (21.7%) overlap with DNA_replication_origin

	./count_breakpoints_with_features.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes_and_origins.gff --bp 25000
	FINAL: Processed 44 breakpoint regions, found 44 (100.0%) overlap with gene
	FINAL: Processed 44 breakpoint regions, found 20 (45.5%) overlap with DNA_replication_origin
	
#### Conclusion 5: Breakpoint regions of 10 Kbp or larger are almost certain to contain at least one gene (in either 2x or 3x data sets)	

#### Conclusion 6: The majority of 3x breakpoint regions do *not* contain replication origins, however they are more likely to contain them than 2x breakpoint regions



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
	
Now combine with main file:

	cat DHS.gff >> all_TAIR10_features.gff
	grep -vE "^#" all_TAIR10_features.gff | sort -k 1,1 -k 4n,4n  all_TAIR10_features.gff > tmp.gff; mv tmp.gff all_TAIR10_features.gff	

Now compare to our master GFF files of breakpoints. Again we may want to try several breakpoint 
region sizes (and do this for 2x and 3x). Will use a wrapper script that will look at 
FRAG00045 and FRAG00062 (2x and 3x) for various sizes. Commands will again be in the form
of:

	./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff DHS.gff --shuffles 1000 --verbose --bp 10000 > FRAG00062_2x_breakpoints_s1000_L10000.tsv


### Now look at Chromatin state information ###

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

Now try a final run with a new set of (simplified) breakpoint region sizes (100, 500,
1 Kbp, 5 Kbp, 10 Kbp, and 20 Kbp). Final set of features that we are looking at is:

CDS
chromosome
DNA_replication_origin
DNAseI_hypersensitive_site
exon
five_prime_UTR
gene
miRNA
mRNA
ncRNA
open_chromatin_state (9 different states)
protein
pseudogene
pseudogenic_exon
pseudogenic_transcript
rRNA
satellite
snoRNA
snRNA
three_prime_UTR
transposable_element
transposable_element_gene
transposon_fragment
tRNA


### Updated results for paper ###

Now want to reproduce results with more sizes of breakpoint regions. Only interested in 
FRAG00062 at first, and will just use percentage of overlapping bp inside/outside breakpoint
regions in the final figure.

Modified overlap_between_two_gff_files.pl script at this point to also include the value
of --bp in the final output file (now as 3rd column of output). Use bash loop script to do
this. First for all breakpoint regions in increasing order of magnitude:

#!/bin/bash

for i in 10 100 1000 10000 100000 1000000;
do
        echo "./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes_and_origins.gff --verbose --bp $i > FRAG00062_2x_${i}bp_for_paper.tsv \& ";
        ./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_2x.gff --feature_gff genes_and_origins.gff --verbose --bp $i > FRAG00062_2x_${i}bp_for_paper.tsv \& ;
        echo "./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes_and_origins.gff --verbose --bp $i > FRAG00062_3x_${i}bp_for_paper.tsv \& ";
        ./overlap_between_two_gff_files.pl --breakpoint_gff FRAG00062_3x.gff --feature_gff genes_and_origins.gff --verbose --bp $i > FRAG00062_3x_${i}bp_for_paper.tsv \& ;
done

Combine results into two files (not properly sorted):

	cat FRAG00062_2x_*for_paper* | sort -ru > 2x_for_paper_results_log.tsv
	cat FRAG00062_3x_*for_paper* | sort -ru > 3x_for_paper_results_log.tsv

Then try with modified loop (100 to 5000, in 100 bp increments. Just change first part of loop:

	for i in {100..5000..100}
	
And then combine results as before (but with better sorting):

	cat FRAG00062_2x_*for_paper* | sort -u | sort -k4,4 -k3n,3n > 2x_for_paper_results_100-5000.tsv
	
