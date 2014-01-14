FRAG_project
============

To contain code and data files for the Fragmentation project, investigating 
shattered chromosome phenotypes in _Arabidopsis thaliana_.


## Data files ##

	CGRs_1kb_plots_forkeith.tsv

This file contains raw and normalized counts of reads mapped to 1 Kbp bins in three lines (FRAG62, FRAG133, and FRAG80 - a control line).
The normalized counts don't just use FRAG80 but are based on a concatenated file of all the F1 diploids that were sequenced. Han calls this the
'SuperF1diploid'.

Earlier work used 2 Kbp bins and there are other files with similar information for th e larger bin sizes but which only normalized using
FRAG80.

The normalized counts come out to approximately 2.0 for duplications, and 3.0 for triplications. Actual values are slightly lower than this.

