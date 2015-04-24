To extract fq reads from a list of chromosome positions using the corresponding SAM file using batch-specific-junction-bin-search.py

#USAGE: batch-specific-junction-bin-search.py -b [binfile.txt] -Q

Options:
  -h, --help            show this help message and exit
  -b B, --startbinfile=B
                        Input bin file list.
  -Q, --fastqmode       Output a fastq instead of fasta.
  
The script to handle multiple putative breakpoint junctions across multiple files. 

It works on a directory of .fq and .sam files.

A couple of points about how to run:

1. It takes a bin.txt file, that has the bin locations in the form "Seq\tBinStart\tBinEnd" with a header line. I have included a sample bin file in the attached folder. Make sure there isn't an extra blank line after the last bin.

2. The only parameters required is the bin input bin file (with a -b), and a -Q if you want the output files to be fastq instead of the default fasta. 

3. In the directory must be a ".fq" and "_aln.sam" file for every library. So if you have a libX1.fq, you must also have a libX1_aln.sam.

4. There will be a file created for every bin and file combination, even if there are no reads for that particular bin. So if there is 2 libraries with three bins in the bin file, there will six result fasta/fastq files even if all are empty.

5. The chromosome/refseq name must exactly match the sam header format. (ie if the refseq was Chr01, chr01 won't work)

7. The result files will automatically have the naming convention libname-junct-[refseq/chrom]-binstart_binend.[fa/fq]

8. There can be multiple bins in the same chromosome/refseq, as long as they don't overlap.
