#! /usr/bin/env python

# Chanlab, UC Davis
# Han Tan 2012 

# This script filters out filtered reads from Illumina that was retained
# after demultiplexing. Also checks the integrity of paired reads.

# USAGE: python Programname.py [m] [-p]
# [m] enter minimum length after trimming 
# [-p] option if running paired ends

import os, sys, math, gzip

minlength = int(sys.argv[1])

pair = 0
pairs = 0
total = 0
filtered = 0
noN = 0

#Adaptors to look for

solAd1 = "AGATCGGAAG"
solB = "AGATCGGAAGAGC"

if '-p' in sys.argv:
    pair = 1 

#prepare Raw_data directory
li = os.listdir(os.getcwd())
dirs = ['Raw_data']
print li
if False in list((x in li for x in dirs)):
   os.system("mkdir Raw_data")

#os.system("mkdir Raw_data")

li = os.listdir(os.getcwd())
todo = filter(lambda x: x[-9:] == '.fastq.gz', li)
todo.sort()
for file in todo:
    print file
    name = file.split('_')[0]
    f = gzip.open(file, 'rb')
    o = open(name+'.fq', 'w')
    total = 0
    filtered = 0
    noN = 0
    pairs = 0

    if pair == 1:
        print 'Running Paired-End Reads'
        while True:
# Open file and read 8 lines at a time 
            name1 = f.readline()
            if name1 == '':
                break
            seq1 = f.readline()
            plus1 = f.readline()
            qual1 = f.readline()

            name2 = f.readline()
            if name2 == '':
                break
            seq2 = f.readline()
            plus2 = f.readline()
            qual2 = f.readline()

            l = name1.split()
            read1 = l[0]
            filter1 = l[1][2]

            m = name2.split()
            read2 = m[0]
            filter2 = m[1][2]

            total += 1
# Remove Cassava 1.8 filtered ie 'Y' reads
            if read1 == read2 and filter1 == 'N' and filter2 == 'N':
                filtered += 1
# Remove sequences with N's
                if 'N' not in seq1 and 'N' not in seq2:
                    noN += 1
# Remove Primary Adaptor Contamination
                    if solB in seq1:
                        t1 = seq1[:seq1.index(solB)]
                        l1 = len(t1)
                        seq1 = t1+'\n'
                        qual1 = qual1[:l1]+'\n'
                        if solB in seq2:
                            t2 = seq2[:seq2.index(solB)]
                            l2 = len(t2)
                            seq2 = t2+'\n'
                            qual2 = qual2[:l2]+'\n'
# Remove Secondary Adaptor Contamination
                    if solAd1 in seq1:
                        t1 = seq1[:seq1.index(solAd1)]
                        l1 = len(t1)
                        seq1 = t1+'\n'
                        qual1 = qual1[:l1]+'\n'
                    if solAd1 in seq2:
                        t2 = seq2[:seq2.index(solAd1)]
                        l2 = len(t2)
                        seq2 = t2+'\n'
                        qual2 = qual2[:l2]+'\n'
# Trim quality
                    quals = map(lambda x: ord(x)-33, qual1[:-1])
                    for x in range(len(quals)-4):
                        cut = quals[x:x+5]
                        ave = float(sum(cut))/5.0
                        if ave < 20.0:
                            qual1 = qual1[:x]+'\n'
                            seq1 = seq1[:x]+'\n'
                            break
                    quals2 = map(lambda x: ord(x)-33, qual2[:-1])
                    for x in range(len(quals2)-4):
                        cut2 = quals2[x:x+5]
                        ave2 = float(sum(cut2))/5.0
                        if ave2 < 20.0:
                            qual2 = qual2[:x]+'\n'
                            seq2 = seq2[:x]+'\n'
                            break
# Check minimum length
                    if len(seq1) > minlength and len(seq2) > minlength and len(seq1) == len(qual1) and len(seq2) == len(qual2):
                        o.write(name1+seq1+plus1+qual1+name2+seq2+plus2+qual2)
                        pairs += 1

    else:
        print "Running Single Reads"
        while True:
# Open file and read 4 lines at a time
            name1 = f.readline()
            if name1 == '':
                break
            seq1 = f.readline()
            plus1 = f.readline()
            qual1 = f.readline()

            l = name1.split()
            filter1 = l[1][2]
            
            total +=1
# Remove Cassava 1.8 filtered ie 'Y' reads                                                                                                                      
            if filter1 == 'N':
                filtered += 1
#Remove sequences with N's                                                                                                                                      
                if 'N' not in seq1:
                    noN += 1
# Remove Primary Adaptor Contamination                                                                                                                          
                    if solB in seq1:
                        t1 = seq1[:seq1.index(solB)]
                        l1 = len(t1)
                        seq1 = t1+'\n'
                        qual1 = qual1[:l1]+'\n'
# Remove Secondary Adaptor Contamination                                                                                                                        
                    if solAd1 in seq1:
                        t1 = seq1[:seq1.index(solAd1)]
                        l1 = len(t1)
                        seq1 = t1+'\n'
                        qual1 = qual1[:l1]+'\n'
# Trim quality                                                                                                                                                  
                    quals = map(lambda x: ord(x)-33, qual1[:-1])
                    for x in range(len(quals)-4):
                        cut = quals[x:x+5]
                        ave = float(sum(cut))/5.0
                        if ave < 20.0:
                            qual1 = qual1[:x]+'\n'
                            seq1 = seq1[:x]+'\n'
                            break
# Check minimum length                                                                                                                                          
                    if len(seq1) > minlength:
                        o.write(name1+seq1+plus1+qual1)
                        pairs += 1

    f.close()
    o.close()

    print 'Total:'
    print total
    print 'Filtered:'
    print filtered
    print 'no Ns:'
    print noN
    print 'Good Reads or Pairs:'
    print pairs

    os.system("mv "+file+" Raw_data/")

