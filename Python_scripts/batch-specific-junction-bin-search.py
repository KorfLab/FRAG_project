#! /usr/bin/env python2.6

import os, sys, math, time
import subprocess
from optparse import OptionParser
from collections import defaultdict

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-b", "--startbinfile", dest="b", help="Input bin file list.")
parser.add_option("-Q", "--fastqmode", dest="fqmode", action="store_true", default = False, help="Output a fastq instead of fasta.")

(opt, args) = parser.parse_args()

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def status(cur, total):
      now = time.time()-start
      perdone = cur/float(total-1)
      percent = (perdone)*100
      emin = now/60
      esec = now%60
      rmin = (total-cur)/(cur/now)/60
      rsec = (total-cur)/(cur/now)%60
      sys.stdout.write('\r')      
      sys.stdout.write(">%-35s< %i%% Elapsed: %im %is Remaining: %im %is     " % ('='*int((perdone*35)),percent,emin,esec,rmin,rsec))
      sys.stdout.flush()

def strtolen(s1):   
   ns = ''
   ls = []
   for c in s1:
      if c.isdigit() == False:
         ns+='.'
         ls.append(c)
      else:
         ns+=c
   ns = ns.split('.')
   ns.remove('')  
   ns = map(lambda x: int(x), ns)
   #print ns, ls
   leng = 0
   for i in range(len(ns)):
      if ls[i] == 'S' or ls[i] == 'D':
         leng -= 0
      else:
         leng += ns[i]
   return leng

def fixseq(ch1, s1, q1):
   ns = ''
   ls = []
   seq = s1[:]
   qual = q1[:]
   for c in ch1:
      if c.isdigit() == False:
         ns+='.'
         ls.append(c)
      else:
         ns+=c
   
   ns = ns.split('.')
   ns.remove('')  
   ns = map(lambda x: int(x), ns)
   fl = zip(ls, ns)
   
   if fl[0][0] == 'S':
      seq = seq[fl[0][1]:]
      qual = qual[fl[0][1]:]
   
   
   if fl[-1][0] == 'S':
      seq = seq[:-fl[-1][1]]
      qual = qual[:-fl[-1][1]]
   
   return seq, qual

###############################   START   #############################


bins = open(opt.b)
#get rid of header
bins.readline()
targets = {}
validchroms = []
fhandle = {}
for l in bins:
   x = l[:-1].split('\t')
   if opt.fqmode == False:
      send = ".fa"
   else:
      send = ".fq"
   fname = "junct-"+x[0]+'-'+x[1]+'_'+x[2]+send
   try:
      targets[x[0]].append([int(x[1]), int(x[2]), fname])
      fhandle[fname] = {}
   except:
      targets[x[0]] = [[ int(x[1]), int(x[2]), fname]]
      validchroms.append(x[0])
      fhandle[fname] = {}

all = os.listdir(os.getcwd())
fqs = filter(lambda x: ".fq" in x or ".fastq" in x, all)
fqs.sort()
print fqs
for fqfile in fqs:
   print '\n'+fqfile
   samname = fqfile.split('.')[0]+'_aln.sam'
   if samname not in all:
      print "WARNING: "+fqfile+" does not have a matching .sam labeled "+samname
      continue
   fqc = file_len(fqfile)/8
   f = open(fqfile)
   s = open(samname)
   
   goodreads = {}
   goodreads[1] = defaultdict(str)
   goodreads[2] = defaultdict(str)

   for fhand in fhandle:
      fhandle[fhand][fqfile] = open(fqfile.split('.')[0]+'-'+fhand, 'w')

   print "Reading in "+ fqfile
   ct2 = 1
   start = time.time()
   while 1:
      if ct2 % 50000 == 8:
         status(ct2, fqc)
      ct2+=1
      
      name = f.readline()
      if name == '':
         break
      seq = f.readline()
      plus = f.readline()
      qual = f.readline()
    
      searchname = name[1:-1].split(' ')[0]
      
      if opt.fqmode == False:
         goodreads[1][searchname] = ('>'+name[1:]+seq)#+plus+qual)
      else:
         goodreads[1][searchname] = (name+seq+plus+qual)
         
      name2 = f.readline()
      seq2 = f.readline()
      plus2 = f.readline()
      qual2 = f.readline() 
      
      if opt.fqmode == False:
         goodreads[2][searchname] = ('>'+name2[1:]+seq2)#+plus2+qual2)
      else:
         goodreads[2][searchname] = (name2+seq2+plus2+qual2)

   ########### reading samfile
   ct = 1
   smc = file_len(samname)
   print "\nReading in "+ samname
   start = time.time()
   while True:
      if ct % 10000 == 8:
         status(ct, smc/2)
      ct+=1
   
      l = s.readline()
      if l == '':
         break
      if l[0] == '@':
         smc -= 1
         continue
      l2 = s.readline()
   
      x1 = l[:-1].split('\t')
      chrom1 = x1[2]
      spos1 = int(x1[3])
      epos1 = int(x1[3])+len(x1[9])+1
      
      x2 = l2[:-1].split('\t')
      chrom2 = x2[2]
      spos2 = int(x2[3])
      epos2 = int(x2[3])+len(x2[9])+1
      
      if chrom1 not in validchroms and chrom2 not in validchroms:
         continue      

      c1 = x1[-1].count(';')+1
      c2 = x2[-1].count(';')+1
      tc = c1 * c2
         
      one = []
      one.append(x1)
      #if opt.seq in x1[-1]:
      aline = x1[-1][5:].split(';')
      for alt in aline:
         val = alt.split(',')
         if val[0] not in validchroms:
            continue
         if 'S' in val[2]:
            break
         #print x1
         nx1 = x1[:]
         nx1[2] = val[0]
         if int(val[1]) < 0:
            nx1[1] = '16'
         else:
            nx1[1] = '0'
         nx1[3] = val[1][1:]
         nx1[5] = val[2]       
         nx1[-1] = nx1[-1].replace(alt+';','')         
         if nx1[-1] == 'XA:Z:':
            nx1 = nx1[:-1]
         if 'XC:i:' in l and 'S' not in val[2]:
            nx1 = nx1[:11]+nx1[12:]
         if len(x1[9]) !=  strtolen(val[2]):       
            nx1[9], nx1[10] = fixseq(x1[5], x1[9], x1[10])  
         nx1[12] = "NM:i:"+val[3]

         if 'S' not in val[2]:
            one.append(nx1)             
      
      two = []
      two.append(x2)
#      if opt.seq != x2[2] and opt.seq in l2:
      aline = x2[-1][5:].split(';')
      for alt in aline:
         val = alt.split(',')
         if val[0] not in validchroms:
            continue
         if 'S' in val[2]:
            break
         #print x2
         nx2 = x2[:]
         nx2[2] = val[0]
         if int(val[1]) < 0:
            nx2[1] = '16'
         else:
            nx2[1] = '0'
         nx2[3] = val[1][1:]
         nx2[5] = val[2]
         nx2[-1] = nx2[-1].replace(alt+';','')
         if nx2[-1] == 'XA:Z:':
            nx2 = nx2[:-1]
         if 'XC:i:' in l2 and 'S' not in val[2]:
            nx2 = nx2[:11]+nx2[12:]
         if len(x2[9]) !=  strtolen(val[2]):       
            nx2[9], nx2[10] = fixseq(x2[5], x2[9], x2[10])               
         nx2[12] = "NM:i:"+val[3]
         if 'S' not in val[2]:
            two.append(nx2)             
      bins = []
      count = 1  
      for line1 in one:    
         for line2 in two:
            #print line1, line2
   
            chrom1 = line1[2]
            spos1 = int(line1[3])
            epos1 = int(line1[3])+len(line1[9])+1
            
            chrom2 = line2[2]
            spos2 = int(line2[3])
            epos2 = int(line2[3])+len(line2[9])+1
            
            if chrom1 not in validchroms and chrom2 not in validchroms:
               continue

            if chrom1 == chrom2: 
               for bin in targets[chrom1]:
                  if spos1 <= bin[1] and epos1 >= bin[0] and spos2 <= bin[1] and epos2 >= bin[0]:
                     continue
                  elif (spos1 <= bin[1] and epos1 >= bin[0]) or (spos2 <= bin[1] and epos2 >= bin[0]):                             
                     bins.append(bin[-1])                       
            else:
               if chrom1 in validchroms:
                  for bin in targets[chrom1]:
                     if spos1 <= bin[1] and epos1 >= bin[0]:
                        bins.append(bin[-1])                 
                  
               if chrom2 in validchroms:
                  for bin in targets[chrom2]:
                     if spos1 <= bin[1] and epos1 >= bin[0]:
                        bins.append(bin[-1]) 
                 
      bins = list(set(bins))
      for fh in bins:
         fhandle[fh][fqfile].write(goodreads[1][x1[0]])
         fhandle[fh][fqfile].write(goodreads[2][x2[0]]) 
       

      
   s.close()
   f.close()

for junct in fhandle:
   for ofh in fhandle[junct]:
      fhandle[junct][ofh].close()      
      
print '\n'
      
      
      
      
      
      
 
 
 
   

      
   
      