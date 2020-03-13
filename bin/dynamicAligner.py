#!/usr/bin/python

#script.py --seqs seqFil --size 300 --smallAln clustalo_msa,mafft_msa --largeAln clustalo_msa,mafft_msa

import sys, getopt, subprocess, os

def countSequences(fileName):
    num = len([1 for line in open(fileName) if line.startswith(">")])
    return num

def convert(string): 
    li = list(string.split(",")) 
    return li 

def countGaps(alnFileName):
    fh = open(alnFileName)
    auxGap = 0
    globalGap = 0
    gap = '-'
    for line in fh:
        if not line.startswith(">"):
            auxGap = line.count(gap)
            globalGap += auxGap
        
    fh.close()
    return globalGap

def main(argv):
   seqs = ''
   size = ''
   smallAln = ''
   largeAln = ''
   numSeqs=''
   try:
      opts, args = getopt.getopt(argv,"hs:",["seqs=","size=","smallAln=","largeAln="])
   except getopt.GetoptError:
      print 'test.py --seqs <seqsFile> --size <sizeInt> --smallAln <aln1,aln2> --largeAln <aln1,aln2>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
        print 'test.py --seqs <seqsFile> --size <sizeInt> --smallAln <aln1,aln2> --largeAln <aln1,aln2>'
        sys.exit()
      elif opt in ("-s","--seqs"):
        seqs = arg
      elif opt in ("--size"):
        size = int(arg)
      elif opt in ("--smallAln"):
        smallAln = arg
      elif opt in ("--largeAln"):
        largeAln = arg

   numSeqs=countSequences(seqs)
   print 'Input file is ', seqs
   print 'Input file has ', numSeqs,' sequences'
   print 'The size is ', size
   print 'The list of small aligners are ', smallAln   
   print 'The list of large aligners are ', largeAln

   if(numSeqs<size):
        dictSmall=dict()

        print 'MODE Small aligners'
        smallAln=convert(smallAln)
        for aln in smallAln:
           cmd = "t_coffee -seq %s -method %s -outfile %s_%s.aln -quiet"%(seqs,aln,seqs,aln)
           #print 'RUNING : ', cmd
           os.system(cmd)

           gapFile = "%s_%s.aln"%(seqs,aln)
           numGaps = countGaps(gapFile)
           #print '**Num gaps : ', numGaps,' aln file: ',gapFile
           dictSmall[aln]=numGaps
       
        minAligner = min(dictSmall.keys(), key=(lambda k: dictSmall[k]))
        minGap=dictSmall[minAligner]
        print '***\n** The minimal alignment has %d gaps and it\'s from %s aligner **\n***'%(minGap,minAligner)
        print 'The complete results are: ',dictSmall

   elif(numSeqs>size):
       dictLarge=dict()

       print 'MODE Large aligners'
       largeAln=convert(largeAln)
       for aln in largeAln:
          cmd = "t_coffee -seq %s -method %s -outfile %s_%s.aln -quiet"%(seqs,aln,seqs,aln)
          #print 'RUNING : ', cmd
          os.system(cmd)

          gapFile = "%s_%s.aln"%(seqs,aln)
          numGaps = countGaps(gapFile)
          #print '**Num gaps : ', numGaps,' aln file: ',gapFile
          dictLarge[aln]=numGaps
       
       minAligner = min(dictLarge.keys(), key=(lambda k: dictLarge[k]))
       minGap=dictLarge[minAligner]
       print '***\n** The minimal alignment has %d gaps and it\'s from %s aligner **\n***'%(minGap,minAligner)
       print 'The complete results are: ',dictLarge

if __name__ == "__main__":
   main(sys.argv[1:])

