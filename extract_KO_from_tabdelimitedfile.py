#!/usr/bin/python                                                                                                                                                       
# Author: Sathish <subramanians@wusm.wustl.edu                                                                                                                               

"""                                                                                                                                                                          
This is a python script that takes as input a tab-delimited file and prints out all lines with the KO of interest                                                                                                 
"""

import sys,os
import glob

## method to extract information from gzipped file
def gziplines(fname):
  from subprocess import Popen, PIPE
  f = Popen(['zcat', fname], stdout=PIPE)
  for line in f.stdout:
    yield line

## method to get full species name from metadata
def meta_lookup(key):
  a = open('/home/comp/jglab/sathish/Bangladesh/128_Genomes_Metadata.txt',"r")
  meta = a.readlines()
  for line in meta:
    if key in line:
      (Genome_ID, Genome_Name, Genus_Species, Genome_Length) = line.split("\t")
      return Genome_Name

K0=sys.argv[2]
print "\nReads extracted that map to "+K0
  
SampleList=[]
SampleList=['Sample']
GutHit_MasterList=[]
GutHitCounts_Master={}
path = sys.argv[1]
count_file = 0

for root, dirs, files in os.walk(path):
  for name in files:
    if 'FastaFlatFile.txt.gz' in name and count_file<8:
      # initializing FastaFlatFile and counters
      sample_id = name.strip('_FastaFlatFile.txt.gz')
      print "Sample is "+sample_id
      SampleList.append(sample_id)
      count_reads = 0
      infile = path+sample_id+"/"+name
      print ("Found %s" % name)
      count_file = count_file + 1
      print "Read file number "+str(count_file)
  
      # initializing dictionaries and lists for each FlatFile
      GutHitList=[] # List of all hits 
      GutHitCounts={}  # dictionary with key being genome hit and value being count
      GutHitCounts_unique={}
      GutHitCounts_unique = {'More_than_one_best_hit':0}
      GutHitCounts_fair={}

      # Read fasta flat file and count number of gut Hits for a given K0           
      for line in gziplines(infile):
        (Header, SequenceLength, TooShort, TooManyN, Replicate, HostGenomeHit, PassFilter, KeggBlastHit, KeggAnnotatedHit, MeropsHit, CogBlastHit, CogAnnotatedHit, GutGenomeHit) = line.split("\t")
        if PassFilter=="1": 
          count_reads = count_reads+1
          if K0 in line:
            gutGenomeHit=GutGenomeHit.strip("\n")
            GutHitList.append(gutGenomeHit)
            if gutGenomeHit not in set(GutHit_MasterList):   # build a masterlist of Gutgenomes encountered during the search of KOs iterated over each sample
              GutHit_MasterList.append(gutGenomeHit)
      #counting occurences of a GutHit in one individual
      for i in GutHitList:
        print i
        GutHitCounts[i] = GutHitList.count(i)
        if " " not in i:
          GutHitCounts_unique[i]=GutHitCounts[i]
          GutHitCounts_fair[i]=GutHitCounts[i]

      print "No. of reads pass filter "+str(count_reads)+"\n"

      # Calculate unique and fair dictionaries that either bin reads with more than one best hit into one category or in the fair scheme allocate 0.5 abundance to each of the best hits - change this later to divide by the number of best blast hits

      for key in GutHitCounts.keys():
        if " " in key:
          GutHitCounts_unique['More_than_one_best_hit']=GutHitCounts_unique['More_than_one_best_hit']+GutHitCounts[key]
          Derep=key.split()
          for i in set(Derep):
            if i in GutHitCounts_fair:
              GutHitCounts_fair[i]=(GutHitCounts_fair[i])+(0.5*GutHitCounts[key])  
            else:
              GutHitCounts_fair[i]=(0.5*GutHitCounts[key])
##
      print GutHitCounts
#      GutHitCounts_Master[i]=

print SampleList
print GutHit_MasterList

#for i in (GutHit_MasterList):
#  GutHitCounts_Master[i]=GutHitCounts[key]

#print GutHitCounts_Master
 
#      outfile = open("GutDB_Hits_"+K0+"_"+sample_id+".txt","w")
      #outfile.write( "\t".join(SampleList))
#      outfile.write("Sample_ID"+"\t"+sample_id)
#      outfile.write("\n")
#      for key in GutHitCounts.keys():
#        if meta_lookup(key):
#          outfile.write(str(meta_lookup(key))+"\t")
#          outfile.write(str(GutHitCounts[key])+"\n")
#        elif " " in key:
#          for i in key.split():
#            outfile.write(str(meta_lookup(i))+"\t")
#          outfile.write(str(GutHitCounts[key])+"\n")
#        else:
#          outfile.write(str(key)+"\t")
#          outfile.write(str(GutHitCounts[key])+"\n")
                                                                                                                                      

# Printing output of counting to screen
                            
#      print "\nGut_Hits_All\n"
  
#      for key in GutHitCounts.keys():  
#        if meta_lookup(key):
#          print str(meta_lookup(key))+"\t"+str(GutHitCounts[key])
#        elif " " in key:
#          for i in key.split():
#            print (meta_lookup(i)),
#          print "\t"+str(GutHitCounts[key])
#        else:
#          print str(key)+"\t"+str(GutHitCounts[key])
      

#      print "\nGut_Hits_Fair\n"
#      for key in GutHitCounts_fair.keys():
#        if meta_lookup(key):
#          print str(meta_lookup(key))+"\t"+str(GutHitCounts_fair[key])
#        else:
#          print str(key)+"\t"+str(GutHitCounts_fair[key])
    
#      print "\nGut_Hits_unique_only\n"

#      for key in GutHitCounts_unique.keys():
#        if meta_lookup(key):
#          print str(meta_lookup(key))+"\t"+str(GutHitCounts_unique[key])
#        else:
#          print str(key)+"\t"+str(GutHitCounts_unique[key])
                     
#      print "\n"
                     

## convert tables to metadata

