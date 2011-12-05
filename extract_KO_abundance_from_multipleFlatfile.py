#!/usr/bin/python                                                             
# Author: Sathish <subramanians@wusm.wustl.edu                                                                                                               
"""                                                                           
This is a python script that takes as input a tab-delimited file and prints out all lines with the KO of interest                                                                                                 
"""

import sys,os

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
    elif key =="????":
      return "Unassigned"
    elif key == "More_than_one_best_hit":
      return key

#Loaded all species genome_ids into GutDB_metadata_list

a = open('/home/comp/jglab/sathish/Bangladesh/128_Genomes_Metadata.txt',"r")
meta = a.readlines()
GutDB_metadata_list=[]
species_number=0
for line in meta:
  (Genome_ID, Genome_Name, Genus_Species, Genome_Length) = line.split("\t")
  if Genome_ID not in GutDB_metadata_list:
    if Genome_ID == "#Genome_ID":
      print "\nReading GenomeIDs from "+"/home/comp/jglab/sathish/Bangladesh/128_Genomes_Metadata.txt\n"
    else:
      GutDB_metadata_list.append(Genome_ID)
      species_number=species_number+1    
print "Loaded "+str(species_number)+" species\n"
print GutDB_metadata_list
GutDB_metadata_list.append('????')
GutDB_metadata_list.append('More_than_one_best_hit')

path = sys.argv[1]
K0=sys.argv[2]
  
SampleList=[]
Sample_Path_Dict={}
count_file = 0

for root, dirs, files in os.walk(path):
  for name in files:
    if 'FastaFlatFile.txt.gz' in name:
      infile = os.path.join(root, name)
      sample_id = (name.split('_FastaFlatFile.txt.gz')[0])
      SampleList.append(sample_id)
      count_file = count_file + 1
      #infile = path+sample_id+"/"+name
      #infile = path+name  - need to change to not matter the naming scheme
      Sample_Path_Dict[sample_id]=infile

print "\nFound "+str(count_file)+" FlatFiles in the path\n"
print sorted(SampleList)
print "\nExtracting all GutDB assignment information for reads that map to "+K0+" from Flatfiles\n"

# initializing dictionaries and lists for each FlatFile                                                                                                      
GutHitList=[] # List of all hits                                               
GutHitCounts={}  # dictionary with key being genome hit and value being count                                                                                 
GutHitCounts_unique={}
GutHitCounts_fair={}  # values are for each flatfile
M_GutHitCounts_fair={} # values are lists in series for all flatfiles
M_GutHitCounts_unique={}

for i in GutDB_metadata_list:
  M_GutHitCounts_fair[i]=[]
  M_GutHitCounts_unique[i]=[]


for file in sorted(Sample_Path_Dict):
  count_reads=0
  path_to_flatfile=Sample_Path_Dict[file]
  for line in gziplines(path_to_flatfile):
    (Header, SequenceLength, TooShort, TooManyN, Replicate, HostGenomeHit, PassFilter, KeggBlastHit, KeggAnnotatedHit, MeropsHit, CogBlastHit, CogAnnotatedHit, GutGenomeHit) = line.split("\t")
    if PassFilter=="1":
      count_reads = count_reads+1
      if K0 in line:
        gutGenomeHit=GutGenomeHit.strip("\n")
        GutHitList.append(gutGenomeHit)

# Count Gut Hits resulting in two type of dictionaries to follow on fair and unique which has an extra label 'More than one best hit'
  for i in GutHitList:
    GutHitCounts[i] = GutHitList.count(i)
    if " " not in i:
      GutHitCounts_unique[i]=GutHitCounts[i]
      GutHitCounts_fair[i]=GutHitCounts[i]
    
  for key in GutHitCounts.keys():
    if " " in key:
      if 'More_than_one_best_hit' in GutHitCounts_unique.keys():
        GutHitCounts_unique['More_than_one_best_hit']=GutHitCounts_unique['More_than_one_best_hit']+GutHitCounts[key]
      else:
        GutHitCounts_unique['More_than_one_best_hit']=GutHitCounts[key]
      Derep=key.split()
      for i in set(Derep):
        if i in GutHitCounts_fair:
          GutHitCounts_fair[i]=(GutHitCounts_fair[i])+(0.5*GutHitCounts[key])
        else:
          GutHitCounts_fair[i]=(0.5*GutHitCounts[key])
  
  for i in range(len(GutDB_metadata_list)):
    s=GutDB_metadata_list[i]
    if s in GutHitCounts_fair.keys():
      M_GutHitCounts_fair[s].append(((GutHitCounts_fair[s])/count_reads)*100) 
    else:
      M_GutHitCounts_fair[s].append(0)  
    if s in GutHitCounts_unique.keys():
      M_GutHitCounts_unique[s].append(((GutHitCounts_unique[s])/count_reads)*100)
    else:
      M_GutHitCounts_unique[s].append(0)


  GutHitList[:]=[]
  GutHitCounts.clear()
  GutHitCounts_unique.clear()
  GutHitCounts_fair.clear()

Final_GutHitsCounts_fair={}
for s in M_GutHitCounts_fair:
  if sum(M_GutHitCounts_fair[s])>0:
    Final_GutHitsCounts_fair[s]=M_GutHitCounts_fair[s]

Final_GutHitsCounts_unique={}
for s in M_GutHitCounts_unique:
  if sum(M_GutHitCounts_unique[s])>0:
    Final_GutHitsCounts_unique[s]=M_GutHitCounts_unique[s]

# From dictionary to tab-delimited file output

path_out=path.strip('/.')
path_out=path_out.replace('/','_')
outfile = open("GutDB_Hits_"+K0+"_"+path_out+"_fair.txt","w")
outfile.write("Sample_ID\t")
for i in sorted(SampleList):
  outfile.write(str(i)+"\t")

for key in Final_GutHitsCounts_fair.keys():
  outfile.write("\n")
  outfile.write(str(meta_lookup(key)))
  for i in range(len(Final_GutHitsCounts_fair[key])):
    outfile.write("\t"+str((Final_GutHitsCounts_fair[key])[i]))
outfile.write("\n")

outfile_u = open("GutDB_Hits_"+K0+"_"+path_out+"_unique.txt","w")
outfile_u.write("Sample_ID\t")
for i in sorted(SampleList):
  outfile_u.write(str(i)+"\t")

for key in Final_GutHitsCounts_unique.keys():
  outfile_u.write("\n")
  outfile_u.write(str(meta_lookup(key)))
  for i in range(len(Final_GutHitsCounts_unique[key])):
    outfile_u.write("\t"+str((Final_GutHitsCounts_unique[key])[i]))
outfile_u.write("\n")
