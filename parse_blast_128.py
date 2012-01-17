#!/usr/bin/python                                                                                                                                                                             
# Author: Sathish <subramanians@wusm.wustl.edu                                                 \
                                                                                              
import sys,os
import warnings

print "\nfirst argument is blast output, second is the filename/directory for output"

blastfile = open(sys.argv[1],"r")

os.makedirs(sys.argv[2]+"_processed")

outfile = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_bychild.txt","w")
outfile2 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_gut_hits","w")
outfile3 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_gut_hits_length_corrected","w")
outfile4 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_genus_species_hits","w")
outfile5 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_genus_species_hits_length_corrected","w")

# unique to 128

metadatafile = open("/home/comp/jglab/sathish/Satz_Scripts/support_files/128_Genomes_Metadata.txt","r")
genome_dict={}
genome_name={}
genus_sp={}
genus_sp_len={}
total_len=0


Gut_dict={}
Gut_len={}

for line in metadatafile:
    if "#" not in line: 
        (Genome_ID,Genome_Name,Genus_Species,Genome_Length)=line.split("\t")
        genome_dict[Genome_ID]=Genus_Species 
        genome_name[Genome_ID]=Genome_Name
        genus_sp[Genus_Species]=0
        genus_sp_len[Genus_Species]=0
        total_len=total_len+float(Genome_Length.strip())
        Gut_dict[Genome_ID]=0
        Gut_len[Genome_ID]=float(Genome_Length.strip())

###########################

print float(total_len/(len(genome_dict.keys())))

ID_dict={}
ID_nonrep={}

lastID=str("last_queryid")
total_fasta=5570500
total_assigned=0

for line in blastfile:
    (queryID, qlen, sseqid, slen, evalue, bitscore, score, length, pident, mismatch, gapopen) = line.split("\t")
    #childID added to dictionary ID_dict
    childID = queryID.split("_")[1]
    childID = childID.strip()
    if childID not in ID_dict.keys(): #initiate ID_dict key if new childID
        ID_dict[childID]=0
        ID_nonrep[childID]=0
    if queryID!=lastID: #if the queryID is not the top hit ignore and do not add to it
        ID_dict[childID]=ID_dict[childID]+1
        if sseqid not in Gut_dict.keys():
            print "sseqid not in Gut_dict "+sseqid
        Gut_dict[sseqid]=float(Gut_dict[sseqid])+1
    ID_nonrep[childID]=ID_nonrep[childID]+1
    lastID=queryID

#print str(len(Gut_dict.keys()))+" should be 128\n"

total_nonrep=0
total_derep=0

for key in sorted(ID_dict.keys()):
    outfile.write(key+"\t"+str(ID_dict[key])+"\t"+str(ID_nonrep[key])+"\n")
    total_derep=total_derep+ID_dict[key]
    total_nonrep=total_nonrep+ID_nonrep[key]

print "\nAverage length of genomes "+str(total_len/len(Gut_dict))

corr_norm_total=0
 
for key in sorted(Gut_dict.keys()):
    unnorm_count=(Gut_dict[key]*100)/total_derep
    uncorr_norm_total=((Gut_dict[key]*100)/total_fasta)
    corr_norm_total=corr_norm_total+(((Gut_dict[key]*100)/total_fasta)*((total_len/len(Gut_dict))/Gut_len[key])) 
    total_assigned=total_assigned+uncorr_norm_total
    outfile2.write(genome_name[key]+"\t"+str(uncorr_norm_total)+"\n")
    genus_sp[(genome_dict[key])]=genus_sp[(genome_dict[key])]+uncorr_norm_total

corr_factor=(total_derep*100/total_fasta)/corr_norm_total
print "Correction factor to make assigned add up to 100 with unassigned) "+str(corr_factor)
corr_totaled = 0

for key in sorted(Gut_dict.keys()):
    corr_total=(((Gut_dict[key]*100)/total_fasta)*((total_len/len(Gut_dict))/Gut_len[key]))*corr_factor
    corr_totaled=corr_totaled+corr_total
    outfile3.write(genome_name[key]+"\t"+str(corr_total)+"\n")
    genus_sp_len[(genome_dict[key])]=genus_sp_len[(genome_dict[key])]+corr_total

for key in sorted(genus_sp.keys()):
    outfile4.write(key+"\t"+str(genus_sp[key])+"\n") 
    outfile5.write(key+"\t"+str(genus_sp_len[key])+"\n")


outfile2.write("Unassigned\t"+str(100-total_assigned)+"\n")
outfile4.write("Unassigned\t"+str(100-total_assigned)+"\n")
outfile3.write("Unassigned\t"+str(100-corr_totaled)+"\n")
outfile5.write("Unassigned\t"+str(100-corr_totaled)+"\n")

print("\n"+"Fraction of reads assigned: "+"\t"+str(float(total_derep)/total_fasta))
print("\n"+"Total number of metagenomes:"+"\t"+str(len(ID_dict)))
print("\n"+"Total number of reads assigned at least once"+"\t"+str(total_derep))
print("\n"+"Average number of assignments considering first hit only:"+"\t"+str(total_derep/len(ID_dict)))                               
print("\n")
