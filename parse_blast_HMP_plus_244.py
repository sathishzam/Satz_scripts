#!/usr/bin/python                                                                                                                                                                             
# Author: Sathish <subramanians@wusm.wustl.edu                                                 \
                                                                                              
import sys,os
import warnings

blastfile = open(sys.argv[1],"r")
outfile = open(sys.argv[2],"w")

os.makedirs(sys.argv[2]+"_processed")

outfile2 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_gut_hits","w")
outfile3 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_gut_hits_length_corrected","w")
outfile4 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_genus_species_hits","w")
outfile5 = open(sys.argv[2]+"_processed/"+sys.argv[2]+"_genus_species_hits_length_corrected","w")

# unique to 128

metadatafile = open("/home/comp/jglab/sathish/Satz_Scripts/support_files/HMP_metadata.txt","r")
metadatafile2 = open("/home/comp/jglab/sathish/Satz_Scripts/support_files/244_custom_mod_metadata.txt","r")
genome_id={}
genome_dict={}
genome_name={}
genome_length={}
genus_sp={}
genus_sp_len={}

for line in metadatafile:
    if "#" not in line: 
        (nonuniqid,Genome_ID,Genome_Name,Genus_Species,Genome_Length)=line.split("\t")
        genome_id[nonuniqid]=Genome_ID
        genome_dict[Genome_ID]=Genus_Species
        genome_name[Genome_ID]=Genome_Name
        genome_length[Genome_ID]=float(Genome_Length.strip())
        genus_sp[Genus_Species]=0
        genus_sp_len[Genus_Species]=0


for line in metadatafile2:
    if "#" not in line: 
        (Genome_Name,Genus_Species,Genome_Length)=line.split("\t")
        genome_id[(Genome_Name.split("_")[0]).strip()]=Genome_Name
        genome_dict[Genome_Name]=Genus_Species.strip()
        genome_name[Genome_Name]=Genome_Name
        genome_length[Genome_Name]=float(Genome_Length.strip())
        genus_sp[Genus_Species.strip()]=0
        genus_sp_len[Genus_Species.strip()]=0

total_len=0

for key in sorted(genome_length.keys()):
    total_len = total_len + genome_length[key]

#print "total length of genomes "+str(total_len) 

###########################

ID_dict={}
ID_nonrep={}
Gut_dict={}
#Gut_len={}

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
        if "gnl|BL_ORD_ID|" in sseqid:
            sseqnc = (sseqid.split("|")[2])
            if sseqnc not in Gut_dict.keys():
                Gut_dict[sseqnc]=0
                genome_length[(genome_id[sseqnc])]=float(slen)
                total_len = total_len + float(slen)
            Gut_dict[sseqnc]=float(Gut_dict[sseqnc])+1
        else:
            sseqnc = (sseqid.split("|")[3]).split(".")[0]
            if sseqnc not in Gut_dict.keys():
                Gut_dict[sseqnc]=0
#            Gut_len[sseqid]=float(slen)
#            total_len=total_len+float(slen)
            Gut_dict[sseqnc]=float(Gut_dict[sseqnc])+1
    ID_nonrep[childID]=ID_nonrep[childID]+1
    lastID=queryID

print "finished reading blastfile"
print "total length of genomes "+str(total_len) 

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
    corr_norm_total=corr_norm_total+(((Gut_dict[key]*100)/total_fasta)*((total_len/len(genome_length))/genome_length[(genome_id[key])])) 
    total_assigned=total_assigned+uncorr_norm_total
    outfile2.write(genome_name[(genome_id[key])]+"\t"+str(uncorr_norm_total)+"\n")
    genus_sp[(genome_dict[(genome_id[key])])]=genus_sp[(genome_dict[(genome_id[key])])]+uncorr_norm_total

corr_factor=(total_derep*100/total_fasta)/corr_norm_total
#print corr_factor
corr_totaled = 0

for key in sorted(Gut_dict.keys()):
    corr_total=(((Gut_dict[key]*100)/total_fasta)*((total_len/len(genome_length))/genome_length[(genome_id[key])]))*corr_factor
    corr_totaled=corr_totaled+corr_total
    outfile3.write(genome_name[(genome_id[key])]+"\t"+str(corr_total)+"\n")
    genus_sp_len[(genome_dict[(genome_id[key])])]=genus_sp_len[(genome_dict[(genome_id[key])])]+corr_total

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
