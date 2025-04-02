# Importing requiered libraries

import sys
from Bio import SeqIO
import os
import pandas as pd
import pysam
import re

# output

outdir = sys.argv[4]

if outdir[-1]!="/":
    outdir = outdir+"/"



# Loading bam file of mapped reads

aln_file = sys.argv[1]

vir_found = sys.argv[2]

sample = sys.argv[3]



# Loading file with mapped genomes

file = pd.read_csv(vir_found)

found_genomes = file["Genome"]

# Importing SAM file

samfile = pysam.AlignmentFile(aln_file)

# Dictionaries

strand = {0:1,16:-1,256:1,272:-1,2048:1,2064:-1}

aln = {0:"primary",16:"primary",256:"not primary",272:"not primary",2048:"supplementary",2064:"supplementary"}

mapped = {"M":1,"D":0,"N":0,"I":2,"S":0,"H":0}

keys = list(mapped.keys())

# Needed function 

def get_mapped_len(x):
    num = re.split(r'\D+',x)
    tmp = re.split(r'\d',x)
    num.remove("")
    
    aligned = []
    
    for i in tmp:
        if i in keys:
            aligned.append(i)
            
    mapped_len = 0
    
    for i in range(0,len(num)):
        
        if mapped[aligned[i]] == 1:
            mapped_len+=int(num[i])       
        else:
            next
            
    return(mapped_len)

# Results files, adding count of reads mapped to viral genome...

if (os.path.getsize(f"viralmap/{sample}_mapped2virus.fasta") != 0):

    print("It's checking for counts")

    genomes = []
    read_count=[]
    genome_len=[]
    read_mapped_len = []
    copy_number = []

    for genome in samfile.get_index_statistics():
        if genome[0] in found_genomes.values:

            print(f"cheking counts for {genome[0]}")

            genomes.append(genome[0])
            read_count.append(genome[1])
            genome_len.append(samfile.get_reference_length(genome[0]))
            coverage = 0
            iter = samfile.fetch(genome[0],0,samfile.get_reference_length(genome[0]))
            for x in iter:
                coverage+=get_mapped_len(str(x).split('\t')[5])
            read_mapped_len.append(coverage)
            copy_number.append(coverage/samfile.get_reference_length(genome[0]))
                
    data = {"Genome":genomes,"ReadCount":read_count,
                "ReadMappedLen":read_mapped_len,"CopyNumber":copy_number}
    data = pd.DataFrame(data,index=None)


    try:
        data = pd.merge(file, data, on='Genome')

        data = data.sort_values(by=['GenomeCoverageIdx'],ascending=False)

        data.to_csv(vir_found,index=False)
    except:
        print("")



# Merging the gff files of the found genomes

gff_file = pd.DataFrame()

if (os.path.getsize(f"viralmap/{sample}_mapped2virus.fasta") != 0):

    for genome in found_genomes:

        gff_temp = pd.read_csv(f"viral_genomes/{genome}.gff",sep="\t")

        gff_file = pd.concat([gff_file, gff_temp])

gff_file.to_csv(f"{outdir}{sample}.gff",index=False,sep="\t")






