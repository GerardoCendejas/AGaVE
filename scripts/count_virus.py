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

print(file)

found_genomes = file["Genome"]

# Importing SAM file

samfile = pysam.AlignmentFile(aln_file)

# Results files, adding count of reads mapped to viral genome...

if (os.path.getsize(f"viralmap/{sample}_mapped2virus.fasta") != 0):

    genomes = []
    read_count=[]
    genome_len=[]
    read_mapped_len = []
    copy_number = []

    for genome in samfile.get_index_statistics():
        if genome[0] in found_genomes:
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

    print(data)

    data = pd.concat(file, data, on='Genome')

    data = data.sort_values(by=['GenomeCoverageIdx'],ascending=False)

    data.to_csv(vir_found,index=False)





# Merging the gff files of the found genomes

gff_file = []

if (os.path.getsize(f"viralmap/{sample}_mapped2virus.fasta") == 0):

    for genome in found_genomes:

        gff_temp = pd.read_csv(f"viral_genomes/{genome}.gff",sep="\t")

        gff_file = pd.concat([gff_file, gff_temp])

data.to_csv(f"{outdir}{sample}.gff",index=False,sep="\t")






