# Importing requiered libraries

import sys
from Bio import SeqIO
import os
import pandas as pd
import pysam
import re

# outdir

outdir = sys.argv[7]

if outdir[-1]!="/":
    outdir = outdir+"/"



# Loading bam file of mapped reads

aln_file = sys.argv[1]

contig_gbk = sys.argv[2]

contig_cds = sys.argv[3]

vir_found = sys.argv[4]

vir_fasta = sys.argv[5]





# Creating fasta files for viral genomes

# Loading fasta file with contig sequences

records = list(SeqIO.parse(vir_fasta, "fasta"))

# Writing only viral contigs into file

sourceFile = open(f'{outdir}{sys.argv[6]}_vir.fasta', 'w')

for record in records:

    if record.id in file[,0]:

        print(f'> {record.id}\n{record.seq}',file=sourceFile)

sourceFile.close()





# Loading file with mapped genomes

file = pd.read_csv(vir_found)

# Importing SAM file

samfile = pysam.AlignmentFile(aln_file)

# Results files, adding count of reads mapped to viral genome...

genomes = []
read_count=[]
genome_len=[]
read_mapped_len = []
copy_number = []

for genome in samfile.get_index_statistics():
    if genome[0] in vir_found[0]:
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


data = pd.merge(file, data, on='Genome')

data = data.sort_values(by=['GenomeCoverageIdx'],ascending=False)

data.to_csv(vir_found,index=False)





# Genebank format to GFF3 of viral contigs (for lifting over)
# And generating fasta file of viral contigs

# Loading file with info of viral contigs

file = pd.read_csv(contig_cds)

# Loading annotation file for contigs

records = list(SeqIO.parse(contig_gbk, "genbank"))

sourceFile = open(f'{outdir}{sys.argv[6]}_vir_contig.gff', 'w')

sourceFile_fasta = open(f'{outdir}{sys.argv[6]}_vir_contig.fasta', 'w')

for record in records:
    if record.id in file[,1]:

        print(f'> {record.id}\n{record.seq}',file=sourceFile_fasta)



        for feature in record.features:

            if feature.type == "CDS":

                try:
                    gene = feature.qualifiers["gene"][0]
                except:
                    gene = feature.qualifiers["product"][0].split(",")[0]

                print(f'{record.id}\tprokka\tCDS\t{feature.location.start}\t{feature.location.end}\t.\t{feature.strand}\t.\tID={feature.qualifiers["locus_tag"]};Name={gene}',file=sourceFile)

sourceFile.close()

sourceFile_fasta.close()