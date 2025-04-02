# Importing requiered libraries

import sys
from Bio import SeqIO
import os
import pandas as pd
import pysam
import re

# outdir

outdir = sys.argv[4]

if outdir[-1]!="/":
    outdir = outdir+"/"


# Input files

vir_found = sys.argv[1]

vir_fasta = sys.argv[2]

sample = sys.argv[3]

cores = sys.argv[5]

# output_file

output = sys.argv[6]



# Loading file with mapped genomes

file = pd.read_csv(vir_found)

genomes = file["Genome"]

# Get files present

files = os.listdir(outdir)

# fasta file of all genomes

records = list(SeqIO.parse(vir_fasta, "fasta"))

# Start

for genome in genomes:

    # Check if genome was previously checked

    if f"{genome}.gff" not in files:

        print(f"Checking for genome {genome}")

        for record in records:

            if record.id == genome:
                sourceFile_fasta = open(f'{outdir}genome.fas', 'w')

                print(f'> {record.id}\n{record.seq}',file=sourceFile_fasta)

                sourceFile_fasta.close()

                print(f"Wrote fasta file for {genome}")

                os.system(f"prokka --addgenes --outdir {outdir}prueba/ --locustag \"{genome}\" --prefix genome --kingdom Viruses --cpus {cores} --norrna --notrna \"{outdir}genome.fas\" --force")

                print(f"Annotated {genome}")

                os.system(f"mv {outdir}prueba/genome.gbk \"{outdir}{genome}.gbk\"")

                os.system(f"mv {outdir}genome.fas \"{outdir}{genome}.fas\"")

                os.system(f"rm -rf {outdir}prueba/")


                annot = list(SeqIO.parse(f"{outdir}{genome}.gbk", "genbank"))


                sourceFile = open(f'{outdir}{genome}.gff', 'w')

                for record in annot:    

                    for feature in record.features:

                        if feature.type == "CDS":

                            try:
                                gene = feature.qualifiers["gene"][0]
                            except:
                                gene = feature.qualifiers["product"][0].split(",")[0]

                            print(f'{record.id}\tprokka\tCDS\t{feature.location.start}\t{feature.location.end}\t.\t{feature.strand}\t.\tID \"{feature.qualifiers["locus_tag"]}\";Name= \"{gene}\";locus \"{feature.qualifiers["locus_tag"]}\";',file=sourceFile)

                sourceFile.close()


                sourceFile = open(output, 'w')
                print(f"Checked for genome: {genome}",file=sourceFile)
                sourceFile.close()


sourceFile = open(output, 'w')
print(f"Checked all genomesfor {sample}",file=sourceFile)
sourceFile.close()





