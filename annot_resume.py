import sys

from Bio import SeqIO

import pandas as pd

import pysam

import re

#import pandas as pd

# Load gbk file

aln_file="sorted.bam"

records = list(SeqIO.parse(sys.argv[1], "genbank"))

# Escribir archivo fasta con secuencias con anotacion diferente a hypothetical protein

sourceFile = open(f'./{sys.argv[2]}_annotated.fasta', 'w')

for record in records:
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers["product"][0]!="hypothetical protein":
            print(f'> {feature.qualifiers["locus_tag"][0]}\n{feature.extract(record.seq)}',file=sourceFile)

sourceFile.close()

# Escribir archivo fasta con secuencias de AMINOACIDOS con anotacion diferente a hypothetical protein

sourceFile = open(f'./{sys.argv[2]}_annotated_aa.fasta', 'w')

for record in records:
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers["product"][0]!="hypothetical protein":
            print(f'> {feature.qualifiers["locus_tag"][0]}\n{feature.qualifiers["translation"][0]}',file=sourceFile)

sourceFile.close()

# Escribir archivo tabular con las anotaciones de cara transcrito y el ID de la proteína a la cuál es idéntico. 

sourceFile = open(f'./{sys.argv[2]}_annotated.csv', 'w')

print("Contig,LocusTag,ContigLen,LocusLen,ProteinLen,Product,Similar2",file=sourceFile)

for record in records:
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers["product"][0]!="hypothetical protein":
            print(f'{record.name},{feature.qualifiers["locus_tag"][0]},{abs(record.features[0].location.end)},{abs(feature.location.end-feature.location.start)},{len(feature.qualifiers["translation"][0])},"{feature.qualifiers["product"][0]}",{feature.qualifiers["inference"][1].split(":")[2]}',file=sourceFile)

sourceFile.close()

### Creando diccionarios necesarios

strand = {0:1,16:-1,256:1,272:-1,2048:1,2064:-1}

aln = {0:"primary",16:"primary",256:"not primary",272:"not primary",2048:"supplementary",2064:"supplementary"}

mapped = {"M":1,"D":0,"N":0,"I":2,"S":0,"H":0}
keys = list(mapped.keys())

### Definiendo función para interpretar CIGAR

def get_cigar(x):
    num = re.split(r'\D+',x)
    tmp = re.split(r'\d',x)
    num.remove("")
    aligned = []
    for i in tmp:
        if i in keys:
            aligned.append(i)
            
    final = [mapped[aligned[0]]]
    
    final.append(int(num[0]))
    
    actual = mapped[aligned[0]]
    
    for i in range(1,len(num)):
        
        if mapped[aligned[i]] == 2:
            next        
        if mapped[aligned[i]] == actual:
            final[-1]+=int(num[i])
        elif mapped[aligned[i]] != actual and int(num[i]) <= 10:
            final[-1]+=int(num[i])
        elif mapped[aligned[i]] != actual and int(num[i]) > 10:
            actual = mapped[aligned[i]]
            final.append(int(num[i]))
            
    final = list(map(str,final))
            
    return(",".join(final))

# Importando archivo en formato SAM

samfile = pysam.AlignmentFile(aln_file)


### Creando archivo que dice el núumero de  contigs alineados a cada reference genome

genomes = []
count=[]
genome_len=[]

for genome in samfile.get_index_statistics():
    if genome[1]>0:
        genomes.append(genome[0])
        count.append(genome[1])
        genome_len.append(samfile.get_reference_length(genome[0]))
        
data = {"Genome":genomes,"MappedContigs":count,"GenomeLength":genome_len}
data = pd.DataFrame(data,index=None)
data=data.sort_values(by=['MappedContigs'],ascending=False)

data.to_csv(f"{sys.argv[2]}_ref_genome_maps.csv",index=False)

### Creando archivo de información de alineamiento.

sourceFile = open(f'./{sys.argv[2]}_aln.csv', 'w')

for i in range(0,len(samfile.get_index_statistics()),1):
    ref = samfile.get_reference_name(i)
    iter = samfile.fetch(ref,0,samfile.get_reference_length(ref))
    for x in iter:
        print('%s,%s,%s,%s,%s,%s,"%s"'%(ref,str(x).split('\t')[0],aln[int(str(x).split('\t')[1])],samfile.get_reference_length(ref),strand[int(str(x).split("\t")[1])],str(x).split('\t')[3],get_cigar(str(x).split('\t')[5])),file=sourceFile)

sourceFile.close()