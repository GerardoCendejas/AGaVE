import sys

from Bio import SeqIO

import pandas as pd

import pysam

import re

import pandas as pd

# Load gbk file

aln_file=sys.argv[3]

viruses_file=sys.argv[4]

records = list(SeqIO.parse(sys.argv[1], "genbank"))

# outdir

outdir = sys.argv[5]

if outdir[-1]!="/":
    outdir = outdir+"/"

# Escribir archivo fasta para todos los ORF detectados

sourceFile = open(f'{outdir}{sys.argv[2]}_all.fasta', 'w')

for record in records:
    for feature in record.features:
        if feature.type == "CDS":
            print(f'> {feature.qualifiers["locus_tag"][0]}\n{feature.extract(record.seq)}',file=sourceFile)

sourceFile.close()

# Escribir archivo fasta con secuencias de AMINOACIDOS de todos los ORF detectados

sourceFile = open(f'{outdir}{sys.argv[2]}_all_aa.fasta', 'w')

for record in records:
    for feature in record.features:
        if feature.type == "CDS":
            print(f'> {feature.qualifiers["locus_tag"][0]}\n{feature.qualifiers["translation"][0]}',file=sourceFile)

sourceFile.close()

# Escribir archivo fasta con secuencias con anotacion diferente a hypothetical protein

sourceFile = open(f'{outdir}{sys.argv[2]}_annotated.fasta', 'w')

for record in records:
    for feature in record.features:
        if feature.type == "CDS" and "product" in feature.qualifiers:
            if feature.qualifiers["product"][0]!="hypothetical protein":
                print(f'> {feature.qualifiers["locus_tag"][0]} {feature.qualifiers["product"][0]}\n{feature.extract(record.seq)}',file=sourceFile)

sourceFile.close()

# Escribir archivo fasta con secuencias de AMINOACIDOS con anotacion diferente a hypothetical protein

sourceFile = open(f'{outdir}{sys.argv[2]}_annotated_aa.fasta', 'w')

for record in records:
    for feature in record.features:
        if feature.type == "CDS" and "product" in feature.qualifiers:
            if feature.qualifiers["product"][0]!="hypothetical protein":
                print(f'> {feature.qualifiers["locus_tag"][0]} {feature.qualifiers["product"][0]}\n{feature.qualifiers["translation"][0]}',file=sourceFile)

sourceFile.close()

# Escribir archivo tabular con las anotaciones de cara transcrito y el ID de la proteína a la cuál es idéntico.

p = re.compile(r'L\wC\wE')

sourceFile = open(f'{outdir}{sys.argv[2]}_annotated.csv', 'w')

print("Contig,LocusTag,ContigLen,LocusLen,ProteinLen,Product,Similar2,LXCXE",file=sourceFile)

for record in records:
    for feature in record.features:
        if feature.type == "CDS" and "product" in feature.qualifiers:
            if feature.qualifiers["product"][0]!="hypothetical protein":
                print(f'{record.name},{feature.qualifiers["locus_tag"][0]},{abs(record.features[0].location.end)},{abs(feature.location.end-feature.location.start)},{len(feature.qualifiers["translation"][0])},"{feature.qualifiers["product"][0]}",{feature.qualifiers["inference"][1].split(":")[2]},{len(p.findall(feature.qualifiers["translation"][0]))!=0}',file=sourceFile)

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

# Importando archivo en formato SAM

samfile = pysam.AlignmentFile(aln_file)

### Creando archivo que dice el número de  contigs alineados a cada reference genome

viruses = pd.read_csv(viruses_file)

viruses = viruses.set_index(viruses["ID"]) 

virus_name = []
genomes = []
count=[]
genome_len=[]
coverages_len = []
coverages_per = []

for genome in samfile.get_index_statistics():
    if genome[1]>0:
        if genome[0].split(".")[0] in viruses.index:
            virus_name.append(viruses.loc[genome[0].split(".")[0]]["Name"])
        else:
            virus_name.append(genome[0].split("|")[0])
        genomes.append(genome[0])
        count.append(genome[1])
        genome_len.append(samfile.get_reference_length(genome[0]))
        coverage = 0
        iter = samfile.fetch(genome[0],0,samfile.get_reference_length(genome[0]))
        for x in iter:
            coverage+=get_mapped_len(str(x).split('\t')[5])
        coverages_len.append(coverage)
        coverages_per.append(coverage/samfile.get_reference_length(genome[0]))
        
data = {"Virus":virus_name,"Genome":genomes,"MappedContigs":count,
        "GenomeLength":genome_len,"GenomeMappedLen":coverages_len,"GenomeCoverageIdx":coverages_per}
data = pd.DataFrame(data,index=None)
data=data.sort_values(by=['GenomeCoverageIdx'],ascending=False)

data.to_csv(f"{outdir}{sys.argv[2]}_ref_genome_maps.csv",index=False)

### Creando archivo de información de alineamiento.

sourceFile = open(f'{outdir}{sys.argv[2]}_aln.csv', 'w')

for i in range(0,len(samfile.get_index_statistics()),1):
    if samfile.get_reference_name(i).split(".")[0] in viruses.index:
            name = viruses.loc[samfile.get_reference_name(i).split(".")[0]]["Name"]
    else:
            name = samfile.get_reference_name(i).split("|")[0]
    ref = samfile.get_reference_name(i)
    iter = samfile.fetch(ref,0,samfile.get_reference_length(ref))
    for x in iter:
        print('%s,%s,%s,%s,%s,"%s"'%(name,str(x).split('\t')[0],aln[int(str(x).split('\t')[1])],samfile.get_reference_length(ref),str(x).split('\t')[3],get_cigar(str(x).split('\t')[5])),file=sourceFile)

sourceFile.close()

# Adding information about annotations

file = pd.read_csv(f'{outdir}{sys.argv[2]}_aln.csv',header=None)

def get_annot_from_contig(record):
    res = []
    for feature in record.features:
        if feature.type=="CDS" and "product" in feature.qualifiers:
            try:
                gene = feature.qualifiers["gene"][0]
            except:
                gene = feature.qualifiers["product"][0].split(",")[0]
            res.extend([gene,str(feature.location.start+1),str(feature.location.end)])
    return(",".join(res))

annot = {}

for record in records:
    annot[record.name]=get_annot_from_contig(record)

annot_col = []

for contig in file[1]:
    annot_col.append(annot[contig.split("|")[2]])

file[6] = annot_col

file.to_csv(f'{outdir}{sys.argv[2]}_aln.csv',index=False,header=None)