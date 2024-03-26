import sys

from Bio import SeqIO

import pandas as pd

import pysam

import re

import pandas as pd

# Load gbk file

aln_file=sys.argv[3]

records = list(SeqIO.parse(sys.argv[1], "genbank"))

# outdir

outdir = sys.argv[4]

if outdir[-1]!="/":
    outdir = outdir+"/"

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

### Creando archivo de información de alineamiento.

sourceFile = open(f'{outdir}{sys.argv[2]}_int_aln.csv', 'w')

for i in range(0,len(samfile.get_index_statistics()),1):
    ref = samfile.get_reference_name(i)
    iter = samfile.fetch(ref,0,samfile.get_reference_length(ref))
    for x in iter:
        print('%s,%s,%s,%s,%s,"%s","%s"'%(ref,str(x).split('\t')[0],aln[int(str(x).split('\t')[1])],samfile.get_reference_length(ref),str(x).split('\t')[3],get_cigar(str(x).split('\t')[5]),str(x).split('\t')[5]),file=sourceFile)

sourceFile.close()

# Adding information about annotations

file = pd.read_csv(f'{outdir}{sys.argv[2]}_int_aln.csv',header=None)

contig_len = {}

for record in records:
    contig_len[record.name]=len(record.seq)

contig_len_fin = []

for contig in file[1]:
    contig_len_fin.append(contig_len[contig.split("|")[2]])

file[7] = contig_len_fin

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

file[8] = annot_col

file_2 = pd.read_csv(f'{outdir}{sys.argv[2]}_aln.csv',header=None)

viral = []

for contig in file[1]:
    if contig in file_2[1].values:
        viral.append(1)
    else:
        viral.append(0)

file[9] = viral

file.to_csv(f'{outdir}{sys.argv[2]}_int_aln.csv',index=False,header=None)


# To get proper chromosome number for GRCh38.p14

def id_chr(id):
    chr = int(id.split("_")[1].split(".")[0])

    if chr<23:
        return(str(chr))
    elif chr==23:
        return("X")
    elif chr==24:
        return("Y")
    else:
        return(None)


sourceFile = open(f'{outdir}{sys.argv[2]}_int_id.csv', 'w')

print("Type,Shape,Chr,Start,End,color,len",file=sourceFile)

for i in range(0,len(file[0]),1):
    if  (id_chr(file.iat[i,0]) is not None) & (file.iat[i,2]=="primary"):
        if (file.iat[i,9]==1) and (file.iat[i,8]!=""):
            print('annot_viral,triangle,"%s",%s,%s,0000ff,%s'%(id_chr(file.iat[i,0]),file.iat[i,4],file.iat[i,4]+100,file.iat[i,6]),file=sourceFile)
        elif (file.iat[i,9]==1) and (file.iat[i,8]==""):
            print('non_annot_viral,triangle,"%s",%s,%s,ffa500,%s'%(id_chr(file.iat[i,0]),file.iat[i,4],file.iat[i,4]+100,file.iat[i,6]),file=sourceFile)
        elif (file.iat[i,9]==0) and file.iat[i,8]!="":
            print('annot_contig,circle,"%s",%s,%s,000000,%s'%(id_chr(file.iat[i,0]),file.iat[i,4],file.iat[i,4]+100,file.iat[i,6]),file=sourceFile)
        else:
            print('non_annot_contig,circle,"%s",%s,%s,ffc0cb,%s'%(id_chr(file.iat[i,0]),file.iat[i,4],file.iat[i,4]+100,file.iat[i,6]),file=sourceFile)

sourceFile.close()
