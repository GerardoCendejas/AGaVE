import sys

from Bio import SeqIO

#import pandas as pd

# Load gbk file

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

for record in records:
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers["product"][0]!="hypothetical protein":
            print(f'{feature.qualifiers["locus_tag"][0]},{feature.qualifiers["product"][0]},{feature.qualifiers["inference"][1].split(":")[2]}',file=sourceFile)

sourceFile.close()



