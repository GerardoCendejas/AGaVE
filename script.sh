#!/bin/bash

# Análisis de calidad de las lecturas usando FastQC

fastqc F73_1.fastq.gz

fastqc F73_2.fastq.gz

# Resumir en un sólo archivo con MultiQC

multiqc .

#Descomprimir archivos fastq.gz

gunzip F73_1.fastq.gz F73_2.fastq.gz

# Separar los archivos en partes más pequeñas con fastq-splitter

n_parts=10

fastq-splitter.pl --n-parts $n_parts --check F73_1.fastq

fastq-splitter.pl --n-parts $n_parts --check F73_2.fastq

# Mover a carpeta nueva

mkdir parts

mv F73_*.part* ./parts/