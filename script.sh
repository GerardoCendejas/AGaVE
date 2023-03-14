#!/bin/bash

# Cargamos las funciones que se utilizarán

MYSELF="$(realpath "$0")"

MYDIR="${MYSELF%/*}"

source "$MYDIR/functions.sh"

####### Nombramiento de variables #######

sample="prueba"  # Muestra a analizar

n_parts=10 # Número de partes en dividir la muestra

####### Análisis de calidad de las lecturas usando FastQC y MultiQC #######

QualityControl $sample

#Descomprimir archivos fastq.gz  #### En el data set de prueba ya esta descomprimido

#gunzip F73_1.fastq.gz F73_2.fastq.gz

### Separar los archivos en partes más pequeñas y mover a propia carpeta ###

SplitFiles $sample $n_parts