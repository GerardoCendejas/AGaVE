#!/bin/bash

# Cargamos las funciones que se utilizarán

MYSELF="$(realpath "$0")"

MYDIR="${MYSELF%/*}"

source "$MYDIR/functions.sh"

####### Nombramiento de variables #######

sample="prueba"  # Muestra a analizar

n_parts=11 # Número de partes en dividir la muestra

min_qual=28 # Calidad mínima de las lecturas para la limpieza

min_len=75 # Longitud mínima de las lecturas para la limpieza

####### Análisis de calidad de las lecturas usando FastQC y MultiQC #######

QualityControl $sample

### Separar los archivos en partes más pequeñas y mover a propia carpeta ###

SplitFiles $sample $n_parts

### Limpieza de las lecturas ###

CleanReads $sample $n_parts $min_qual $min_len