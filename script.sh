#!/bin/bash

# Cargamos las funciones que se utilizarán

MYSELF="$(realpath "$0")"

MYDIR="${MYSELF%/*}"

source "$MYDIR/functions.sh"

####### Nombramiento de variables #######

sample="SRR9298859"  # Muestra a analizar

n_parts=11 # Número de partes en dividir la muestra

min_qual=28 # Calidad mínima de las lecturas para la limpieza

min_len=75 # Longitud mínima de las lecturas para la limpieza

genome="/labs/csbig/gerardo/genomes/human/STAR" # Path to reference genome

### Limpieza de las lecturas ###

CleanReads $sample $min_qual $min_len $MYDIR

### Mapeo a host

HostMapping $sample $genome


