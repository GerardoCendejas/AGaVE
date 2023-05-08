#!/bin/bash

source /labs/csbig/gerardo/miniconda3/etc/profile.d/conda.sh

# Cargamos las funciones que se utilizarán

MYSELF="$(realpath "$0")"

MYDIR="${MYSELF%/*}"

source "$MYDIR/functions.sh"

####### Nombramiento de variables #######

sample="SRR9298859"  # Muestra a analizar

min_qual=28 # Calidad mínima de las lecturas para la limpieza

min_len=75 # Longitud mínima de las lecturas para la limpieza

genome="/labs/csbig/gerardo/genomes/human/STAR" # Path to reference genome

virusdb="/labs/csbig/gerardo/genomes/virus/virus.mmi" # Path to virus database minimap2 index

### Limpieza de las lecturas ###

CleanReads $sample $min_qual $min_len $MYDIR

### Mapeo a host

HostMapping $sample $genome

### Ensamblaje de novo de lecturas no mapeadas

Assembly $sample

### Mapear los contigs a la base de datos de genomas virales

VirusMapping $sample $virusdb

### Anotando los contigs encontrados

AnnotatingContigs $sample


