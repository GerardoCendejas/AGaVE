# (Name)


## Introduction

(Name) is a pipeline designed for the search of viral sequences in RNAseq samples from human tissues or human related samples. This pipeline was designed used previously developed tools that are available via conda and organized using the snakemake workflow management system. A basic diagram of the pipeline is shown in the next image.

<p align="center" width="100%">
    <img width="50%" src="https://github.com/GerardoCMM/RNAseq-Retinoblastoma/blob/main/dag.svg"> 
</p>

## Requirements

### Hardware

- 16GB of memory

### Software

- Linux operating system

- Some version of conda, we recommend the miniconda distribution.

You can use the next commands to install miniconda on a linux computer:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh
```

After running the previous commands, follow the instructions for a full installation.

## Installation

Clone repository:

```
git clone https://github.com/GerardoCMM/RNAseq-Retinoblastoma.git
```

Create conda environment:

```
cd RNAseq-Retinoblastoma

conda env create -f rna_seq_analysis.yml
```

## Usage

### Activating the environment

This is neccesary every time you want to run the pipeline.

```
conda activate rna_seq_analysis
```

### Databases

It is recommended that the databases and their indexes are outside of the repo directory; this should be done only one time.

The human genome can be downloaded from the [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/) and indexed with the next commands:

- STAR
```
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genomic/fasta.fna --sjdbGTFfile /path/to/genomic/annotations.gtf --sjdbOverhang 150
```

Where `/path/to/genomic/fasta.fna` and `/path/to/genomic/annotations.gtf` are the locations for the downloaded files for the genome and `/path/to/genomeDir` is the directory of the index to be created; `--runThreadN` sets the number of threads to be used in the indexing process.

- Minimap2
```
minimap2 -d /path/to/genomic/index.mmi /path/to/genomic/fasta.fna
```

The viral database is present in the repo as `virus.fa` and can be indexed with the following command:
```
minimap2 -d /path/to/viral/index.mmi virus.fa
```

### Config file

The `config.yml` file should be modified to have the proper database directories and number of cores to be used.
```
genome_dir : "/path/to/genomeDir"
virus_db : "/path/to/viral/index.mmi"
human_minimap : "/path/to/genomic/index.mmi"
cores : "8"
```

### Samples

Your samples should be located at the samples/ directory in fastq format and the naming:
- *sample*_1.fastq
- *sample*_2.fastq
 
where _sample_ can be any given name for your sample.

If your samples are in the SRA database you can use [SRA Toolkit](https://github.com/ncbi/sra-tools):
```
mkdir samples && cd samples

fasrterq-dump --split-files -p SRRXXXXXX
```

Where `SRRXXXXXX` is the run ID.

### For running the full pipeline

Checking the steps to be run (from the base of the repo):

```
snakemake -np sample_log.a humanmap/sample_mapped2human.fastq
```

Running the pipeline with a maximum of 8 cores:

```
snakemake --cores 8 --use-conda sample_log.a humanmap/sample_mapped2human.fastq
```
