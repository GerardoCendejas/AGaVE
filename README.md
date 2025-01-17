# AGaVe 

## Introduction

AGaVE (Analysis of Genomic and Viral Elements) is a pipeline designed for the search of viral sequences in RNAseq samples from human tissues or human related samples. This pipeline was designed using previously developed tools that are available via conda and organized using the snakemake workflow management system. A basic diagram of the pipeline is shown in the next image.

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
git clone https://github.com/GerardoCendejas/AGaVE.git
```

Create conda environment:

```
cd AGaVE

conda env create -f agave.yml
```

## Usage

### Activating the environment

This is neccesary every time you want to run the pipeline.

```
conda activate agave
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

The viral database used in this case is the 28 version (Nov, 2023) of the RVDB clustered database release, this is a Reference Viral Database developed for enhancing virus detection using high-throughput/next-generation sequencing (HTS/NGS) technologies, and it can be found in this [link](https://rvdb.dbi.udel.edu/download/C-RVDBv28.0.fasta.gz).

Once you download the fasta.gz file you can unzip the file with the next command:

```
gzip -d file.fasta.gz
```

You can index the database with minimap2 with this command:

```
minimap2 -d /path/to/viral/index.mmi file.fasta
```

Where `file.fasta` is the fasta file containing the reference sequences and `/path/to/viral/index.mm` is the final location of the indexed database. Please note how the reference database can be easily modified and replaced just by indexing a different fasta file that may be more relevant to you. 

### Config file

The `config.yml` file should be modified to have the proper database directories and number of cores to be used.
```
genome_dir : "/path/to/genomeDir"
virus_db : "/path/to/viral/index.mmi"
human_minimap : "/path/to/genomic/index.mmi"
cores : "8"
kmers : "31,53,75,91"
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
snakemake -np sample_log.a sample_int_log.a
```

Running the pipeline with a maximum of 8 cores:

```
snakemake --cores 8 --use-conda --conda-frontend conda sample_log.a sample_int_log.a
```

## Results

A new directory named results will be created after running the pipeline. It will contain two subfolders named contig and cluster, which contain the results for the assembles contigs and the clustered contigs respectively and have similar file organization. The files found in these directories are:

| File | Format | Description  |
| ---- | --------- | ------------ |
| sample_all.fasta | fasta | Contains all the coding sequences found in all assembled contigs or clustered sequences |
| sample_all_aa.fasta | protein fasta | Contains the translation of all the coding sequences found in all assembled contigs or clustered sequences |
| sample_annotated.fasta | fasta | Contains all the coding sequences found in all assembled contigs or clustered sequences and that where annotated to some known gene or protein product |
| sample_annotated_aa.fasta | protein fasta | Contains the translation of all the coding sequences found in all assembled contigs or clustered sequences and that where annotated to some known gene or protein product |
| sample_annotated.csv | csv | Contains the gene identity for the annotated CDSs and shows if the LxCxE motif is present |
| sample_ref_genome_maps.csv | csv | Shows if any sequence was mapped to the viral reference database, to which reference sequence and gives a coverage index |
| sample_aln.csv | csv | Gives information about which contigs or cluster sequences mapped to the reference virus database and in which genomes and positions |
| sample_int_aln.csv | csv | Gives information about which contigs or cluster sequences mapped to the host genome and if they code for any protein |
| sample_int_id.csv | csv | File neccesary for the plotting function of the contigs or cluster sequences mapped to the host genome |






















