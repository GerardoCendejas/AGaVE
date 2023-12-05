# (Name)


## Introduction

(Name) is a pipeline designed for the search of viral sequences in RNAseq samples from human tissues or human related samples. This pipeline was designed used previously developed tools that are available via conda and organized using the snakemake workflow management system. A basic diagram of the pipeline is shown in the next image.

<p align="center" width="100%">
    <img width="50%" src="https://github.com/GerardoCMM/RNAseq-Retinoblastoma/blob/main/dag.svg"> 
</p>

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

>Your samples should be located at the samples/ directory in fastq format and the naming:
> - *sample*_1.fastq
> - *sample*_2.fastq
> where *sample* can be any given name for your sample.


Activating the environment:

```
conda activate rna_seq_analysis
```

### For running the full pipeline

Checking the steps to be run:

```
snakemake -np *sample*_log.a humanmap/*sample*_mapped2human.fastq
```

Running the pipeline with a maximum of 8 cores

```
snakemake --cores 8 --use-conda *sample*_log.a humanmap/*sample*_mapped2human.fastq
```