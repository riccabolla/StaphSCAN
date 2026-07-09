# StaphSCAN
[![Anaconda-Server Badge](https://anaconda.org/bioconda/staphscan/badges/version.svg)](https://anaconda.org/bioconda/staphscan)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/staphscan/badges/latest_release_date.svg)](https://anaconda.org/bioconda/staphscan)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/staphscan/badges/license.svg)](https://anaconda.org/bioconda/staphscan)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/staphscan)](https://anaconda.org/bioconda/staphscan)

StaphSCAN is a genome-based tool for **S**urveillance through **C**omprehensive **A**nalysis and sta**N**dardized reporting of ***Staph**ylococcus aureus*.

It integrates the following steps:
* Species identification
* Assembly Qcheck
* MLST typing
* *spa* typing
* SCCmec typing
* Capsular typing
* Detection of virulence genes (i.e. PVL)
* Detection of biofilm-related genes
* Detection of antimicrobial resistance genes (i.e. mecA)
* Detection of clinically-relevant mutations (i.e. involved in AMR development)

<br>

![StaphSCAN workflow](docs/images/pipeline.svg)

## Documentation

A documentation is available [here](https://staphscan.readthedocs.io/)

## Citation

You can cite StaphSCAN using:

> Bollini, R., Cento, V. StaphSCAN: a genomic surveillance framework for Staphylococcus aureus. 
> BMC Genomics (2026. https://doi.org/10.1186/s12864-026-13165-0


## Requirements
StaphSCAN has been built to optimize dependecies. 
It requires:
* python (v3.10 or greater)
* mash
* blast

And the following packages:
* pandas
* biopython
* rauth
* requests

## Usage

```bash
conda create -n staphscan -c bioconda staphscan -y

conda activate staphscan

staphscan -h
```
### Basic usage

```bash
staphscan -i /path/to/genomes/*.fasta -o results_directory
```
### Options

```bash
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input FASTA files (supports *.fasta)
  -o OUTDIR, --outdir OUTDIR
                        Output directory
  -m MODULES, --modules MODULES
                        Comma-separated list of modules to run (default: "all")
  --list-modules        Generate the list of available modules
  --mlst_update         Update mlst database
  --db_mlst             Custom folder for mlst db                        
  --report REPORT_NAME  Generate a report with a custom name
  --version             Print current version and exit
```
