# StaphScan
StaphScan is a tool that allow rapid and standardized genotyping of Staphylococcus aureus

It integrates the following steps:
* Species identification
* MLST typing
* *spa* typing
* SCCmec typing
* Detection of virulence genes (i.e PVL; Exofoliatin; tsst)
* Detection of antimicrobial resistance genes (mecA/C; blaZ; parC mutations)
* Quality control of genome assemblies

## Requirements
StaphSCAN has been built to optimize dependecies. 
For this reason, it only requires:
* python (v3.10 or greater)
* mash
* blast

The following python packages are also required:
* pandas
* biopython

## Installation
Currently, it can only be installed cloning this repo via:
```bash
git clone https://github.com/riccabolla/StaphSCAN.git

cd StaphSCAN/

python main.py -h
```
## Usage

### Basic usage

```bash
python main.py -i /path/to/genomes/*.fasta -o results_directory
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
                        Options: assembly, mlst, spa, sccmec, agr, virulence, resistance
  --polish              Generate a simplified, clinical-style summary report (recommended)
```