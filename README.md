# StaphScan
StaphScan is a clinical-oriented tool for the genomic surveillance of *Staphylococcus aureus*.

It integrates the following steps:
* Species identification
* Assembly Qcheck
* MLST typing
* *spa* typing
* SCCmec typing
* Capsular typing
* Detection of virulence genes (i.e. PVL)
* Detection of antimicrobial resistance genes (i.e. mecA)
* Detection of clinically-relevant mutations (i.e. involved in AMR development)

## Recommendation
⚠️ StaphSCAN is currently in its beta stage. ⚠️ <br>

During this phase, you may encounter bugs, have access to a limited set of features, and experience frequent updates and changes.. <br>

Despite this, we strongly encourage you to use it, helping to accelerate its improvement and its transition to a stable release. <br>

A more structured and comprehensive documentation page will be also deployed soon

## Requirements
StaphSCAN has been built to optimize dependecies. 
It requires:
* python (v3.10 or greater)
* mash
* blast

And the following packages:
* pandas
* biopython

## Usage

```bash
conda create -n staphscan -c bioconda python=3.10 mash blast biopython pandas -y

conda activate staphscan

git clone https://github.com/riccabolla/StaphSCAN.git

cd StaphSCAN/

python main.py -h
```
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
                        Options: assembly, mlst, spa, capsule, sccmec, agr, virulence, resistance
  --list-modules        Generate the list of available modules                        
  --polish              Generate a simplified, clinical-style summary report (recommended)
  --version             Print current version and exit
```
