# Usage

## Input files

Genome assemblies in FASTA format. 
Can be either draft or completed assemblies (completed is better because it reduces the risk of fragmented genes/loci).

## Basic usage

```python
staphscan -i *.fasta -o staphscan_output --polish
```
`-i *.fasta` : Specifies the input files (assemblies) to be analyzed (.fasta). 


`-o` : Specifies the directory where the output files will be saved. 


`--complete` : Specifies the output format (i.e. detailed output)

## Modules selection

Module can be run separately, by using the `-m` parameter:

```python
# select one module

staphscan -i *.fasta -o staphscan_output --polish -m mlst

# select more modules

staphscan -i *.fasta -o staphscan_output --polish -m mlst,spa
```
To see all modules:
```python
staphscan --list_modules
```

Check available modules, check version, print help:

```python
staphscan [--list_modules] [--version] [-h]
```

## Output files

Output files are tab-delimited (.tsv) files. Columns included in each output file will depend on the modules that are run, and on the report version selected (simplified or detailed).

## Parameters

**Input/Output**

`-i ASSEMBLIES [ASSEMBLIES ...], --input ASSEMBLIES [ASSEMBLIES ...]` <br>
FASTA file(s) for assemblies
<br>

`-o OUTDIR, --outdir OUTDIR` <br>
Directory for storing output files (default: staphscan_results)
<br>

`--complete` <br>
Create a complete final report based on the modules run

**Modules**

`--list_modules` <br>
Print a list of all available modules and then quit (default: False)

`-m MODULES, --modules MODULES`: <br>

Comma-delimited list of staphscan modules to use (default: all)

**Help**

`-h, --help`: Show a help message and exit <br>


`-v, --version`: Show program's version number and exit

