# Usage

## Input files

StaphSCAN accepts two input formats:

* **Genome assemblies (FASTA)**: `.fasta`, `.fa`, `.fna` (can be either draft or completed assemblies; completed assemblies reduce the risk of fragmented genes/loci).
* **Raw sequencing reads (FASTQ)**: uncompressed or gzip-compressed FASTQ files (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`) from Illumina short-read (`sr`), Oxford Nanopore (`ont`), or PacBio HiFi (`pb`) platforms.

## Basic usage

### FASTA assemblies

```python
staphscan -i *.fasta -o staphscan_output
```
`-i *.fasta` : Specifies the input files (assemblies) to be analyzed (.fasta). 

`-o` : Specifies the directory where the output files will be saved. 

### FASTQ raw reads

**Auto-paired paired-end short reads via `-i`:**

```bash
# Pair R1 and R2 files for a single sample
staphscan -i sample1_R1.fastq.gz sample1_R2.fastq.gz -o staphscan_output

# Batch process multiple samples (or mix FASTA and FASTQ inputs in the same run)
staphscan -i *.fastq.gz -o staphscan_output
```

When FASTQ files are passed to `-i`, StaphSCAN automatically pairs forward and reverse reads by sample name (stripping extensions and trailing `_R1`/`_R2` or `_1`/`_2`, with optional `_001`). You can also mix FASTA and FASTQ files in the same `-i` invocation.

**Explicit paired-end or single-end short reads via `-1` / `-2`:**

```bash
# Paired-end short reads
staphscan -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o staphscan_output

# Single-end short reads
staphscan -1 sample.fastq.gz -o staphscan_output
```

**Long reads via `-lr` / `--long_reads`:**

```bash
# Oxford Nanopore reads (default read type for -lr)
staphscan -lr sample_ont.fastq.gz -o staphscan_output

# PacBio HiFi reads
staphscan -lr sample_hifi.fastq.gz --read_type pb -o staphscan_output

# Multiple long-read samples
staphscan -lr sample1_ont.fastq.gz sample2_ont.fastq.gz -o staphscan_output
```

**Custom consensus depth and thread count:**

```bash
staphscan -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --min_depth 10 -t 8 -o staphscan_output
```

## FASTQ workflow

When FASTQ files are provided, StaphSCAN runs an integrated read-based pipeline before executing downstream typing and characterization modules:

1. **Species identification & read QC (`assembly` module)**: Reads are streamed into `mash sketch -r` with adaptive minimum k-mer copy filtering (`-m 4` up to `-m 6` based on estimated coverage) and compared with `mash dist` against the StaphSCAN reference database to estimate species, Mash distance, genome size, sequencing coverage, and intra-genus contamination (`Mixed`). Samples not identified as *S. aureus* skip downstream analyses.
2. **Reference-guided consensus generation**: StaphSCAN builds a non-redundant combined target reference for active downstream modules (`mlst`, `capsule`, `agr`, `biofilm`, `virulence`, `resistance`, `sccmec`), aligns reads with `minimap2` (using preset `sr`, `map-ont`, or `map-hifi` according to `--read_type`), sorts alignments with `samtools sort`, and generates consensus sequences with `samtools consensus` at `--min_depth` (default: `5`).
3. **Targeted locus reconstruction & accuracy safeguards**:
    * Performs an insertion-aware pileup check (`--show-ins yes`) to detect genuine high-depth frameshift truncations, preserves in-frame whole-codon deletions, and applies strictly capped (≤ 10% of gene length) reference healing only when unhealed high-depth coverage is ≥ 83% (excluding `cna` and `agr` operons).
    * For short-read (`sr`) inputs, reconstructs hypervariable biofilm MSCRAMM genes (`clfA`, `clfB`, `fnbA`, `fnbB`) using two-round baited local multi-k de Bruijn assembly with paired-end repeat resolution (bridging shared internal paralog repeats such as `fnbA`/`fnbB`) to avoid reference-template bias, and estimates `cna` B-repeat copy number from normalized B-domain-to-A-domain aligned read depth (`3 × mean_B_depth / mean_A_depth`, masking B2–B3 when `< 2.5` copies).
    * Maps reads against SCCmec reference regions to recover the best-covered region for SCCmec subtyping.
    * For `spa` typing, extracts primer-spanning reads directly when available (e.g., long reads) or baits repeat/primer k-mers and walks a local de Bruijn graph across the four standard `spa` primer sets to reconstruct the `spa` amplicon.
4. **Standardized reporting**: The temporary consensus FASTA is passed to all selected downstream modules, producing the exact same `.tsv` summary report format as FASTA mode.

!!! note
    `minimap2` and `samtools` (v1.16 or later) must be installed and available in your `PATH` when analyzing FASTQ inputs.

## Modules selection

Module can be run separately, by using the `-m` parameter:

```python
# select one module

staphscan -i *.fasta -o staphscan_output -m mlst

# select more modules

staphscan -i *.fasta -o staphscan_output -m mlst,spa
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

Output files are tab-delimited (.tsv) files. Columns included in each output file will depend on the modules that are run.

To use a custom report name:

```python
staphscan --report custom_name
```

## Parameters

**Input/Output**

`-i INPUT [INPUT ...], --input INPUT [INPUT ...]`<br>
Input FASTA or FASTQ file(s) (paired-end FASTQ files are auto-paired by sample name)
<br>

`-1 READ1, --read1 READ1`<br>
Forward FASTQ read (R1) or single-end short reads
<br>

`-2 READ2, --read2 READ2`<br>
Reverse FASTQ read (R2) for paired-end reads
<br>

`-lr LONG_READS [LONG_READS ...], --long_reads LONG_READS [LONG_READS ...]`<br>
Single-end long-read FASTQ file(s) (e.g., Nanopore/PacBio). Defaults `--read_type` to `ont` unless `pb` is set
<br>

`--read_type {sr,ont,pb}`<br>
Read type for FASTQ mapping: `sr` (Illumina), `ont` (Nanopore), or `pb` (PacBio) (default: `sr`)
<br>

`--min_depth MIN_DEPTH`<br>
Minimum read depth to call a base in FASTQ consensus (default: `5`)
<br>

`-t THREADS, --threads THREADS`<br>
Number of threads for FASTQ mapping (default: `4`)
<br>

`-o OUTDIR, --outdir OUTDIR`<br>
Directory for storing output files
<br>

`--report REPORT`<br>
Custom filename for the output report (default: `staphscan_summary.tsv`)

**Modules**

`--list_modules` <br>
Print a list of all available modules and then quit (default: False)

`-m MODULES, --modules MODULES`: <br>
Comma-delimited list of staphscan modules to use (default: all)

`--mlst_update` <br>
Download latest database from PubMLST

`--db_mlst` <br>
Location of the alternative folder containing the mlst database (useful when working in read-only env)

Module-specific parameters are available in the [Modules](modules.md) section.

**Help**

`-h, --help` <br>
Show a help message and exit <br>


`--version` <br>
Show program's version number and exit

