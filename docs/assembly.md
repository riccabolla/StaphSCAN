# Species detection & contig stats


``-m assembly``


This module will attempt to identify the species of each input genome assembly (FASTA) or raw sequencing read dataset (FASTQ), and generate basic quality control statistics. 

## Species identification


Species identification is performed using Mash ([Ondov et al., 2016](https://pubmed.ncbi.nlm.nih.gov/27323842/)) to compare the input assembly against a database of reference genomes.
The database includes the following high-quality genomes:

|Specie|Genome|
|------|---------|
|*Staphylococcus argenteus*|[GCF_000236925.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000236925.1/)|
|*Staphylococcus aureus*|[GCF_000013425.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/)|
|*Staphylococcus capitis*|[GCF_040739365.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_040739365.1/)|
|*Staphylococcus epidermidis*|[GCF_006094375.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_006094375.1/)|
|*Staphylococcus haemolyticus*|[GCF_006094395.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_006094395.1/)|
|*Staphylococcus lugdunensis*|[GCF_001558775.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001558775.1/)|
|*Staphylococcus schweitzeri*|[GCF_900636685.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900636685.1/)|


Currently, it uses a mash distance threshold of 0.04. 

To report the species, the following criteria are used:

* The species with the lowest distance is reported, if lower than the threshold. 
* If the lowest Mash distance is > 0.04, the result is reported as “No match found”.

For FASTQ inputs, reads are sketched directly using `mash sketch -r` with adaptive minimum k-mer copy filtering (starting at `-m 4` and scaling up to `-m 6` for higher-coverage datasets) before computing Mash distances against the same reference database using the same `0.04` distance threshold.

!!! note
    All genomes identified as not *Staphylococcus aureus* are skipped for downstream analysis. 

## Assembly stats

Quality control checks depend on whether the input is a genome assembly (FASTA) or raw sequencing reads (FASTQ):

### FASTA inputs

* Total assembly size (compared to expected size for *Staphylococcus aureus* of 2.6 - 3.1 Mbp; flagged `Too short` or `Too long`)
* Contamination: If more than one reference genome returns a Mash distance of **≤ 0.04**, the assembly is flagged (`Mixed`). This indicates an intra-genus mixed culture (e.g., *S. aureus* + *S. epidermidis*).
* N50 (`>= 10 kbp`; otherwise flagged `Low_N50`)
* Presence of ambiguous bases (`Ns`; flagged `Ambiguous_Bases`)

### FASTQ inputs

* Estimated genome size (from Mash): flagged `Too short` if `< 2.5 Mbp`; flagged `Too long` if `> 3.3 Mbp` (for short-read `sr` inputs only).
* Estimated sequencing coverage (from Mash): flagged `Low_Coverage` if `< 10x`.
* Contamination: If more than one reference genome returns a Mash distance of **≤ 0.04**, the sample is flagged (`Mixed`).
* Contig metrics (`N_contig` and `N50`) are not applicable to raw reads and are reported as `-`.

## Outputs

The assembly module generates the following output columns in the report:

|Field|Description|
|-----|-----------|
|`Species`|Detected specie|
|`Mash_distance`|Calculated Mash distance to the best-matching reference genome|
|`Total_size`|Assembly size (bp), or Mash-estimated genome size for FASTQ input|
|`N_contig`|Total number of contigs in the assembly (`-` for FASTQ input)|
|`N50`|Assembly N50 value (`-` for FASTQ input)|
|`QC`|Overall QC status (`PASS`  or `FAILED`)|


All results from failed QC checks should be treated with caution.
