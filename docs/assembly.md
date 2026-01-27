# Species detection & contig stats


``-m assembly``


This module will attempt to identify the species of each input assembly, and generate some basic assembly statistics. 

## Species identification


Species identification is performed using Mash (Ondov et al., 2016) to compare the input assembly against a database of reference genomes.
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

!!! note
    All genomes identified as not *Staphylococcus aureus* are skipped for downstream analysis. 

## Contigs stats

For assembly quality the following parameters are considered:

* Total assembly size (compared to expected size for *Staphylococcus aureus* of 2.6 - 3.1 Mbp)
* N50 (>=10 kbp)
* Presence of ambiguous bases (Ns)

## Outputs

The assembly module generates the following output columns in the default report:

|Field|Description|
|-----|-----------|
|`Species`|Detected specie|
|`Total_size`|Assembly size (bp)|
|`QC`|Overall QC status (`PASS`  or `FAILED`)|


All results from failed QC checks should be treated with caution.
