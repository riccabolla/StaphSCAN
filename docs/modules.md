# Modules

StaphSCAN includes a range of modules for typing bacterial genomes. 
These modules will be run by default, unless specified otherwise by the user (using the `-m` parameter). All modules support both FASTA assemblies and FASTQ raw reads (via read-based QC in `assembly` and reference-guided consensus reconstruction for downstream modules).

Here a list of the available modules and their functions:

* **assembly**: Quality control of genome assemblies or raw sequencing reads (FASTQ), including species identification
* **mlst**: Multi-locus sequence typing (MLST) of *Staphylococcus aureus*
* **spa**: spa typing of *Staphylococcus aureus*
* **capsule**: Capsular polysaccharide typing of *Staphylococcus aureus*
* **sccmec**: SCCmec typing of *Staphylococcus aureus*
* **agr**: Accessory gene regulator (agr) typing of *Staphylococcus aureus*
* **resistance**: Detection of antimicrobial resistance genes (including clinically relevant mutations)
* **virulence**: Detection of virulence genes (including clinically relevant mutations)
* **biofilm**: Detection of biofilm associated genes

For more information on each module, see the respective sections in this documentation.