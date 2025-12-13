# Species detection & contig stats


``-m assembly``


This module will attempt to identify the species of each input assembly, and generate some basic assembly statistics. 

## Species identification


Species identification is performed using Mash (Ondov et al., 2016) to compare the input assembly against a database of reference genomes.
The database includes high-quality complete genomes of several *Staphylococcus* species downloaded from NCBI RefSeq.

Currently, it uses two thresholds:

* Mash distance threshold for a strong species match (default: 0.02)

* Mash distance threshold for a weak species match (default: 0.04)

To report the species, the following criteria are used:

* The species with the lowest distance is reported. 
* If *Staphylococcus aureus* is identified with a distance below the strong threshold, it is reported as a "Strong match". 
* If *Staphylococcus aureus* is identified with a distance between the strong and weak threshold, it is reported as a "Weak match".
* If *Staphyloccocus aureus* is not the best match, but one of the listed genomes is, it is reported.
* If the lowest Mash distance is > 0.04, the result is reported as “No match found”.

## Contigs stats

For assembly quality the following parameters are considered:

* Total assembly size (compared to expected size for *Staphylococcus aureus* of 2.6 - 3.1 Mbp)
* N50 (>=10 kbp)
* Presence of ambiguous bases (Ns)

## Outputs

The assembly module generates the following output columns in the default report:

* **Species**: <br>
Detected species based on Mash distance to reference genomes (the strongest match is considered). 

* **Total_size**: <br>
Total assembly size (in bp).

* **QC**: <br>
Overall QC status of the assembly: `PASS`  or `FAILED` . 

All results from failed QC checks should be treated with caution.
