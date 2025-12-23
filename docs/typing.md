# Typing

StaphSCAN includes several typing methods, each of which can be run as a stand-alone module.

For clarity, all these methods and their respective parameters are presented together on this page

## Multi Locus Sequence Typing

``-m mlst`` 

All genomes identified as *Staphylococcus aureus* are subject to MLST using the seven-locus typing scheme described [here](https://pubmlst.org/organisms/staphylococcus-aureus) 

A copy of the MLST alleles and ST definitions is stored in the ``/data`` directory of this module.

### Parameters

Locus detection is filtered by minimum alignment identity and coverage:

``--min_id_mlst `` : Minimum alignment percentage identity (default: 95)

``--min_cov_mlst `` : Minimum alignment percent coverage (default: 95)

### Output

Both the ST and the allelic profile are reported. 

For each locus, the following annotations may be reported:

* Exact matches: allele sequences that exactly match a known MLST allele are reported using the corresponding allele number.
* Putative novel alleles: loci with full-length sequences that do not exactly match any known allele are reported as the closest known allele followed by a `` * ``
* Partial loci: loci detected but not covering the full reference length are reported as Partial.
* Missing loci: loci not detected in the assembly are reported as `` - ``.

Imprecise or incomplete allelic profiles result in approximate ST assignments. In these cases, StaphSCAN reports the closest matching ST followed by the number of differing loci (n-locus variants, up to two). Example: `` ST1-1LV `` (closest match is ST1 with one differing allele)

## spa typing

``-m spa``

The spa-typing is a method based on the characterization of the repeat regions of Staphylococcus protein A gene (spa).This method is widely used for rapid typing of MRSA, particularly in hospital and surveillance settings.

For more information visit [here](https://spa.ridom.de/index.shtml). 

A local copy of [Ridom database](https://spa.ridom.de/spatypes.shtml) is distributed with this module and stored the modules's ``/data`` directory.

Genome assemblies are screened for the presence of the *spa* gene X-region by simulating PCR amplification with multiple published primer sets. Each primer set is tested against all contigs, and both forward and reverse-complement orientations are evaluated. 

When multiple primer sets yield a valid amplicon, the first detected product is used for downstream analysis. Assemblies in which no valid amplicon is detected are reported as spa-negative.

The amplified X-region is scanned to identify spa repeat units using a curated database of known repeat sequences. Detected repeats are recorded sequentially to generate a repeat pattern, which is then compared against the reference spa type database.

### Output

* If the identified repeat patterns match a known spa type, the corresponding type is reported.

* Patterns not present in the reference database were classified as novel and reported together with their repeat composition.

* Assemblies in which an X-region was amplified but no known repeat units were detected were reported as “Unknown".

### Limitation

Spa typing is dependent on genome assembly quality, and fragmentation or sequencing errors within the spa X-region may result in spa-negative or Unknown calls. Novel or divergent repeat patterns not present in the reference database are reported as Novel. As with all in silico typing approaches, results may differ from laboratory-based spa typing in cases of mixed populations or incomplete assemblies.

## Accessory Gene Regulator (agr) Typing

``-m agr``

The agr module identifies the Staphylococcus aureus accessory gene regulator (agr) type, a quorum-sensing system involved in virulence regulation and commonly classified into four major groups (I–IV).

A curated FASTA file containing representative agr group target sequences is bundled with the module and stored in the ``/data`` directory.

### Typing Strategy

Agr typing is performed using nucleotide BLAST against a local database of agr target sequences.

Hits are filtered by minimum sequence identity (90.0%).

Candidate matches are ranked by:

* Percent identity

* Alignment coverage

* BLAST bitscore

The top-ranking hit is used to assign the agr group.

Agr group identifiers are mapped to standard agr types as follows:

| Internal ID | Reported agr type |
|-------------|-------------------|
| gp1         | agr I             |
| gp2         | agr II            |
| gp3         | agr III           |
| gp4         | agr IV            |


### Output

The agr module reports:

| Field               | Description                                  |
|---------------------|----------------------------------------------|
| agr_type            | Assigned agr group (agr I–IV)                 |
| agr_match_confidence| Percent identity of the best BLAST hit        |


* agr I–IV

    A valid agr group was identified and reported.

* Negative

    No BLAST hits passed the identity threshold.

* Error
 
    An unexpected error occurred during analysis.

### Notes and limitations
* Only the top-ranked BLAST hit is used for agr assignment.
* Assemblies with fragmented or highly divergent agr loci may yield Negative results.
* Reported confidence reflects percent identity, not overall locus completeness.
* The module assumes one dominant agr group per genome.