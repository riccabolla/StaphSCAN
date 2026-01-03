# Typing

StaphSCAN includes several typing methods, each of which can be run as a stand-alone module.

## Multi Locus Sequence Typing

``-m mlst`` 

All genomes identified as *Staphylococcus aureus* are subject to MLST using the seven-locus typing scheme described [here](https://pubmlst.org/organisms/staphylococcus-aureus) (1). 

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

### Citations

1) Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018;3:124. Published 2018 Sep 24. doi:10.12688/wellcomeopenres.14826.1

## spa typing

``-m spa``

The spa-typing is a method based on the characterization of the repeat regions of Staphylococcus protein A gene (spa).This method is widely used for rapid typing of MRSA, particularly in hospital and surveillance settings.

For more information visit [here](https://spa.ridom.de/index.shtml). 

A local copy of [Ridom database](https://spa.ridom.de/spatypes.shtml) is distributed with this module and stored the modules's ``/data`` directory.

### Output

* If the identified repeat patterns match a known spa type, the corresponding type is reported.

* Patterns not present in the reference database were classified as novel and reported together with their repeat composition.

* Assemblies in which an X-region was amplified but no known repeat units were detected were reported as “Unknown".

### Limitation

Spa typing is dependent on genome assembly quality, and fragmentation or sequencing errors within the spa X-region may result in spa-negative or Unknown calls. Novel or divergent repeat patterns not present in the reference database are reported as Novel. As with all in silico typing approaches, results may differ from laboratory-based spa typing in cases of mixed populations or incomplete assemblies.

## Accessory Gene Regulator (agr) Typing

``-m agr``

The agr module identifies the *Staphylococcus aureus* accessory gene regulator (agr) type, a quorum-sensing system involved in virulence regulation and commonly classified into four major groups (I–IV) (1).

A curated FASTA file containing representative agr group target sequences is bundled with the module and stored in the ``/data`` directory.


### Output

Agr group identifiers are mapped to standard agr types as follows:

| Internal ID | Reported agr type |
|-------------|-------------------|
| gp1         | agr I             |
| gp2         | agr II            |
| gp3         | agr III           |
| gp4         | agr IV            |

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

### Citations

1) Raghuram V, Alexander AM, Loo HQ, Petit RA, Goldberg JB, Read TD.2022.Species-Wide Phylogenomics of the Staphylococcus aureus Agr Operon Revealed Convergent Evolution of Frameshift Mutations. Microbiol Spectr10:e01334-21.https://doi.org/10.1128/spectrum.01334-21

## Capsule

``-m capsule``

The capsule module identifies the *Staphylococcus aureus* capsular polysaccharide operon and assigns the predominant capsule serotype (Type 5 or Type 8). Capsular polysaccharides are major virulence determinants involved in immune evasion and are encoded by the cap operon (capA–P), with serotype specificity driven by the H–K loci (1).

A curated FASTA file containing representative capsule gene sequences is bundled with the module and stored in the ``/data`` directory.

### Parameters

Hits are filtered by minimum alignment identity and coverage:

``--min_id_capsule`` : Minimum alignment percentage identity (default: 90)

``--min_cov_capsule``: Minimum alignment percentage coverage (default: 80)

### Output

Capsule serotype is inferred based on the presence of serotype-specific loci:

* **Type 5**: cap5H, cap5I, cap5J, cap5K
* **Type 8**: cap8H, cap8I, cap8J, cap8K

Once a serotype is assigned, operon completeness is evaluated by checking for the presence of all expected genes. 

The operon is classified as:

* **Complete** : all genes detected
* **Incomplete** : at least one gene missing

| Field            | Description                                        |
| ---------------- | -------------------------------------------------- |
| cap_type         | Assigned capsule serotype (Type 5, Type 8, or -)   |
| cap_completeness | Capsule operon status (Complete, Incomplete, or -) |
| cap_genes        | Semicolon-separated list of detected capsule genes |

### Citations

1) Cocchiaro, J.L., Gomez, M.I., Risley, A., Solinga, R., Sordelli, D.O. and Lee, J.C. (2006), Molecular characterization of the capsule locus from non-typeable Staphylococcus aureus. Molecular Microbiology, 59: 948-960. https://doi.org/10.1111/j.1365-2958.2005.04978.x

