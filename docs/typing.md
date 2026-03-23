# Typing

StaphSCAN includes several typing methods, each of which can be run as a stand-alone module.

## Multi Locus Sequence Typing

``-m mlst`` 

All genomes identified as *Staphylococcus aureus* are subject to MLST using the seven-locus typing scheme described [here](https://pubmlst.org/organisms/staphylococcus-aureus). When you use this module, please remember to cite [Jolley et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6192448/). 

A copy of the MLST alleles and ST definitions is stored in the ``/data`` directory of this module.

!!! warning
    Due to new PubMLST policies (read [this](https://pubmlst.org/terms-conditions)) it is no more possible to redistribute data published after 31/12/2024. For this reason, all the MLST data bundled with StaphSCAN are up to that date. To update the db follow the guide in the [installation section](installation.md).

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

### Output

* If the identified repeat patterns match a known spa type, the corresponding type is reported.

* Patterns not present in the reference database were classified as novel and reported together with their repeat composition.

* Assemblies in which an X-region was amplified but no known repeat units were detected were reported as “Unknown".

### Limitation

Spa typing is dependent on genome assembly quality, and fragmentation or sequencing errors within the spa X-region may result in spa-negative or Unknown calls. Novel or divergent repeat patterns not present in the reference database are reported as Novel. As with all in silico typing approaches, results may differ from laboratory-based spa typing in cases of mixed populations or incomplete assemblies.

## Accessory Gene Regulator (agr) Typing

``-m agr``

The agr module identifies the *Staphylococcus aureus* accessory gene regulator (agr) type, a quorum-sensing system involved in virulence regulation and commonly classified into four major groups (I–IV) ([Raghuram V et al. 2022](https://doi.org/10.1128/spectrum.01334-21)).

It is an adaption of the tool [agrVATE](https://github.com/VishnuRaghuram94/AgrVATE/tree/main).

A curated set of reference sequences is bundled with the module and stored in the ``/data`` directory.

This module tries to

1) Identify agr type: 
    
The assembly is queried against group-specific probes in `targets.fasta` file. The *agr* type is assigned to the group with the highest count of unique matching probes.

2) Evaluate operon functionality: 

The full *agr* operon is extracted from the assembly using the specific reference sequence (`.gbk`) for the identified group

### Output

Agr group identifiers are mapped to standard agr types as follows:

| Internal ID | Reported agr type |
|-------------|-------------------|
| `gp1`         | `agr I`             |
| `gp2`         | `agr II`            |
| `gp3`         | `agr III`           |
| `gp4`         | `agr IV`            |

The agr module reports:

| Field  | Description|
|--------|------------|
|`agr_type`| Assigned agr group (agr I–IV)|
|`agr_ confidence`| N of probes matching the assigned group |
|`agr_frameshifts`|Report eventual gene defections|
|`agr_operon_status`|Operon functionality assessment|

The following criteria are used to report the `agr_operon-status`:

* `Intact`: *agrC* and *agrA* coding sequences are complete and functional.

* `Pseudogene`: Frameshifts, or premature stop codons detected in *agrC* or *agrA*.

* `Assembly Gap`: The gene appears truncated but ends precisely at the edge of a contig. This suggests the gene might be intact but wasn't fully assembled, rather than a biological mutation.

* `Missing/Fragmented`: The operon could not be extracted or is too fragmented to analyze.

* `Ref Missing`: Reference data for the identified group is unavailable.

### Notes and limitations

* While the module attempts to distinguish assembly breaks (`Assembly Gap`) from true mutations (`Pseudogene`), highly fragmented assemblies may still result in ambiguous functionality calls. For this reason, positive `Pseudogene` results, should be treated with caution, and investigate properly using a read-based method (i.e. [Snippy](https://github.com/tseemann/snippy))

* Typing uses relaxed BLAST parameters (90% identity) to correctly classify divergent lineages, but relies on a strict count of unique probes to ensure specificity.

* The module is tuned to ignore natural allelic variation (SNPs) while catching structural defects (indels/stops) that destroy protein function.

## Capsule Typing

``-m capsule``

The capsule module identifies the *Staphylococcus aureus* capsular polysaccharide operon and assigns the predominant capsule serotype (Type 5 or Type 8). Capsular polysaccharides are major virulence determinants involved in immune evasion and are encoded by the cap operon (capA–P), with serotype specificity driven by the H–K loci ([Cocchiaro et al. 2006](https://doi.org/10.1111/j.1365-2958.2005.04978.x)).

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
| `cap_type`         | Assigned capsule serotype (Type 5, Type 8, or -)   |
| `cap_completeness` | Capsule operon status (Complete, Incomplete, or -) |
| `cap_genes`        | Semicolon-separated list of detected capsule genes |

## SCCmec Typing 

``-m sccmec``

The sccmec module detects and classifies *Staphylococcus aureus* SCCmec elements, which carry methicillin resistance determinants and are defined by combinations of the mec gene complex and ccr recombinase genes. 

It is adapted from the tool [sccmec](https://github.com/rpetit3/sccmec).

A curated FASTA file containing representative SCCmec-associated target genes is bundled with the module and stored in the ``/data`` directory.

The module supports the classification of the following `types`:

| Type | Reference |
|------|----------|
| I    | [Katayama et al. 2000](https://doi.org/10.1128/aac.44.6.1549-1555.2000) |
| II   | [Katayama et al. 2000](https://doi.org/10.1128/aac.44.6.1549-1555.2000), [Ito et al. 2001](https://doi.org/10.1128%2FAAC.45.5.1323-1336.2001) |
| III  | [Katayama et al. 2000](https://doi.org/10.1128/aac.44.6.1549-1555.2000) |
| IV   | [Ma et al. 2002](https://doi.org/10.1128%2FAAC.46.4.1147-1152.2002) |
| V    | [Ito et al. 2004](https://doi.org/10.1128/aac.48.7.2637-2651.2004) |
| VI   | [Oliveira et al. 2006](https://doi.org/10.1128%2FAAC.00629-06) |
| VII  | [Berglund et al. 2008](https://doi.org/10.1128%2FAAC.00087-08) |
| VIII | [Zhang et al. 2009](https://doi.org/10.1128%2FAAC.01118-08) |
| IX   | [Li et al. 2011](https://doi.org/10.1128%2FAAC.01475-10) |
| X    | [Li et al. 2011](https://doi.org/10.1128%2FAAC.01475-10) |
| XI   | [García-Álvarez et al. 2011](https://doi.org/10.1128/aac.01692-15) |
| XII  | [Wu et al. 2015](https://doi.org/10.1128/JCM.39.2.607-612.2001) |
| XIII | [Baig et al. 2018](https://doi.org/10.1016/j.meegid.2018.03.013) |
| XIV  | [Urushibara et al. 2020](https://doi.org/10.1093/jac/dkz406) |
| XV   | [Wang et al. 2022](https://doi.org/10.1093/jac/dkab500) |

And of the following `subtype`: 

| SubType | Reference |
|------|----------|
| Ia   | [Ito et al. 2001](https://doi.org/10.1128/AAC.45.5.1323-1336.2001) |
| Ib   | [Han et al. 2009](https://doi.org/10.1128/AAC.00772-08), [Oliveira et.al. 2006](https://doi.org/10.1093/jac/dkl208) |
| IIa  | [Katayama et al. 2000](https://doi.org/10.1128/aac.44.6.1549-1555.2000), [Ito et al. 2001](https://doi.org/10.1128/AAC.45.5.1323-1336.2001) |
| IIb  | [Hisata et al. 2005](https://doi.org/10.1128/JCM.43.7.3364-3372.2005) |
| IIc  | [Shore et al. 2005](https://doi.org/10.1128/AAC.49.5.2070-2083.2005) |
| IId  | [Kondp et al. 2007](https://doi.org/10.1128/AAC.00165-06) |
| IIe  | [Han et al. 2009](https://doi.org/10.1128/AAC.00772-08) |
| IVa  | [Ma et al. 2002](https://doi.org/10.1128/AAC.46.4.1147-1152.2002) |
| IVb  | [Ma et al. 2002](https://doi.org/10.1128/AAC.46.4.1147-1152.2002) |
| IVc  | [Ma et al. 2006](https://doi.org/10.1128/JCM.00985-06) |
| IVd  | [Ma et al. 2006](https://doi.org/10.1128/JCM.00985-06) |
| IVg  | [Kwon et al. 2005](https://doi.org/10.1093/jac/dki306) |
| IVh  | [Milheirico et al. 2007](https://doi.org/10.1093/jac/dkm112) |
| IVi  | [Berglund et al. 2009](https://doi.org/10.1093/jac/dkn435) |
| IVj  | [Berglund et al. 2009](https://doi.org/10.1093/jac/dkn435)  |
| IVk  | - |
| IVl  | [Iwao et al. 2012](https://doi.org/10.1007/s10156-012-0379-6) |
| IVm  | [Hosoya et al. 2014](https://doi.org/10.11150/kansenshogakuzasshi.88.840) |
| IVn  | - |
| Va   | [Ito et al. 2004](https://doi.org/10.1128/AAC.48.7.2637-2651.2004) |
| Vb   | [Hisata et al. 2011](https://doi.org/10.1007/s10156-011-0223-4) |
| Vc   | [Li et al. 2011](https://doi.org/10.1128/AAC.01475-10) |

### Parameters

Hits are filtered based on the following parameters:

* Minimum alignment percentage of 90
* Minimum coverage percentage of 80

### Output

The mec gene complex is classified as follows:

| mec class | Required components                 |
| --------- | ----------------------------------- |
| A         | mecI + mecR1 + mecA                 |
| B         | IS1272 + mecA                       |
| C         | IS431 + mecA                        |
| Unknown   | mecA or mecC present but incomplete |
| None      | No mec genes detected               |

Both mecA and mecC are supported. 

Detected recombinase complexes include:

| ccr complex | Required genes |
| ----------- | -------------- |
| 1           | ccrA1 + ccrB1  |
| 2           | ccrA2 + ccrB2  |
| 3           | ccrA3 + ccrB3  |
| 4           | ccrA4 + ccrB4  |
| C1          | ccrC1          |
| C2          | ccrC2          |
| A1B6        | ccrA1 + ccrB6  |
| A1B3        | ccrA1 + ccrB3  |

SCCmec types are inferred by combining the detected mec class and ccr complex(es), following established nomenclature where possible:

| mec class | ccr complex | Assigned type      |
| --------- | ----------- | ------------------ |
| B         | 1           | Type I (1B)        |
| A         | 2           | Type II (2A)       |
| B         | 2           | Type IV (2B)       |
| A         | 3           | Type III (3A)      |
| B         | 4           | Type VI (4B)       |
| A         | 4           | Type VIII (4A)     |
| C         | C1 (5)      | Type V (5C)        |
| C + IS12960D         | C1 (5)      | Type VII (5C + IS12960D)        |
| C         | 1           | Type IX (1C)       |
| C/B       | A1B6 (7)    | Type X (A1B6)      |
| A/E       | A1B3 (8)    | Type XI (mecC-associated)       |
| C         | C2 (9)      | Type XII (9C)      |
| A         | C2 (9)      | Type XIII (9A)     |
| A         | C1 (5)      | Type XIV (5A)      | 
| A         | A1B6 (7)    | Type XV (A1B6)     |

If multiple compatible SCCmec types are detected, a Composite SCCmec assignment is reported.

If mec genes are detected but no ccr genes are found, the element is reported as an "orphan" cassette.

The module reports:

| Field          | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `sccmec_type `   | Assigned SCCmec type                                         |
| `sccmec_subtype` | Assigned subtype                                             |
| `sccmec_genes`   | Semicolon-separated list of detected SCCmec-associated genes |