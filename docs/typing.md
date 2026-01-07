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

## sccmec 

``-m sccmec``

The sccmec module detects and classifies *Staphylococcus aureus* SCCmec elements, which carry methicillin resistance determinants and are defined by combinations of the mec gene complex and ccr recombinase genes.

A curated FASTA file containing representative SCCmec-associated target genes is bundled with the module and stored in the ``/data`` directory.

The module supports the classification of the following SCC*mec* types:

| Type | Reference |
| :---  | :--- |
| **I** |(1) |
| **II** |(1, 2) |
| **III** |(1) |
| **IV** |(3) |
| **V** |(4) |
| **VI** |(5) |
| **VII** |(6) |
| **VIII** |(7) |
| **IX** | (8) |
| **X** | (8) |
| **XI** | (9) |
| **XII** | (10)|
| **XIII** | (11) |
| **XIV** | (12) |
| **XV** | (13) |

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
| sccmec_type    | Assigned SCCmec type                                         |
| sccmec_genes   | Semicolon-separated list of detected SCCmec-associated genes |

### Citations

1) Katayama Y, Ito T, Hiramatsu K.2000.A New Class of Genetic Element, Staphylococcus Cassette Chromosome mec, Encodes Methicillin Resistance in Staphylococcus aureus. Antimicrob Agents Chemother44:.https://doi.org/10.1128/aac.44.6.1549-1555.2000

2) Ito T, Katayama Y, Asada K, Mori N, Tsutsumimoto K, Tiensasitorn C, Hiramatsu K.2001.Structural Comparison of Three Types of Staphylococcal Cassette Chromosome mec Integrated in the Chromosome in Methicillin-Resistant Staphylococcus aureus. Antimicrob Agents Chemother45:.https://doi.org/10.1128/aac.45.5.1323-1336.2001

3) Ma XX, Ito T, Tiensasitorn C, Jamklang M, Chongtrakool P, Boyle-Vavra S, Daum RS, Hiramatsu K2002.Novel Type of Staphylococcal Cassette Chromosome mec Identified in Community-Acquired Methicillin-Resistant Staphylococcus aureus Strains. Antimicrob Agents Chemother46:.https://doi.org/10.1128/aac.46.4.1147-1152.2002

4) Ito TMa XX, Takeuchi F, Okuma K, Yuzawa H, Hiramatsu K.2004.Novel Type V Staphylococcal Cassette Chromosome mec Driven by a Novel Cassette Chromosome Recombinase, ccrC. Antimicrob Agents Chemother48:.https://doi.org/10.1128/aac.48.7.2637-2651.2004

5) Oliveira DCMilheiriço C, de Lencastre H2006.Redefining a Structural Variant of Staphylococcal Cassette Chromosome mec, SCCmec Type VI. Antimicrob Agents Chemother50:.https://doi.org/10.1128/aac.00629-06

6) Berglund CIto TIkeda M, Ma XX, Söderquist B, Hiramatsu K2008.Novel Type of Staphylococcal Cassette Chromosome mec in a Methicillin-Resistant Staphylococcus aureus Strain Isolated in Sweden. Antimicrob Agents Chemother52:.https://doi.org/10.1128/aac.00087-08

7) Zhang KMcClure J, Elsayed S, Conly JM2009.Novel Staphylococcal Cassette Chromosome mec Type, Tentatively Designated Type VIII, Harboring Class A mec and Type 4 ccr Gene Complexes in a Canadian Epidemic Strain of Methicillin-Resistant Staphylococcus aureus. Antimicrob Agents Chemother53:.https://doi.org/10.1128/aac.01118-08

8) Li SSkov RL, Han X, Larsen AR, Larsen J, Sørum M, Wulf M, Voss A, Hiramatsu K, Ito T2011.Novel Types of Staphylococcal Cassette Chromosome mec Elements Identified in Clonal Complex 398 Methicillin-Resistant Staphylococcus aureus Strains . Antimicrob Agents Chemother55:.https://doi.org/10.1128/aac.01475-10

9) García-Álvarez L, Holden MT, Lindsay H, et al. Meticillin-resistant Staphylococcus aureus with a novel mecA homologue in human and bovine populations in the UK and Denmark: a descriptive study. Lancet Infect Dis. 2011;11(8):595-603. doi:10.1016/S1473-3099(11)70126-8

10) Wu Z, Li F, Liu D, Xue H, Zhao X. 2015. Novel type XII staphylococcal
cassette chromosome mec harboring a new cassette chromosome recombinase,
CcrC2. Antimicrob Agents Chemother 59:7597–7601. doi:10.1128/AAC.01692-15.

11) Sharmin Baig, Thor Bech Johannesen, Søren Overballe-Petersen, Jesper Larsen, Anders Rhod Larsen, Marc Stegger, Novel SCCmec type XIII (9A) identified in an ST152 methicillin-resistant Staphylococcus aureus,
Infection, Genetics and Evolution, Volume 61, 2018, Pages 74-76,
ISSN 1567-1348, https://doi.org/10.1016/j.meegid.2018.03.013.

12) Noriko Urushibara, Meiji Soe Aung, Mitsuyo Kawaguchiya, Nobumichi Kobayashi, Novel staphylococcal cassette chromosome mec (SCCmec) type XIV (5A) and a truncated SCCmec element in SCC composite islands carrying speG in ST5 MRSA in Japan, Journal of Antimicrobial Chemotherapy, Volume 75, Issue 1, January 2020, Pages 46–50, https://doi.org/10.1093/jac/dkz406

13) Wei Wang, Yue Hu, Michelle Baker, Tania Dottorini, Hui Li, Yinping Dong, Yao Bai, Séamus Fanning, Fengqin Li, Novel SCCmec type XV (7A) and two pseudo-SCCmec variants in foodborne MRSA in China, Journal of Antimicrobial Chemotherapy, Volume 77, Issue 4, April 2022, Pages 903–909, https://doi.org/10.1093/jac/dkab500