## Biofilm

``-m biofilm``

The biofilm module identifies key *Staphylococcus aureus* genes involved in biofilm formation, including surface proteins (*clfA*, *clfB*, *fnbA*, *fnbB*), polysaccharide intercellular adhesin genes (*icaA–D*), and the *icaR* regulator. These genes contribute to adhesion, biofilm maturation, and regulation of the biofilm phenotype (1).

A curated FASTA file containing representative biofilm gene sequences is bundled with the module and stored in the ``/data`` directory.

All the genes found are annotated as follow:

* No annotation: an exact nucelotide match is found

* ``^``: inexact nucleotide match but perfect amino acidic match

* ``*``: inexact nucleotide and inexact amino acidic match

* ``?``: incomplete match

* ``-X%``: truncated amino acidic sequence

### Parameters

Hits are filtered by minimum alignment identity and coverage:

``--min_id_biofilm`` : Minimum alignment percentage identity (default: 90)

``--min_cov_biofilm``: Minimum alignment percentage coverage (default: 80)

### Biofilm score

A biofilm score is calculated, based on the genes detected (2) (3). 

It follows this criteria:

* 3: fnbAB and icaADBC complete

* 2: clfAB and icaADBC complete

* 1: A single group complete

* 0: no complete groups detected

### Output

Results are grouped for each gene-family, and reported individually if any imperfect match is found.

| Field                          | Description                                                   |
| ------------------------------ | ------------------------------------------------------------- |
| biofilm_score                  | Biofilm score (0–3)                                           |
| biofilm_truncated_hits         | Genes with truncations or premature stop codons               |
| clfAB                          | Status of clfA and clfB genes (Complete/Incomplete/-)         |
| fnbAB                          | Status of fnbA and fnbB genes (Complete/Incomplete/-)         |
| icaADBC                        | Status of icaA–D genes (Complete/Incomplete/-)                |
| clf_genes                      | Detected clf genes with annotation (truncated, partial, etc.) |
| fnb_genes                      | Detected fnb genes with annotation                            |
| ica_genes                      | Detected ica genes with annotation                            |
| icaR_mutations (4)             | Truncated icaR regulator, if any                              |


### Citations

1) Idrees M, Sawant S, Karodia N, Rahman A. Staphylococcus aureus Biofilm: Morphology, Genetics, Pathogenesis and Treatment Strategies. Int J Environ Res Public Health. 2021 Jul 16;18(14):7602. doi: 10.3390/ijerph18147602. 

2) O'Neill E, Pozzi C, Houston P, et al. A novel Staphylococcus aureus biofilm phenotype mediated by the fibronectin-binding proteins, FnBPA and FnBPB. J Bacteriol. 2008;190(11):3835-3850. doi:10.1128/JB.00167-08

3) Foster TJ. Surface Proteins of Staphylococcus aureus. Microbiol Spectr. 2019;7(4):10.1128/microbiolspec.gpp3-0046-2018. doi:10.1128/microbiolspec.GPP3-0046-2018

4) Schwartbeck B, Rumpf CH, Hait RJ, Janssen T, Deiwick S, Schwierzeck V, Mellmann A, Kahl BC. Various mutations in icaR, the repressor of the icaADBC locus, occur in mucoid Staphylococcus aureus isolates recovered from the airways of people with cystic fibrosis. Microbes Infect. 2024 May-Jun;26(4):105306. doi: 10.1016/j.micinf.2024.105306.