## Biofilm

``-m biofilm``

The biofilm module identifies key *Staphylococcus aureus* genes involved in biofilm formation, including surface proteins (*clfA*, *clfB*, *fnbA*, *fnbB*), polysaccharide intercellular adhesin genes (*icaA–D*), and the *icaR* regulator. These genes contribute to adhesion, biofilm maturation, and regulation of the biofilm phenotype.

A curated FASTA file containing representative biofilm gene sequences is bundled with the module and stored in the ``/data`` directory.

All the genes found are annotated as follow:
* No annotation: an exact nucelotide match is found

* ``^``: inexact nucleotide match but perfect amino acidic match

* ``*``: inexact nucleotide and inexact amino acid match

* ``?``: incomplete match

* ``-X%``: truncated amino acid sequence

### Parameters

Hits are filtered by minimum alignment identity and coverage:

``--min_id_biofilm`` : Minimum alignment percentage identity (default: 90)

``--min_cov_biofilm``: Minimum alignment percentage coverage (default: 80)

### Biofilm score

A biofilm score is calculated, based on the genes detected. 

It follows this criteria:
* 3 – fnbAB and icaADBC complete
* 2 – clfAB and icaADBC complete
* 1 – A single group complete
* 0 – no complete groups detected

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
| icaR_mutations                 | Truncated icaR regulator, if any                              |
