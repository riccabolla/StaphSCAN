## Resistance 

``-m resistance``

The resistance module detects acquired antimicrobial resistance genes and chromosomal point mutations associated with reduced antimicrobial susceptibility in *Staphylococcus aureus*.

It screens for resistance determinants across multiple antibiotic classes, including:

- Aminoglycosides  
- Beta-lactams  
- Fluoroquinolones  
- Linezolid
- MLSB  
- Rifampicin  
- Tetracyclines  
- Vancomycin  

A curated FASTA file containing representative resistance gene sequences is bundled with the module and stored in the `/data` directory.

All detected genes are annotated as follows:

- No annotation**: exact nucleotide match  
- ``^``: inexact nucleotide match but perfect amino acid match (synonymous mutation)  
- ``*``: inexact nucleotide and inexact amino acid match  
- ``?``: incomplete match (low coverage)  
- ``-X%``: truncated amino acid sequence  

The module screens specific chromosomal loci for point mutations known to confer resistance.

| Locus | Mutation targets | Associated phenotype |
|------|------------------|----------------------|
| gyrA | S84L, S88P | Fluoroquinolones-R |
| gyrB | T451S | Fluoroquinolones-R |
| parC | S80F, S80Y, E84K, E84G, E84V | Fluoroquinolones-R |
| rpoB | H481Y | Rifampicin-R |
| 23S rRNA | G2576T, G2447T, T2500A | Linezolid-R |

### Parameters

Hits are filtered using minimum alignment identity and coverage thresholds:

- ``--min_id_res``  
  Minimum alignment percentage identity (default: 90)

- ``--min_cov_res``  
  Minimum alignment percentage coverage (default: 80)

All hits with 80–90% identity or 40–80% coverage are reported as spurious hits.

### Resistance Score

A resistance score is calculated to categorize isolates based on a hierarchy of clinical resistance.

Scoring criteria:

- **3**: `vanA` detected (Vancomycin resistance)  
- **2**: `mecA` or `mecC` detected (Methicillin resistance – MRSA)  
- **1**: `blaZ` detected (Penicillin/Beta-lactamase resistance)  
- **0**: No key resistance drivers detected  

### Output

Results are grouped by antibiotic class and reported individually.

| Field | Description |
|------|------------|
| `res_score` | Resistance score (0–3) |
| `res_gene_count` | Number of genes conferring resistance |
| `res_class_count` | Number of classes with at least one resistance determinants |
| `Aminoglycosides` | Detected aminoglycoside resistance genes (e.g., AAC(6')-Ie-APH(2'')-Ia) |
| `Mec_RES` | Detected methicillin resistance genes (`mecA`, `mecC`) |
| `Beta_lactamases` | Detected beta-lactamase genes (`blaZ`) |
| `Fluoroquinolones` | Genes/mutations conferring fluoroquinolone resistance |
| `MLSB` | Detected macrolide/streptogramin B/lincosamide resistance genes |
| `Rifampicin` | Genes/mutations conferring rifampicin resistance |
| `Tetracyclines` | Detected tetracycline resistance genes (tet family) |
| `Linezolid` | Genes/mutations conferring linezolid resistance (`cfrA`, 23S rRNA) |
| `Vancomycin` | Detected vancomycin resistance genes (`vanA`) |
| `truncated_resistance_hits` | Genes with truncations |
| `spurious_resistance_hits` | Genes with weak hits |

For `res_gene_count` and `res_class_count` the following criteria is used:

* The presence of resistance mutations do not contribute to the resistance gene count

* Mutations do contribute to the drug class count

* Genes reported in the truncated_resistance_hits and spurious_resistance_hits columns do not contribute to the counts