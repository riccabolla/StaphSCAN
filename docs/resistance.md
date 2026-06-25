## Resistance 

``-m resistance``

The resistance module detects acquired antimicrobial resistance genes and chromosomal point mutations associated with reduced antimicrobial susceptibility in *Staphylococcus aureus*.

It screens for resistance determinants across multiple antibiotic classes, including:

- Aminoglycosides  
- Beta-lactams  
- Fluoroquinolones  
- Oxazolinidones
- MLSB  
- Rifampicin  
- Tetracyclines  
- Glycopeptides  

A curated set of FASTA files containing representative resistance gene sequences are bundled with the module and stored in the `/data` directory.

All detected genes are annotated as follows:

- No annotation**: exact nucleotide match  
- ``^``: inexact nucleotide match but perfect amino acid match
- ``*``: inexact nucleotide and inexact amino acid match  
- ``?``: incomplete match  
- ``-X%``: truncated amino acid sequence  

The module screens specific chromosomal loci for point mutations known to confer resistance.

| Locus | Mutation targets | Associated phenotype |
|------|------------------|----------------------|
| `gyrA` | `S84L`, `E88K`, `E88G` | `Fluoroquinolones-R` |
| `gyrB` | `T451S` | `Fluoroquinolones-R` |
| `parE` | `D432N` | `Fluoroquinolones-R` |
| `parC` | `S80F`, `S80Y`, `E84K`, `E84G`, `E84V` | `Fluoroquinolones-R` |
| `rpoB` | `H481Y`, `H481N`, `L466N`, `A473T`, `A477T` | `Rifampicin-R` |
| `23S rRNA` | `G2576T`, `G2447T`, `T2500A`, `A2058G`, `A2059G` | `Linezolid-R` |

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

- **3**: `vanA` or `vanB` detected (Vancomycin resistance)  
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
| `Amino_res` | Detected aminoglycoside resistance genes |
| `Bla_res` | Detected beta-lactamase genes |
| `Flq_res` | Genes/mutations conferring fluoroquinolone resistance |
| `Gly_res` | Detected glycopeptide resistance genes |
| `Mec_res` | Detected methicillin resistance genes |
| `MLSB_res` | Detected macrolide/streptogramin B/lincosamide resistance genes |
| `Oxa_res` | Genes/mutations conferring oxazolinidones resistance |
| `Rif_res` | Genes/mutations conferring rifampicin resistance |
| `Tet_res` | Detected tetracycline resistance genes (tet family) |
| `spurious_resistance_hits` | Genes with weak hits |
| `truncated_resistance_hits` | Genes with truncations |

For `res_gene_count` and `res_class_count` the following criteria is used:

* The presence of resistance mutations do not contribute to the resistance gene count

* Mutations do contribute to the drug class count

* Genes reported in the truncated_resistance_hits and spurious_resistance_hits columns do not contribute to the counts

### List of genes

| Gene | Drug Class | Variants^ | Reference number^^ | Citation |
|------| -----------|----------|------------------|----------|
|AAC(6')-Ie-APH(2'')-Ia|aminoglycoside antibiotic|1|[3002597](https://card.mcmaster.ca/ontology/38997)|[Daigle DM, et al. 1999](https://pubmed.ncbi.nlm.nih.gov/10021417/)
|ANT(4')-Ia|aminoglycoside antibiotic|1|[3002623](https://card.mcmaster.ca/ontology/39023)|[Santanam P and Kayser FH. 1978](https://pubmed.ncbi.nlm.nih.gov/659332/)|
|ANT(4')-Ib|aminoglycoside antibiotic|1|[3003905](https://card.mcmaster.ca/ontology/40608)|[McDougal LK, et al. 2010](http://www.ncbi.nlm.nih.gov/pubmed/20585117)|
|APH(3')-IIIa|aminoglycoside antibiotic|1|[3002647](https://card.mcmaster.ca/ontology/39047)|[Trieu-Cuot P and Courvalin P. 1983](http://www.ncbi.nlm.nih.gov/pubmed/6313476)|
|ANT(6)-Ia|aminoglycoside antibiotic|1|[3002626](https://card.mcmaster.ca/ontology/39026)|[Gill SR, et al. 200](http://www.ncbi.nlm.nih.gov/pubmed/15774886)|
|ANT(9)-Ia|aminoglycoside antibiotic|1|[3002630](https://card.mcmaster.ca/ontology/39030)|[Murphy et al. 1985](http://www.ncbi.nlm.nih.gov/pubmed/3004956)|
|mecA|penicillin beta-lactam|291|[3000617](https://card.mcmaster.ca/ontology/36911)|[Ubukata K, et al. 1989](http://www.ncbi.nlm.nih.gov/pubmed/2708325)|
|mecC|penicillin beta-lactam|6|[3001209](https://card.mcmaster.ca/ontology/37590)|[Garcia-Alvarez L, et al. 2011](http://www.ncbi.nlm.nih.gov/pubmed/21641281)|
|blaZ|penicillin beta-lactam|473|[3000621](https://card.mcmaster.ca/ontology/36963)|[McLaughlin JR, et al. 1981](http://www.ncbi.nlm.nih.gov/pubmed/6793593)|
|blaZ mecC-type|penicillin beta-lactam|5|[3005097](https://card.mcmaster.ca/ontology/43312)|[Shore AC, et al. 2011](http://www.ncbi.nlm.nih.gov/pubmed/21636525)|
|cfrA|phenicol antibiotic|1|[3003441](https://card.mcmaster.ca/ontology/40028)|[Schwarz S, et al. 2000](http://www.ncbi.nlm.nih.gov/pubmed/10952608)|
|poxtA|oxazolinidone antibiotic|1|[3004470](https://card.mcmaster.ca/ontology/41688)|[Antonelli A, et al. 2018](http://www.ncbi.nlm.nih.gov/pubmed/29635422)|
|ermA|streptogramin antibiotic, streptogramin B antibiotic, streptogramin A antibiotic, lincosamide antibiotic, macrolide antibiotic|1|[3000347](https://card.mcmaster.ca/ontology/36486)|[Malhotra-Kumar S, et al. 2008](http://www.ncbi.nlm.nih.gov/pubmed/18952616)|
|ermB|streptogramin antibiotic, streptogramin B antibiotic, streptogramin A antibiotic, lincosamide antibiotic, macrolide antibiotic|1|[3000375](https://card.mcmaster.ca/ontology/36514)|[Yu L, et al. 1997](http://www.ncbi.nlm.nih.gov/pubmed/9187657)|
|ermC|streptogramin antibiotic, streptogramin B antibiotic, streptogramin A antibiotic, lincosamide antibiotic, macrolide antibiotic|1|[3000250](https://card.mcmaster.ca/ontology/36389)|[Shivakumar AG and Dubnau D. 1981](http://www.ncbi.nlm.nih.gov/pubmed/6792593)|
|ermT|streptogramin antibiotic, streptogramin B antibiotic, streptogramin A antibiotic, lincosamide antibiotic, macrolide antibiotic|1|[3000595](https://card.mcmaster.ca/ontology/36734)|[Tannock GW, et al. 1994](http://www.ncbi.nlm.nih.gov/pubmed/8171126)|
|msrA|streptogramin antibiotic, streptogramin B antibiotic, macrolide antibiotic|1|[3000251](https://card.mcmaster.ca/ontology/36390)|[Poole K. 2005](http://www.ncbi.nlm.nih.gov/pubmed/15914491)|
|tetM|tetracycline antibiotic|3|[3000186](https://card.mcmaster.ca/ontology/36325); [SAV0398](https://www.genome.jp/dbget-bin/www_bget?sav:SAV0398); [M21136.1](https://www.ncbi.nlm.nih.gov/nuccore/M21136)|[Akhtar M, et al. 2009](http://www.ncbi.nlm.nih.gov/pubmed/19475445)|
|tetK|tetracycline antibiotic|1|[3000178](https://card.mcmaster.ca/ontology/36317)|[Roberts MC. 2005](http://www.ncbi.nlm.nih.gov/pubmed/15837373)|
|vanA|glycopeptide antibiotic|1|[3000010](https://card.mcmaster.ca/ontology/36019)|[Marshall CG, et al. 1997](http://www.ncbi.nlm.nih.gov/pubmed/9177243)|
|vanB|glycopeptide antibiotic|1|[3000013](https://card.mcmaster.ca/ontology/36022)|[Marshall CG, et al. 1997](http://www.ncbi.nlm.nih.gov/pubmed/9177243)|

*Note*:  

^ :  Number of variants in the database;

^^: Reference sequence of each variant. If multiple variants are present, but only one reference is reported, it means that all variants have the same reference from the same database.