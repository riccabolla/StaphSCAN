## Virulence

``-m virulence``

The virulence module screens *Staphylococcus aureus* genomes for major toxin genes associated with specific clinical syndromes.

The Panton-Valentine Leukocidin (PVL), encoded by the subunits *lukF-PVL* and *lukS-PVL*, is an exotoxin responsible of skin and soft-tissue infection, and in some cases of life-threatening infections (ex- necrotizing pneumonia) ([Shallcross et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3530297/); [Hoppe et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6756729/); [ Labandeira-Rey et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17234914/)). 
Other two important leukocidins, are *lukD* and *lukE*, which enhance virulence and contribute to *Staphylococcus aureus* pathogenesis ([Alonzo et al. 2012](https://pubmed.ncbi.nlm.nih.gov/22142035/)).

The Toxic Shock Syndrome Toxin-1 (TSST-1) is an enterotoxin associated with the development of the toxic shock syndrome, a severe life-threatening condition ([Manosur et al. 2021](https://www.mdpi.com/2072-6651/13/10/677)) 

The exofoliative toxin *etA*, *etB*, *etD*, and *etE* can cause the staphylococcal scaled skin syndrome (SSSS), a severe infection ([Nishifuji et al. 2008](https://pubmed.ncbi.nlm.nih.gov/17582744/)).

The *Staphylococcal* enterotoxins (SE) are mainly responsible of food poisoning, and are classified as *superantigens*, due to the ability to stimulate large populations of T cells ([Pinchuk et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC3153290/); [Choi et al. 1989](https://pubmed.ncbi.nlm.nih.gov/2479030/)).   


A curated FASTA file containing representative virulence gene sequences is bundled with the module and stored in the `/data` directory.

All detected genes are annotated as follows:

- **No annotation**: exact nucleotide match
- ``^``: inexact nucleotide match but perfect amino acid match
- ``*``: inexact nucleotide and inexact amino acid match
- ``?``: incomplete match
- ``-X%``: truncated amino acid sequence

### Parameters

Hits are filtered using minimum alignment identity and coverage thresholds:

- ``--min_id_vir``  
  Minimum alignment percentage identity (default: 90)

- ``--min_cov_vir``  
  Minimum alignment percentage coverage (default: 80)

All hits with 80–90% identity or 40–80% coverage are reported as spurious hits.

### Virulence Score

A virulence score is calculated to categorize isolates based on a hierarchy of clinical severity and epidemic potential. The score reflects the highest-risk toxin detected.

Scoring criteria:

- **3**: TSST-1 detected and/or PVL detected
- **2**: Exfoliative toxins detected
- **1**: *Staphylococcal* enterotoxins (SE) OR LukED detected
- **0**: No major accessory virulence factors detected

The following criteria also apply:

* For lukFS and lukED both subunits must be found to contribute to the score

* Genes reported in the `truncated_virulence_hits` and `spurious_virulence_hits` columns do not contribute to the virulence score.

* Only the "most relevant" gene is considered, the cumulative presence of multiple genes do not contribute to the score. 

### Output

Results are reported as follows:

| Field | Description |
|------|------------|
| `vir_score` | Virulence score (0–3) |
| `vir_pvl` | Detection of Panton-Valentine Leukocidin (Positive, Partial, or -) |
| `vir_tsst` | Detected Toxic Shock Syndrome Toxin genes |
| `vir_et` | Detected exfoliative toxin genes |
| `vir_lukED` | Detection of Leukocidin ED (Positive, Partial, or -) |
| `vir_se` | Detected *Staphylococcal* enterotoxin genes |
| `truncated_virulence_hits` | Genes with truncations or premature stop codons |
| `spurious_virulence_hits` | Genes with weak hits |

### List of genes

| Gene |Accession^ |Citation^^|
|------|----------|--------|
|etA|[AAA17490](https://www.ncbi.nlm.nih.gov/protein/AAA17490)|[Sakurai et al. 1988](https://pubmed.ncbi.nlm.nih.gov/3183619/)|
|etB|[ WP_010994026](https://www.ncbi.nlm.nih.gov/protein/WP_010994026)|[Lee at al. 1987](https://pubmed.ncbi.nlm.nih.gov/3040666/)|
|etD|[QNW86327.1](http://www.ncbi.nlm.nih.gov/protein/QNW86327.1)|[Caban et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7670629/)|
|etE|[AUG74510.1](https://www.ncbi.nlm.nih.gov/protein/AUG74510.1)|[Le Marechal et al. 2011](https://pubmed.ncbi.nlm.nih.gov/21398544/)|
|lukD|[WP_000782464](https://www.ncbi.nlm.nih.gov/protein/WP_000782464); [WP_000782469](https://www.ncbi.nlm.nih.gov/protein/WP_000782469); [YP_500453](https://www.ncbi.nlm.nih.gov/protein/YP_500453); [WP_000782463](https://www.ncbi.nlm.nih.gov/protein/WP_000782463)|[Gravet et al. 1998](https://pubmed.ncbi.nlm.nih.gov/9781679/)|
|lukE| [WP_000473596](https://www.ncbi.nlm.nih.gov/protein/WP_000473596); [WP_000473600](https://www.ncbi.nlm.nih.gov/protein/WP_000473600); [YP_500454](https://www.ncbi.nlm.nih.gov/protein/YP_500454); [WP_000473593](https://www.ncbi.nlm.nih.gov/protein/WP_000473593) |[Gravet et al. 1998](https://pubmed.ncbi.nlm.nih.gov/9781679/)|
|lukF|[WP_024937002](https://www.ncbi.nlm.nih.gov/protein/WP_024937002) |[Kaneko et al. 1997](https://pubmed.ncbi.nlm.nih.gov/9404084/); [Badarau et al 2017](https://pubmed.ncbi.nlm.nih.gov/28455832/)|
|lukS|[WP_000239545](https://www.ncbi.nlm.nih.gov/protein/WP_000239545) |[Kaneko et al. 1997](https://pubmed.ncbi.nlm.nih.gov/9404084/); [Badarau et al 2017](https://pubmed.ncbi.nlm.nih.gov/28455832/)|
|sea|[WP_000750406](https://www.ncbi.nlm.nih.gov/protein/WP_000750406); [WP_000750412](https://www.ncbi.nlm.nih.gov/protein/WP_000750412) |[Betley et al. 1988](https://pubmed.ncbi.nlm.nih.gov/3335483/)|
|sec|[WP_001043550](https://www.ncbi.nlm.nih.gov/protein/WP_001043550); [WP_000278088](https://www.ncbi.nlm.nih.gov/protein/WP_000278088); [WP_001043551](https://www.ncbi.nlm.nih.gov/protein/WP_001043551) |[Bohach et al. 1987](https://pubmed.ncbi.nlm.nih.gov/2823067/); [Chen et al. 2001](https://pubmed.ncbi.nlm.nih.gov/11764893/)|
|seh|[WP_000608674](https://www.ncbi.nlm.nih.gov/protein/WP_000608674)|[Nilsson et al. 1999](https://pubmed.ncbi.nlm.nih.gov/10586065/)|
|sell|[WP_000746599](https://www.ncbi.nlm.nih.gov/protein/WP_000746599); [WP_000746597](https://www.ncbi.nlm.nih.gov/protein/WP_000746597)|[Baba et al. 2008](https://pubmed.ncbi.nlm.nih.gov/17951380/) |
|selk|[WP_000733775](https://www.ncbi.nlm.nih.gov/protein/WP_000733775); [WP_000733771](https://www.ncbi.nlm.nih.gov/protein/WP_000733771); [WP_000734020](https://www.ncbi.nlm.nih.gov/protein/WP_000734020)|-|
|selq|[WP_001033316](https://www.ncbi.nlm.nih.gov/protein/WP_001033316);[WP_001033320](https://www.ncbi.nlm.nih.gov/protein/WP_001033320);[ WP_001033321](https://www.ncbi.nlm.nih.gov/protein/WP_001033321) |[McClure et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30042755/)|
|tsst-1|[WP_001035599](https://www.ncbi.nlm.nih.gov/protein/WP_001035599); [WP_001035596](https://www.ncbi.nlm.nih.gov/protein/WP_001035596) |[Blomster-Hautamaa et al. 1986](https://pubmed.ncbi.nlm.nih.gov/3782090/)|

*Note*: <br>
^ One accession per variant is incuded. <br>
^^ The citation is the one provided at the accession page