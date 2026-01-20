## Virulence

``-m virulence``

The virulence module screens *Staphylococcus aureus* genomes for major toxin genes associated with specific clinical syndromes, ranging from food poisoning to necrotizing pneumonia and toxic shock syndrome.

It detects virulence determinants across the following categories:

- **Panton-Valentine Leukocidin (PVL)** (necrotizing pneumonia, recurrent abscesses)
- **Toxic Shock Syndrome Toxin-1** (TSST-1)
- **Exfoliative Toxins** (Scalded Skin Syndrome)
- **Leukocidins** (Skin and soft tissues infection)
- **Staphylococcal Enterotoxins** (Food poisoning)

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

- **4**: PVL detected (*lukS-PV* + *lukF-PV*)
- **3**: TSST-1 detected (*tst* or *tsst1*)
- **2**: Exfoliative toxins detected (*eta*, *etb*)
- **1**: Enterotoxins (*sea*—*set*) OR LukED detected
- **0**: No major accessory virulence factors detected

The following criteria also apply:

* For lukFS and lukED both subunits must be found to be contribute to the score

*  Genes reported in the `truncated_virulence_hits` and `spurious_virulence_hits` columns do not contribute to the virulence score.

### Output

Results are reported as follows:

| Field | Description |
|------|------------|
| `vir_score` | Virulence score (0–4) |
| `vir_pvl` | Detection status of Panton-Valentine Leukocidin (Positive, Partial, or -) |
| `vir_tsst` | Detected Toxic Shock Syndrome Toxin genes |
| `vir_et` | Detected exfoliative toxin genes (`eta`, `etb`) |
| `vir_lukED` | Detection status of Leukocidin ED (Positive, Partial, or -) |
| `vir_se` | Detected enterotoxin genes (e.g., `sea`, `sec`, `seh`) |
| `truncated_virulence_hits` | Genes with truncations or premature stop codons |
| `spurious_virulence_hits` | Genes with weak hits (low identity/coverage) |