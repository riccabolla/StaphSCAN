import pandas as pd
import subprocess
import io
import tempfile
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align


@dataclass
class GeneHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    slen: int
    qlen: int
    sstart: int
    send: int
    bitscore: float

    @property
    def coverage(self) -> float:
        return (self.length / self.qlen) * 100 if self.qlen else 0.0

    @property
    def family(self) -> str:
        return self.qseqid.split('_')[0] if self.qseqid else ''

    @property
    def aro_id(self) -> str:
        parts = self.qseqid.split('_')
        if len(parts) > 1 and parts[1].isdigit():
            return parts[1]
        return 'Unknown'

    @property
    def clean_name(self) -> str:
        parts = self.qseqid.split('_')
        if len(parts) > 2:
            return '_'.join(parts[:-1])
        return self.qseqid

    def annotation_tag(self) -> str:
        return f"(Id:{int(self.pident)}% Cov:{int(self.coverage)}%)"


class Module:
    def __init__(self):
        self.name = "resistance"
        self.module_dir = Path(__file__).resolve().parent
        self.db_fasta = self.module_dir / "data" / "targets.fasta"
        self.aro_csv = self.module_dir / "data" / "aro_index.csv"
        self.aro_tsv = self.module_dir / "data" / "aro_categories_index.tsv"

        self.card_data: Dict[str, Dict[str, str]] = {}
        self.ref_prot_dict: Dict[str, str] = {}

        self.mutation_targets = ["gyrA", "parC", "grlA", "gyrB", "rpoB"]
        self.known_mutations = {
            "gyrA": {84: ('S', ['L']), 88: ('S', ['P'])},
            "parC": {80: ('S', ['F', 'Y']), 84: ('E', ['K', 'G', 'V'])},
            "grlA": {80: ('S', ['F', 'Y'])},
            "rpoB": {481: ('H', ['N', 'Y']), 527: ('I', ['M'])}
        }

        self.mutation_phenotypes = {
            "gyrA": "Fluoroquinolones-R",
            "parC": "Fluoroquinolones-R",
            "grlA": "Fluoroquinolones-R",
            "rpoB": "Rifampicin-R",
        }

        self.category_map = {
            "mecA": ("Mec_RES", "MRSA"),
            "mecC": ("Mec_RES", "MRSA"),
            "blaZ": ("Beta_lactamases", "Penicillin-R"),
            "vanA": ("Other_RES", "VRSA"),
        }

        self.aligner_dna = Align.PairwiseAligner()
        self.aligner_dna.mode = 'global'
        self.aligner_dna.match_score = 1
        self.aligner_dna.mismatch_score = -1
        self.aligner_dna.open_gap_score = -2
        self.aligner_dna.extend_gap_score = -0.5

        self.aligner_prot = Align.PairwiseAligner()
        self.aligner_prot.mode = 'global'
        self.aligner_prot.match_score = 5
        self.aligner_prot.mismatch_score = -4

        self.aligner_mut = Align.PairwiseAligner()
        self.aligner_mut.mode = 'global'
        self.aligner_mut.match_score = 5
        self.aligner_mut.mismatch_score = -4
        self.aligner_mut.open_gap_score = -10
        self.aligner_mut.extend_gap_score = -0.5

        self.load_card_data()
        self.load_ref_seqs()

    def load_card_data(self) -> None:
        file_to_read = None
        sep = ','
        if self.aro_tsv.exists():
            file_to_read = self.aro_tsv
            sep = '\t'
        elif self.aro_csv.exists():
            file_to_read = self.aro_csv
            sep = ','

        if not file_to_read:
            return

        try:
            df_meta = pd.read_csv(file_to_read, sep=sep)
            for _, row in df_meta.iterrows():
                acc = str(row.get('ARO Accession', '')).replace('ARO:', '')
                self.card_data[acc] = {
                    'drug_class': row.get('Drug Class', '-'),
                    'mechanism': row.get('Resistance Mechanism', '-')
                }
        except Exception:
            return

    def load_ref_seqs(self) -> None:
        if not self.db_fasta.exists():
            return

        try:
            for record in SeqIO.parse(self.db_fasta, 'fasta'):
                prot_raw = str(record.seq.translate(table=11))
                start_idx = prot_raw.find('M')
                if start_idx != -1:
                    self.ref_prot_dict[record.id] = prot_raw[start_idx:].strip('*')
                else:
                    self.ref_prot_dict[record.id] = prot_raw.strip('*')
        except Exception:
            return

    def check_db(self) -> bool:
        if not self.db_fasta.exists():
            print(f"Error: Database fasta not found at {self.db_fasta}")
            return False
        return True

    def extract_gene_from_assembly(self, assembly_seqs: Dict[str, Seq], contig_id: str, start: int, end: int) -> Seq:
        if contig_id not in assembly_seqs:
            return Seq("")
        full_seq = assembly_seqs[contig_id].seq
        if start < end:
            return full_seq[start-1:end]
        return full_seq[end-1:start].reverse_complement()

    def get_best_frame_translation(self, dna_seq: Seq, ref_prot: str) -> str:
        best_score = -float('inf')
        best_prot = ''

        for frame in range(3):
            sliced = dna_seq[frame:]
            trim = len(sliced) % 3
            if trim > 0:
                sliced = sliced[:-trim]
            if not sliced:
                continue

            cand_prot = str(sliced.translate(table=11)).strip('*')
            try:
                score = self.aligner_dna.score(ref_prot, cand_prot)
            except Exception:
                score = 0

            if score > best_score:
                best_score = score
                best_prot = cand_prot

        return best_prot

    def trim_to_reference_alignment(self, prot_found: str, prot_ref: str) -> str:
        if not prot_ref or not prot_found:
            return prot_found

        aligner = self.aligner_prot
        try:
            aligner.end_insertion_score = 0.0
            aligner.end_deletion_score = 0.0
        except AttributeError:
            try:
                aligner.target_end_gap_score = 0.0
                aligner.query_end_gap_score = 0.0
            except Exception:
                pass

        try:
            alignment = aligner.align(prot_ref, prot_found)[0]
            aln_ref = alignment[0]
            aln_found = alignment[1]

            found_idx = 0
            found_start_trim = 0
            match_started = False

            for r_char, f_char in zip(aln_ref, aln_found):
                if r_char == '-':
                    if f_char != '-':
                        found_idx += 1
                    continue

                found_start_trim = found_idx
                match_started = True
                break

            if match_started:
                return prot_found[found_start_trim:]
        except Exception:
            pass

        return prot_found

    def check_mutations(self, fam: str, prot_found: str, prot_ref: str) -> List[str]:
        if fam not in self.known_mutations:
            return []

        try:
            alignment = self.aligner_mut.align(prot_ref, prot_found)[0]
        except Exception:
            return []

        aligned_ref = alignment[0]
        aligned_found = alignment[1]

        found_muts: List[str] = []
        target_loci = self.known_mutations[fam]

        ref_pos_counter = 0
        for ref_char, found_char in zip(aligned_ref, aligned_found):
            if ref_char != '-':
                ref_pos_counter += 1

            if ref_pos_counter in target_loci:
                expected_ref, resistant_list = target_loci[ref_pos_counter]
                if found_char != '-' and found_char != expected_ref:
                    if found_char in resistant_list:
                        found_muts.append(f"{expected_ref}{ref_pos_counter}{found_char}")

        return found_muts

    def run_blast(self, assembly_path: Path) -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_db = Path(temp_dir) / 'assembly_db'
            cmd_makedb = [
                'makeblastdb', '-in', str(assembly_path), '-dbtype', 'nucl', '-out', str(temp_db)
            ]
            try:
                subprocess.run(cmd_makedb, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                return pd.DataFrame()

            cmd_blast = [
                'blastn', '-query', str(self.db_fasta), '-db', str(temp_db),
                '-outfmt', '6 qseqid sseqid pident length slen qlen sstart send bitscore',
                '-max_target_seqs', '1000'
            ]

            res = subprocess.run(cmd_blast, capture_output=True, text=True)
            if not res.stdout:
                return pd.DataFrame()

            df = pd.read_csv(io.StringIO(res.stdout), sep='\t',
                             names=["qseqid", "sseqid", "pident", "length", "slen", "qlen", "sstart", "send", "bitscore"])
            return df

    def prepare_blast_hits(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        df = df[df['pident'] >= 40.0].copy()
        if df.empty:
            return df
        df['coverage'] = (df['length'] / df['qlen']) * 100
        df['family'] = df['qseqid'].apply(lambda x: x.split('_')[0])
        df = df.sort_values('bitscore', ascending=False)
        return df

    def make_default_output(self) -> Dict[str, str]:
        headers = [
            "Input_sequence_ID", "Input_gene_start", "Input_gene_stop", "Input_gene_length",
            "Reference_accession", "Sequence_identity", "Coverage",
            "Drug_class", "Resistance_mechanism",
            "truncated_resistance_hits", "spurious_resistance_hits",
            "res_genes", "res_phenotype", "res_score", "res_mutations",
            "Mec_RES", "Beta_lactamases", "Fluoroquinolones", "Other_RES",
            "Mec_AA_Found", "Mec_AA_Ref"
        ]
        defaults = {k: '-' for k in headers}
        defaults['res_score'] = '0'
        return defaults

    def build_output_from_containers(self, defaults: Dict[str, str], best_df: pd.DataFrame,
                                     strong_genes: List[str], mutations_found: List[str],
                                     trunc_list: List[str], spur_list: List[str],
                                     cat_mec: List[str], cat_bla: List[str], cat_fq: List[str], cat_other: List[str],
                                     mec_aa_found: List[str], mec_aa_ref: List[str], phenotypes: List[str]) -> Dict[str, str]:
        out = defaults.copy()
        if not best_df.empty:
            top = best_df.iloc[0]
            out['Input_sequence_ID'] = str(top['sseqid'])
            out['Input_gene_start'] = str(top['sstart'])
            out['Input_gene_stop'] = str(top['send'])
            out['Input_gene_length'] = str(top['length'])

        mapping = {
            'res_genes': strong_genes,
            'res_mutations': mutations_found,
            'truncated_resistance_hits': trunc_list,
            'spurious_resistance_hits': spur_list,
            'Mec_RES': cat_mec,
            'Beta_lactamases': cat_bla,
            'Fluoroquinolones': cat_fq,
            'Other_RES': cat_other,
            'Mec_AA_Found': mec_aa_found,
            'Mec_AA_Ref': mec_aa_ref,
        }

        for key, val in mapping.items():
            out[key] = '; '.join(val) if val else '-'

        uniq_pheno = sorted(list(set(phenotypes)))
        out['res_phenotype'] = ' + '.join(uniq_pheno) if uniq_pheno else 'Susceptible'

        res_score = 0
        if 'VRSA' in phenotypes:
            res_score = 3
        elif 'MRSA' in phenotypes:
            res_score = 2
        elif 'Penicillin-R' in phenotypes:
            res_score = 1
        out['res_score'] = str(res_score)

        return out

    def run(self, assembly_path: str) -> Dict[str, str]:
        defaults = self.make_default_output()

        if not self.check_db():
            return defaults

        try:
            assembly_seqs = SeqIO.to_dict(SeqIO.parse(assembly_path, 'fasta'))
        except Exception:
            return defaults

        df = self.run_blast(Path(assembly_path))
        if df.empty:
            return defaults

        df = self.prepare_blast_hits(df)
        if df.empty:
            return defaults

        best_hits = df.drop_duplicates(subset=['family'])

        # containers
        strong_genes: List[str] = []
        mutations_found: List[str] = []
        trunc_list: List[str] = []
        spur_list: List[str] = []
        phenotypes: List[str] = []

        cat_mec: List[str] = []
        cat_bla: List[str] = []
        cat_fq: List[str] = []
        cat_other: List[str] = []

        mec_aa_found: List[str] = []
        mec_aa_ref: List[str] = []

        for _, row in best_hits.iterrows():
            hit = GeneHit(
                qseqid=row['qseqid'], sseqid=row['sseqid'], pident=row['pident'], length=row['length'],
                slen=row['slen'], qlen=row['qlen'], sstart=int(row['sstart']), send=int(row['send']), bitscore=row['bitscore']
            )

            fam = hit.family
            aro = hit.aro_id
            meta = self.card_data.get(aro, {})
            d_class = meta.get('drug_class', '-')

            dna_found = self.extract_gene_from_assembly(assembly_seqs, hit.sseqid, hit.sstart, hit.send)
            prot_ref = self.ref_prot_dict.get(hit.qseqid, '')
            prot_found_raw = self.get_best_frame_translation(dna_found, prot_ref)
            prot_found = self.trim_to_reference_alignment(prot_found_raw, prot_ref)

            # Mutation targets handling
            if fam in self.mutation_targets and prot_ref:
                fam_muts = self.check_mutations(fam, prot_found, prot_ref)
                if fam_muts:
                    mut_str = ','.join(fam_muts)
                    mutations_found.append(f"{fam} [{mut_str}]")
                    pheno = self.mutation_phenotypes.get(fam)
                    if pheno:
                        phenotypes.append(pheno)
                        if fam in ["gyrA", "parC", "grlA"]:
                            cat_fq.append(f"{fam} [{mut_str}]")
                        else:
                            cat_other.append(f"{fam} [{mut_str}]")
                continue

            display_str = f"{hit.clean_name}"

            if hit.pident < 100.0:
                if prot_ref and prot_found == prot_ref:
                    pass
                else:
                    display_str += "*"

            info_tag = hit.annotation_tag()

            is_spurious = (hit.pident < 90.0)
            is_truncated = (hit.coverage < 90.0)
            is_pseudogene = ('*' in prot_found)

            # pseudogenes for selected families
            if is_pseudogene and fam in ["mecA", "mecC", "blaZ"] and hit.coverage > 90.0:
                is_pseudogene = False
                display_str += "(fs)"

            if is_spurious:
                spur_list.append(f"{display_str} {info_tag} [Low Identity]")
                continue
            if is_pseudogene:
                trunc_list.append(f"{display_str} {info_tag} [Pseudogene]")
                continue
            if is_truncated:
                trunc_list.append(f"{display_str} {info_tag} [Truncated]")
                continue

            strong_genes.append(display_str)
            polished_str = display_str.replace('*', '')

            # category assignment
            if fam in ("mecA", "mecC"):
                cat_mec.append(fam)
                phenotypes.append('MRSA')
                mec_aa_found.append(f"{fam}_Found:{prot_found}")
                mec_aa_ref.append(f"{fam}_Ref:{prot_ref}")
            elif fam == 'blaZ':
                cat_bla.append(polished_str)
                phenotypes.append('Penicillin-R')
            elif fam == 'vanA':
                cat_other.append(polished_str)
                phenotypes.append('VRSA')
            else:
                cat_other.append(polished_str)
                if d_class and d_class != '-':
                    phenotypes.append(d_class)

        final_out = self.build_output_from_containers(
            defaults, df, strong_genes, mutations_found,
            trunc_list, spur_list, cat_mec, cat_bla, cat_fq, cat_other,
            mec_aa_found, mec_aa_ref, phenotypes
        )

        return final_out
