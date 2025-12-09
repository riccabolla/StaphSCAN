import pandas as pd
import subprocess
import io
import tempfile
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align

@dataclass
class BiofilmHit:
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
    def gene_name(self) -> str:
        return self.qseqid.split('_')[0]

class Module:
    def __init__(self):
        self.name = "biofilm"
        self.module_dir = Path(__file__).resolve().parent
        self.db_fasta = self.module_dir / "data" / "biofilm.fasta"
        
        self.clf_targets = ['clfA', 'clfB']
        self.fnb_targets = ['fnbA', 'fnbB']
        self.ica_targets = ['icaA', 'icaB', 'icaC', 'icaD']
        self.regulator = 'icaR'

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
        
        self.ref_prot_dict = {}
        self.load_ref_seqs()

    def check_db(self) -> bool:
        if not self.db_fasta.exists():
            print(f"Error: Database fasta not found at {self.db_fasta}")
            return False
        return True

    def load_ref_seqs(self) -> None:
        if not self.db_fasta.exists(): return
        try:
            for record in SeqIO.parse(self.db_fasta, 'fasta'):
                prot_raw = str(record.seq.translate(table=11))
                start_idx = prot_raw.find('M')
                if start_idx != -1:
                    self.ref_prot_dict[record.id] = prot_raw[start_idx:].strip('*')
                else:
                    self.ref_prot_dict[record.id] = prot_raw.strip('*')
        except Exception as e:
            print(f"Warning: Error loading biofilm reference FASTA: {e}")

    def extract_gene(self, assembly_seqs, contig_id, start, end):
        if contig_id not in assembly_seqs: return Seq("")
        full_seq = assembly_seqs[contig_id].seq
        if start < end:
            return full_seq[start-1:end]
        return full_seq[end-1:start].reverse_complement()

    def get_best_frame_translation(self, dna_seq, ref_prot):
        best_score = -float('inf')
        best_prot = ''
        for frame in range(3):
            sliced = dna_seq[frame:]
            trim = len(sliced) % 3
            if trim > 0: sliced = sliced[:-trim]
            if not sliced: continue
            cand = str(sliced.translate(table=11)).strip('*')
            try:
                score = self.aligner_dna.score(ref_prot, cand)
            except: score = 0
            if score > best_score:
                best_score = score
                best_prot = cand
        return best_prot

    def trim_to_ref(self, prot_found, prot_ref):
        if not prot_ref or not prot_found: return prot_found
        aligner = self.aligner_prot
        try:
            aligner.end_insertion_score = 0.0
            aligner.end_deletion_score = 0.0
        except AttributeError:
            try:
                aligner.target_end_gap_score = 0.0
                aligner.query_end_gap_score = 0.0
            except: pass
        
        try:
            aln = aligner.align(prot_ref, prot_found)[0]
            aln_ref, aln_found = aln[0], aln[1]
            found_idx = 0
            start_trim = 0
            started = False
            for r, f in zip(aln_ref, aln_found):
                if r == '-': 
                    if f != '-': found_idx += 1
                    continue
                start_trim = found_idx
                started = True
                break
            if started: return prot_found[start_trim:]
        except: pass
        return prot_found

    def get_group_status(self, found_genes: set, target_genes: list) -> str:
        present_count = sum(1 for g in target_genes if g in found_genes)
        if present_count == len(target_genes):
            return "Complete"
        elif present_count > 0:
            return "Incomplete"
        else:
            return "-"

    def run(self, assembly_path: Path) -> Dict[str, str]:
        defaults = {
            "biofilm_score": "0",
            "biofilm_truncated_hits": "-",
            "clfAB": "-", "fnbAB": "-", "icaADBC": "-",
            "clfA": "0", "clfB": "0",
            "fnbA": "0", "fnbB": "0",
            "icaA": "0", "icaB": "0", "icaC": "0", "icaD": "0",
            "biofilm_genes": "-",
            "icaR_mutations": "-"
        }

        if not self.check_db(): return defaults

        try:
            assembly_seqs = SeqIO.to_dict(SeqIO.parse(assembly_path, 'fasta'))
        except: return defaults

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_db = Path(temp_dir) / 'assembly_db'
            subprocess.run(['makeblastdb', '-in', str(assembly_path), '-dbtype', 'nucl', '-out', str(temp_db)], 
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            cmd = ['blastn', '-query', str(self.db_fasta), '-db', str(temp_db),
                   '-outfmt', '6 qseqid sseqid pident length slen qlen sstart send bitscore',
                   '-max_target_seqs', '1000']
            
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: return defaults
            
            df = pd.read_csv(io.StringIO(res.stdout), sep='\t', 
                             names=["qseqid", "sseqid", "pident", "length", "slen", "qlen", "sstart", "send", "bitscore"])

        df = df[df['pident'] >= 90.0].copy()
        if df.empty: return defaults
        df['coverage'] = (df['length'] / df['qlen']) * 100
        df = df[df['coverage'] >= 80.0].copy()
        
        df['gene'] = df['qseqid'].apply(lambda x: x.split('_')[0])
        df = df.sort_values('bitscore', ascending=False).drop_duplicates(subset=['gene'])

        found_genes_functional = set()
        
        clf_list, fnb_list, ica_list, truncated_list = [], [], [], []
        icar_truncation = None

        for _, row in df.iterrows():
            hit = BiofilmHit(row['qseqid'], row['sseqid'], row['pident'], row['length'], row['slen'], 
                             row['qlen'], int(row['sstart']), int(row['send']), row['bitscore'])
            gene = hit.gene_name
            dna = self.extract_gene(assembly_seqs, hit.sseqid, hit.sstart, hit.send)
            ref_prot = self.ref_prot_dict.get(hit.qseqid, '')
            prot_raw = self.get_best_frame_translation(dna, ref_prot)
            prot_clean = self.trim_to_ref(prot_raw, ref_prot)
            
            is_pseudogene = ('*' in prot_clean)
            
            if is_pseudogene:
                stop_pos = prot_clean.find('*')
                pct = int((stop_pos / len(ref_prot)) * 100) if ref_prot else 0
                trunc_str = f"{gene}-{pct}%"
                truncated_list.append(trunc_str)
                if gene == 'icaR': icar_truncation = trunc_str
                continue 

            display_str = gene
            if hit.pident < 100.0:
                if ref_prot and prot_clean == ref_prot: display_str += "^"
                else: display_str += "*"
            if hit.coverage < 100.0: display_str += "?"

            found_genes_functional.add(gene)
            
            if gene in self.clf_targets: clf_list.append(display_str)
            elif gene in self.fnb_targets: fnb_list.append(display_str)
            elif gene in self.ica_targets: ica_list.append(display_str)

        status_clf = self.get_group_status(found_genes_functional, self.clf_targets)
        status_fnb = self.get_group_status(found_genes_functional, self.fnb_targets)
        status_ica = self.get_group_status(found_genes_functional, self.ica_targets)

        clf_ok = (status_clf == "Complete")
        fnb_ok = (status_fnb == "Complete")
        ica_ok = (status_ica == "Complete")

        score = 0
        
        if fnb_ok and ica_ok:
            score = 3
        elif clf_ok and ica_ok:
            score = 2
        elif clf_ok or fnb_ok or ica_ok:
            score = 1
        else:
            score = 0
        out = defaults.copy()
        
        out["biofilm_score"] = str(score)
        
        out["clfAB"] = status_clf
        out["clf_genes"] = "; ".join(sorted(clf_list)) if clf_list else "-"
        
        out["fnbAB"] = status_fnb
        out["fnb_genes"] = "; ".join(sorted(fnb_list)) if fnb_list else "-"
        
        out["icaADBC"] = status_ica
        out["ica_genes"] = "; ".join(sorted(ica_list)) if ica_list else "-"
        
        out["biofilm_truncated_hits"] = "; ".join(sorted(truncated_list)) if truncated_list else "-"
        
        out["clfA"] = "1" if 'clfA' in found_genes_functional else "0"
        out["clfB"] = "1" if 'clfB' in found_genes_functional else "0"
        out["fnbA"] = "1" if 'fnbA' in found_genes_functional else "0"
        out["fnbB"] = "1" if 'fnbB' in found_genes_functional else "0"
        out["icaA"] = "1" if 'icaA' in found_genes_functional else "0"
        out["icaB"] = "1" if 'icaB' in found_genes_functional else "0"
        out["icaC"] = "1" if 'icaC' in found_genes_functional else "0"
        out["icaD"] = "1" if 'icaD' in found_genes_functional else "0"

        if icar_truncation: out["icaR_mutations"] = icar_truncation

        return out