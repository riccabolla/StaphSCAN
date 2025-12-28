import pandas as pd
import subprocess
import io
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple
from Bio.Seq import Seq
from Bio import SeqIO

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
        return self.qseqid.split("_")[0]

class Module:
    def __init__(self, min_id=90.0, min_cov=80.0):
        self.name = "virulence"
        self.module_dir = Path(__file__).resolve().parent
        self.data_dir = self.module_dir / "data"
        self.db_fasta = self.data_dir / "targets.fasta"
        
        self.min_id = min_id
        self.min_cov = min_cov

    def check_db(self) -> bool:
        return self.db_fasta.exists()

    def extract_and_translate(self, seqs: Dict[str, Seq], hit: GeneHit) -> Tuple[str, int]:
        """
        Extracts DNA based on BLAST coords and translates immediately.
        """
        if hit.sseqid not in seqs: return "", 0
        contig_seq = seqs[hit.sseqid].seq
        
        if hit.sstart < hit.send:
            dna = contig_seq[hit.sstart-1 : hit.send]
        else:
            dna = contig_seq[hit.send-1 : hit.sstart].reverse_complement()

        dna_str = str(dna).replace("-", "")
        
        remainder = len(dna_str) % 3
        if remainder > 0:
            dna_str = dna_str[:-remainder]
            
        if not dna_str: return "", 0

        try:
            prot = str(Seq(dna_str).translate(table=11))
            if "*" in prot:
                prot = prot.split("*")[0]
            return prot, len(prot)
        except Exception:
            return "X", 0

    def run(self, assembly_path: Path) -> Dict[str, str]:
        results = {"vir_pvl": "-", "vir_tsst": "-", "vir_et": "-", "spurious_virulence_hits": "-"}
        
        if not self.check_db():
            return results

        cmd = [
            "blastn", "-task", "blastn",
            "-query", str(self.db_fasta),
            "-subject", str(assembly_path),
            "-outfmt", "6 qseqid sseqid pident length slen qlen sstart send bitscore"
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: return results
            
            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", 
                             names=["qseqid", "sseqid", "pident", "length", "slen", "qlen", "sstart", "send", "bitscore"])
            
            df['coverage'] = (df['length'] / df['qlen']) * 100
            df = df[(df['pident'] >= 80.0) & (df['coverage'] >= 40.0)]
            
            if df.empty: return results
            
            seqs = SeqIO.to_dict(SeqIO.parse(assembly_path, "fasta"))
            
            strong_hits = []
            spurious_hits = []
            
            df['family'] = df['qseqid'].apply(lambda x: x.split('_')[0])
            best_hits = df.sort_values("bitscore", ascending=False).drop_duplicates("family")

            for _, row in best_hits.iterrows():
                hit_data = row.drop(['family', 'coverage'], errors='ignore').to_dict()
                hit = GeneHit(**hit_data)
                
                is_strong = (hit.pident >= self.min_id) and (hit.coverage >= self.min_cov)
                display_name = hit.family 
                
                if not is_strong:
                    tag = f"(Id:{int(hit.pident)}% Cov:{int(hit.coverage)}%)"
                    spurious_hits.append(f"{display_name}{tag}")
                    continue
                
                prot_seq, prot_len = self.extract_and_translate(seqs, hit)
                
                annotation = ""
                if hit.pident < 100.0: annotation += "*"
                
                expected_len = hit.qlen / 3
                if prot_len < (expected_len * 0.9):
                    annotation += "^"
                
                strong_hits.append(f"{display_name}{annotation}")

            def find_genes(search_terms):
                found = []
                for term in search_terms:
                    found.extend([x for x in strong_hits if term.lower() in x.lower()])
                return sorted(list(set(found)))

            lukS = find_genes(["luks"])
            lukF = find_genes(["lukf"])
            
            if lukS and lukF:
                results["vir_pvl"] = "Positive"
            elif lukS or lukF:
                present = lukS + lukF
                results["vir_pvl"] = f"Partial ({', '.join(present)})"
            
            tsst = find_genes(["tst"])
            results["vir_tsst"] = "; ".join(tsst) if tsst else "-"
            
            et = find_genes(["eta", "etb"])
            results["vir_et"] = "; ".join(et) if et else "-"
            
            results["spurious_virulence_hits"] = "; ".join(spurious_hits) if spurious_hits else "-"
            
            return results

        except Exception as e:
            print(f"Error in virulence module: {e}", file=sys.stderr)
            results["vir_pvl"] = "Error"
            return results