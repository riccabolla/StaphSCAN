import pandas as pd
import subprocess
import io
import sys
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
        return self.qseqid.split("_")[0]

class Module:
    def __init__(self, min_id=90.0, min_cov=80.0):
        self.name = "virulence"
        self.module_dir = Path(__file__).resolve().parent
        self.data_dir = self.module_dir / "data"
        self.db_fasta = self.data_dir / "targets.fasta"
        
        self.min_id = min_id
        self.min_cov = min_cov

        self.aligner_prot = Align.PairwiseAligner()
        self.aligner_prot.mode = "global"
        self.aligner_prot.match_score = 5
        self.aligner_prot.mismatch_score = -4

        self.ref_prot_dict: Dict[str, str] = {}
        self.load_ref_seqs()

    def load_ref_seqs(self):
        if not self.db_fasta.exists(): return
        for rec in SeqIO.parse(self.db_fasta, "fasta"):
            seq = rec.seq
            remainder = len(seq) % 3
            if remainder > 0:
                seq = seq[:-remainder]
            prot = str(seq.translate(table=11)).strip("*")
            self.ref_prot_dict[rec.id] = prot

    def check_db(self) -> bool:
        return self.db_fasta.exists()

    def extract_gene(self, seqs: Dict[str, Seq], contig: str, start: int, end: int) -> Seq:
        if contig not in seqs: return Seq("")
        s = seqs[contig].seq
        return s[start - 1:end] if start < end else s[end - 1:start].reverse_complement()

    def best_translation(self, dna: Seq, ref_prot: str) -> str:
        best_cand = ""
        best_score = -float("inf")
        
        for frame in range(3):
            sub = dna[frame:]
            trim = len(sub) % 3
            if trim > 0: sub = sub[:-trim]
            if not sub: continue
            
            cand = str(sub.translate(table=11)).strip("*")
            score = self.aligner_prot.score(ref_prot, cand) if ref_prot else 0
            
            if score > best_score:
                best_score = score
                best_cand = cand
        return best_cand

    def trim_to_ref(self, found: str, ref: str) -> str:
        if not found or not ref: return found
        aln = self.aligner_prot.align(ref, found)[0]
        r, f = aln[0], aln[1]
        idx = 0
        for rc, fc in zip(r, f):
            if rc != "-": break
            if fc != "-": idx += 1
        return found[idx:]

    def run(self, assembly_path: Path) -> Dict[str, str]:
        results = {
            "vir_score": "-", "vir_pvl": "-", "vir_tsst": "-", "vir_et": "-", "vir_lukED": "-", "vir_se": "-",
            "spurious_virulence_hits": "-", "truncated_virulence_hits": "-"
        }
        
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
            
            if df.empty: 
                results["vir_score"] = 0
                return results
            
            seqs = SeqIO.to_dict(SeqIO.parse(assembly_path, "fasta"))
            
            strong_hits = []
            spurious_hits = []
            truncated_hits = []
            
            df['family'] = df['qseqid'].apply(lambda x: x.split('_')[0])

            best_hits = df.sort_values(["pident", "bitscore"], ascending=[False, False]).drop_duplicates("family")

            for _, row in best_hits.iterrows():
                hit_data = row.drop(['family', 'coverage'], errors='ignore').to_dict()
                hit = GeneHit(**hit_data)
                
                dna = self.extract_gene(seqs, hit.sseqid, hit.sstart, hit.send)
                ref = self.ref_prot_dict.get(hit.qseqid, "")
                prot = self.trim_to_ref(self.best_translation(dna, ref), ref)

                is_truncated = False
                trunc_pct = 0
                
                if "*" in prot:
                    is_truncated = True
                    trunc_pct = int((prot.find("*") / len(ref)) * 100) if ref else 0
                else:
                     expected_len = hit.qlen / 3
                     if len(prot) < (expected_len * 0.9):
                         is_truncated = True
                         trunc_pct = int((len(prot) / expected_len) * 100)

                if is_truncated:
                    truncated_hits.append(f"{hit.family}-{trunc_pct}%")
                    continue 

                display_str = hit.family
                
                if hit.pident < 100.0:
                    if ref and prot == ref:
                        display_str += "^"
                    else:
                        display_str += "*"
                
                if hit.coverage < 100.0:
                    display_str += "?"

                is_strong = (hit.pident >= self.min_id) and (hit.coverage >= self.min_cov)
                
                if is_strong:
                    strong_hits.append(display_str)
                else:
                    spurious_hits.append(display_str)

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
            
            tsst = find_genes(["tsst1"])
            results["vir_tsst"] = "; ".join(tsst) if tsst else "-"
            
            et = find_genes(["eta", "etb"])
            results["vir_et"] = "; ".join(et) if et else "-"

            lukE = find_genes(["luke"])
            lukD = find_genes(["lukd"])
            if lukE and lukD:
                results["vir_lukED"] = "Positive"
            elif lukE or lukD:
                present = lukE + lukD
                results["vir_lukED"] = f"Partial ({', '.join(present)})"

            #added se genes block    
            se_genes = ["sea", "sec", "seh", "selk", "sell", "selq"]
            se_found = []
            for gene in se_genes:
                se_found.extend(find_genes([gene]))
            se_found = sorted(list(set(se_found)))
            results["vir_se"] = "; ".join(se_found) if se_found else "-"    
            #end of se genes block

            results["spurious_virulence_hits"] = "; ".join(spurious_hits) if spurious_hits else "-"
            results["truncated_virulence_hits"] = "; ".join(truncated_hits) if truncated_hits else "-"
            
            #score logic block
            clean_hits = set()
            for h in strong_hits:
                clean_name = h.replace("*", "").replace("^", "").replace("?", "").lower()
                clean_hits.add(clean_name)

            score = 0
            if "lukf" in clean_hits and "luks" in clean_hits:
                score = 4
            elif "tsst1" in clean_hits:
                score = 3
            elif "eta" in clean_hits or "etb" in clean_hits:
                score = 2
            else:
                has_se = not clean_hits.isdisjoint(se_genes)
                has_luk = "lukd" in clean_hits and "luke" in clean_hits
                
                if has_se or has_luk:
                    score = 1

            results["vir_score"] = score
            #end of score logic block

            return results

        except Exception as e:
            print(f"Error in virulence module: {e}", file=sys.stderr)
            results["vir_pvl"] = "Error"
            return results