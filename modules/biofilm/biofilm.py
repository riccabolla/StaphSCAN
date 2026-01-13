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

# @dataclass
# class BiofilmHit:
#     qseqid: str
#     sseqid: str
#     pident: float
#     length: int
#     qlen: int
#     sstart: int
#     send: int
#     bitscore: float

#     @property
#     def coverage_pct(self) -> float:
#         return (self.length / self.qlen) * 100 if self.qlen else 0.0

#     @property
#     def gene_name(self) -> str:
#         return self.qseqid.split("_")[0]

class Module:
    def __init__(self, min_id=90.0, min_cov=80.0):
        self.name = "biofilm"
        self.module_dir = Path(__file__).resolve().parent
        self.db_fasta = self.module_dir / "data" / "biofilm.fasta"

        self.clf_targets = ["clfA", "clfB"]
        self.fnb_targets = ["fnbA", "fnbB"]
        self.ica_targets = ["icaA", "icaB", "icaC", "icaD"]
        self.regulator = "icaR"  

        self.min_identity = float(min_id)
        self.min_coverage = float(min_cov)

        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = "global"
        self.aligner.match_score = 1
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5

        self.ref_prot_dict: Dict[str, str] = {}
        self.load_ref_seqs()

    def load_ref_seqs(self):
        if not self.db_fasta.exists(): return
        for rec in SeqIO.parse(self.db_fasta, "fasta"):
            prot = str(rec.seq.translate(table=11)).rstrip("*")
            self.ref_prot_dict[rec.id] = prot

    def extract_gene(self, seqs, contig, start, end):
        if contig not in seqs: return Seq("")
        s = seqs[contig].seq
        return s[start - 1:end] if start < end else s[end - 1:start].reverse_complement()


    def translate_frames(self, dna: Seq) -> List[str]:
        prots = []
        for frame in range(3):
            sub = dna[frame:]
            trim = len(sub) % 3
            if trim > 0: sub = sub[:-trim]
            if sub:
                prots.append(str(sub.translate(table=11)))
        return prots

    def best_translation(self, dna: Seq, ref_prot: str) -> str:
        best_cand = ""
        best_score = -float("inf")
        
        check_len = min(len(ref_prot), 15)
        ref_head = ref_prot[:check_len]

        for prot in self.translate_frames(dna):
            cand = prot.rstrip("*")
            if not cand: continue

            try:
                score = self.aligner.score(ref_prot, cand)
            except: score = 0
            
            start_bonus = 0
            cand_head = cand[:check_len]
            if len(cand_head) == check_len:
                matches = sum(1 for r, c in zip(ref_head, cand_head) if r == c)
                if matches / check_len >= 0.6: 
                    start_bonus = 1000  

            final_score = score + start_bonus
            
            if final_score > best_score:
                best_score = final_score
                best_cand = cand

        return best_cand


    def detect_truncation(self, ref: str, qry: str):
        if not qry: return None

        try:
            aln = self.aligner.align(ref, qry)[0]
        except: return None
        
        r_aln, q_aln = aln[0], aln[1]

        aa_pos = 0 
        
        for ra, qa in zip(r_aln, q_aln):
            if ra != "-": 
                aa_pos += 1
            
            if qa == "*":
                if aa_pos > (len(ref) * 0.90):
                    return None
                pct = int((aa_pos / len(ref)) * 100)
                return pct if pct > 0 else 1 

        return None

    def run(self, assembly_path: Path) -> Dict[str, str]:
        out = {k: "0" for k in [
            "clfA", "clfB", "fnbA", "fnbB", "icaA", "icaB", "icaC", "icaD"
        ]}
        out.update({
            "biofilm_score": "0",
            "biofilm_truncated_hits": "-",
            "biofilm_spurious_hits": "-",
            "clfAB": "-", "fnbAB": "-", "icaADBC": "-",
            "icaR_mutations": "-" 
        })

        if not assembly_path.exists(): return out
        
        try:
            seqs = SeqIO.to_dict(SeqIO.parse(assembly_path, "fasta"))
        except: return out

        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "db"
            subprocess.run(
                ["makeblastdb", "-in", str(assembly_path), "-dbtype", "nucl", "-out", str(db)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
            )

            cmd = [
                "blastn", "-query", str(self.db_fasta), "-db", str(db),
                "-outfmt", "6 qseqid sseqid pident length qlen sstart send bitscore"
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: return out

        try:
            df = pd.read_csv(
                io.StringIO(res.stdout), sep="\t",
                names=["qseqid", "sseqid", "pident", "length", "qlen", "sstart", "send", "bitscore"]
            )
        except: return out

        df["coverage_pct"] = (df["length"] / df["qlen"]) * 100

        df = df[
            (df["pident"] >= 80) & # added min id for spurious hits
            (df["coverage_pct"] >= 40) # added min cov for spurious hits    
        ]

        if df.empty: return out

        df["gene"] = df["qseqid"].str.split("_").str[0]
        df = df.sort_values(["pident", "bitscore"], ascending=[False, False]) #identity first
        df = df.drop_duplicates("gene")

        truncated = []
        spurious = []
        found_genes = set()
        gene_display_map = {}
        icar_trunc_str = None

        for _, r in df.iterrows():
            gene = r["gene"]
            pid = r["pident"]
            cov = r["coverage_pct"]

            is_strong = (pid >= self.min_identity) and (cov >= self.min_coverage)

            if not is_strong:
                tag = f"{gene}"
                spurious.append(tag)
                continue
            
            dna = self.extract_gene(seqs, r["sseqid"], r["sstart"], r["send"])
            ref = self.ref_prot_dict.get(r["qseqid"], "")
            prot = self.best_translation(dna, ref)
            
            trunc_pct = self.detect_truncation(ref, prot)
            
            if trunc_pct is not None:
                t_str = f"{gene}-{trunc_pct}%"
                truncated.append(t_str)
                
                # icaR specific truncations
                if gene == self.regulator:
                    icar_trunc_str = t_str
                
                continue

            found_genes.add(gene)
            out[gene] = "1"

            display_str = gene

            if pid < 100.0:
                if ref and prot == ref:
                    display_str += "^"
                else:
                    display_str += "*"
            
            if cov < 100.0:
                display_str += "?"

            gene_display_map[gene] = display_str

        out["biofilm_truncated_hits"] = "; ".join(truncated) if truncated else "-"
        out["biofilm_spurious_hits"] = "; ".join(spurious) if spurious else "-"
        
        if icar_trunc_str:
            out["icaR_mutations"] = icar_trunc_str

        # genes
        out["clf_genes"] = ";".join(sorted([gene_display_map[g] for g in found_genes if g in self.clf_targets])) or "-"
        out["fnb_genes"] = ";".join(sorted([gene_display_map[g] for g in found_genes if g in self.fnb_targets])) or "-"
        out["ica_genes"] = ";".join(sorted([gene_display_map[g] for g in found_genes if g in self.ica_targets])) or "-"
        #completeness    
        def status(targets):
            n = sum(g in found_genes for g in targets)
            return "Complete" if n == len(targets) else "Incomplete" if n else "-"

        out["clfAB"] = status(self.clf_targets)
        out["fnbAB"] = status(self.fnb_targets)
        out["icaADBC"] = status(self.ica_targets)

        has_ica = all(g in found_genes for g in self.ica_targets)
        
        has_clf = any(g in found_genes for g in self.clf_targets)
        has_fnb = any(g in found_genes for g in self.fnb_targets)

        # new score system
        # it now considers only the completeness of icaAD relevant.
        # the presence of just one clf or fnb gene is sufficient to increase the score
        # the new score reflects just the biofilm potential, not prediction on the strength
        score = 0

        if has_ica:
            score = 1            
            if has_clf and has_fnb:
                score = 3
            elif has_clf or has_fnb:
                score = 2
        
        out["biofilm_score"] = str(score)

        return out
    
        #old score
        #score = 0
        #if out["fnbAB"] == "Complete" and out["icaADBC"] == "Complete":
        #    score = 3
        #elif out["clfAB"] == "Complete" and out["icaADBC"] == "Complete":
        #    score = 2
        #elif any(out[x] == "Complete" for x in ("clfAB", "fnbAB", "icaADBC")):
        #    score = 1    
        #out["biofilm_score"] = str(score)
        #return out