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

    def annotation_tag(self) -> str:
        return f"(Id:{int(self.pident)}% Cov:{int(self.coverage)}%)"


class Module:
    def __init__(self):
        self.name = "resistance"
        self.module_dir = Path(__file__).resolve().parent

        self.db_files = {
            "bla": self.module_dir / "data" / "bla_res.fasta",
            "flq": self.module_dir / "data" / "flq_res.fasta",
            "tet": self.module_dir / "data" / "tet_res.fasta",
            "van": self.module_dir / "data" / "van_res.fasta",
        }

        self.ref_prot_dict: Dict[str, str] = {}

        self.mutation_targets = ["gyrA", "parC", "grlA", "gyrB", "rpoB"]
        self.known_mutations = {
            "gyrA": {84: ('S', ['L']), 88: ('S', ['P'])},
            "parC": {80: ('S', ['F', 'Y']), 84: ('E', ['K', 'G', 'V'])},
            "grlA": {80: ('S', ['F', 'Y'])},
            "rpoB": {481: ('H', ['N', 'Y']), 527: ('I', ['M'])},
            "gyrB": {451: ('T', ['S'])},
        }

        self.mutation_phenotypes = {
            "gyrA": "Fluoroquinolones-R",
            "gyrB": "Fluoroquinolones-R",
            "parC": "Fluoroquinolones-R",
            "grlA": "Fluoroquinolones-R",
            "rpoB": "Rifampicin-R",
        }

        self.aligner_dna = Align.PairwiseAligner()
        self.aligner_dna.mode = "global"
        self.aligner_dna.match_score = 1
        self.aligner_dna.mismatch_score = -1
        self.aligner_dna.open_gap_score = -2
        self.aligner_dna.extend_gap_score = -0.5

        self.aligner_prot = Align.PairwiseAligner()
        self.aligner_prot.mode = "global"
        self.aligner_prot.match_score = 5
        self.aligner_prot.mismatch_score = -4

        self.aligner_mut = Align.PairwiseAligner()
        self.aligner_mut.mode = "global"
        self.aligner_mut.match_score = 5
        self.aligner_mut.mismatch_score = -4
        self.aligner_mut.open_gap_score = -10
        self.aligner_mut.extend_gap_score = -0.5

        self.load_ref_seqs()

    def load_ref_seqs(self) -> None:
        for fpath in self.db_files.values():
            if not fpath.exists():
                continue
            for record in SeqIO.parse(fpath, "fasta"):
                seq = record.seq
                remainder = len(seq) % 3
                if remainder > 0:
                    seq = seq[:-remainder]
                
                prot = str(seq.translate(table=11)).strip("*")
                mpos = prot.find("M")
                self.ref_prot_dict[record.id] = prot[mpos:] if mpos != -1 else prot

    def check_db(self) -> bool:
        return all(f.exists() for f in self.db_files.values())

    def extract_gene(self, seqs: Dict[str, Seq], contig: str, start: int, end: int) -> Seq:
        if contig not in seqs:
            return Seq("")
        seq = seqs[contig].seq
        return seq[start - 1:end] if start < end else seq[end - 1:start].reverse_complement()

    def best_translation(self, dna: Seq, ref: str) -> str:
        best_score = -1e9
        best = ""
        for frame in range(3):
            frag = dna[frame:]
            frag = frag[:len(frag) - len(frag) % 3]
            if not frag:
                continue
            prot = str(frag.translate(table=11)).strip("*")
            score = self.aligner_dna.score(ref, prot) if ref else 0
            if score > best_score:
                best_score = score
                best = prot
        return best

    def trim_to_ref(self, found: str, ref: str) -> str:
        if not found or not ref:
            return found
        aln = self.aligner_prot.align(ref, found)[0]
        r, f = aln[0], aln[1]
        idx = 0
        for rc, fc in zip(r, f):
            if rc != "-":
                break
            if fc != "-":
                idx += 1
        return found[idx:]

    def check_mutations(self, fam: str, found: str, ref: str) -> List[str]:
        if fam not in self.known_mutations:
            return []
        aln = self.aligner_mut.align(ref, found)[0]
        r, f = aln[0], aln[1]
        pos = 0
        muts = []
        targets = self.known_mutations[fam]
        for rc, fc in zip(r, f):
            if rc != "-":
                pos += 1
            if pos in targets and fc != "-" and fc != rc:
                ref_aa, allowed = targets[pos]
                if fc in allowed:
                    muts.append(f"{ref_aa}{pos}{fc}")
        return muts

    def run_blast(self, assembly: Path) -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "db"
            subprocess.run(
                ["makeblastdb", "-in", str(assembly), "-dbtype", "nucl", "-out", str(db)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )

            query = Path(td) / "query.fasta"
            with open(query, "w") as out:
                for f in self.db_files.values():
                    out.write(f.read_text() + "\n")

            cmd = [
                "blastn",
                "-query", str(query),
                "-db", str(db),
                "-outfmt", "6 qseqid sseqid pident length slen qlen sstart send bitscore",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout:
                return pd.DataFrame()

            return pd.read_csv(
                io.StringIO(res.stdout),
                sep="\t",
                names=["qseqid", "sseqid", "pident", "length", "slen", "qlen", "sstart", "send", "bitscore"],
            )

    def filter_hits(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        df = df.copy()
        df["coverage"] = (df["length"] / df["qlen"]) * 100
        
        df = df[(df["pident"] >= 90) & (df["coverage"] >= 80)].copy()
        if df.empty:
            return df
            
        df["family"] = df["qseqid"].str.split("_").str[0]
        return df.sort_values("bitscore", ascending=False)

    def make_output(self) -> Dict[str, str]:
        keys = [
            "res_genes",
            "res_mutations",
            "truncated_resistance_hits",
            "spurious_resistance_hits",
            "res_score",
            "Mec_RES",
            "Beta_lactamases",
            "Fluoroquinolones",
            "Tetracyclines",
            "Vancomycin",
            "Other_RES",
        ]
        return {k: "-" for k in keys}

    def run(self, assembly: Path) -> Dict[str, str]:
        out = self.make_output()
        if not self.check_db():
            return out

        seqs = SeqIO.to_dict(SeqIO.parse(assembly, "fasta"))
        df = self.filter_hits(self.run_blast(assembly))
        if df.empty:
            return out

        best = df.drop_duplicates("family")

        strong, muts, trunc, spur = [], [], [], []
        cat_mec, cat_bla, cat_fq, cat_tet, cat_van, cat_other = [], [], [], [], [], []
        mec_aa_found, mec_aa_ref = [], []

        for _, r in best.iterrows():
            hit = GeneHit(
                qseqid=r.qseqid, sseqid=r.sseqid, pident=r.pident,
                length=r.length, slen=r.slen, qlen=r.qlen,
                sstart=r.sstart, send=r.send, bitscore=r.bitscore
            )
            dna = self.extract_gene(seqs, hit.sseqid, hit.sstart, hit.send)
            ref = self.ref_prot_dict.get(hit.qseqid, "")
            found = self.trim_to_ref(self.best_translation(dna, ref), ref)

            if hit.family in self.mutation_targets and ref:
                mm = self.check_mutations(hit.family, found, ref)
                if mm:
                    mut_str = f"{hit.family} [{','.join(mm)}]"
                    muts.append(mut_str)
                    if hit.family in ["gyrA", "parC", "grlA", "gyrB"]:
                        cat_fq.append(mut_str)
                    else:
                        cat_other.append(mut_str)
                continue

            if "*" in found:
                pct = int((found.find("*") / len(ref)) * 100) if ref else 0
                trunc.append(f"{hit.family}-{pct}%")
                continue

            is_synonymous = False
            display_str = hit.family

            if hit.pident < 100:
                if ref and found == ref:
                    is_synonymous = True
                    display_str += "^"
                else:
                    display_str += "*"
            
            if hit.coverage < 100:
                display_str += "?"

            is_mec = hit.family in ["mecA", "mecC"]
            is_perfect = (hit.pident == 100 and hit.coverage == 100)
            is_valid = is_perfect or is_synonymous or is_mec

            if not is_valid:
                spur.append(display_str)
                continue

            strong.append(display_str)

            if is_mec:
                cat_mec.append(display_str)
                mec_aa_found.append(f"{hit.family}_Found:{found}")
                mec_aa_ref.append(f"{hit.family}_Ref:{ref}")
            elif hit.family == "blaZ":
                cat_bla.append(display_str)
            elif hit.family.startswith("tet"):
                cat_tet.append(display_str)
            elif hit.family == "vanA":
                cat_van.append(display_str)
            else:
                cat_other.append(display_str)

        out["res_genes"] = "; ".join(strong) if strong else "-"
        out["res_mutations"] = "; ".join(muts) if muts else "-"
        out["truncated_resistance_hits"] = "; ".join(trunc) if trunc else "-"
        out["spurious_resistance_hits"] = "; ".join(spur) if spur else "-"
        out["Mec_RES"] = "; ".join(cat_mec) if cat_mec else "-"
        out["Beta_lactamases"] = "; ".join(cat_bla) if cat_bla else "-"
        out["Fluoroquinolones"] = "; ".join(cat_fq) if cat_fq else "-"
        out["Tetracyclines"] = "; ".join(cat_tet) if cat_tet else "-"
        out["Vancomycin"] = "; ".join(cat_van) if cat_van else "-"
        out["Other_RES"] = "; ".join(cat_other) if cat_other else "-"
        
        out["Mec_AA_Found"] = " | ".join(mec_aa_found) if mec_aa_found else "-"
        out["Mec_AA_Ref"] = " | ".join(mec_aa_ref) if mec_aa_ref else "-"

        score = 0
        if any("vanA" in g for g in cat_van):
            score = 3
        elif cat_mec:
            score = 2
        elif cat_bla:
            score = 1
            
        out["res_score"] = str(score)

        return out