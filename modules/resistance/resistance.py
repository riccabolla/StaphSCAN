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
    def __init__(self, min_id=90.0, min_cov=80.0):
        self.name = "resistance"
        self.module_dir = Path(__file__).resolve().parent

        self.db_files = {
            "amino": self.module_dir / "data" / "amino_res.fasta",
            "bla": self.module_dir / "data" / "bla_res.fasta",
            "flq": self.module_dir / "data" / "flq_res.fasta",
            "oxa": self.module_dir / "data" / "oxa_res.fasta", 
            "mlsb": self.module_dir / "data" / "mlsb_res.fasta", 
            "rif": self.module_dir / "data" / "rif_res.fasta", 
            "tet": self.module_dir / "data" / "tet_res.fasta",
            "gly": self.module_dir / "data" / "gly_res.fasta",
        }

        self.ref_prot_dict: Dict[str, str] = {}

        self.mutation_targets = ["gyrA", "parC", "gyrB", "23S", "rpoB"]
        self.known_mutations = {
            "gyrA": {84: ('S', ['L']), 88: ('S', ['P'])},
            "parC": {80: ('S', ['F', 'Y']), 84: ('E', ['K', 'G', 'V'])},
            "rpoB": {481: ('H', ['Y'])}, # look at comibinatorial mutations (ex. L466S+H481N) https://doi.org/10.1016/j.jgar.2021.12.005
            "gyrB": {451: ('T', ['S'])},
            "23S": {2576: ('G', ['T']), 2447: ('G', ['T']), 2500: ('T', ['A'])}, # added for linezolid res
        }

        self.mutation_phenotypes = {
            "gyrA": "Fluoroquinolones-R",
            "gyrB": "Fluoroquinolones-R",
            "parC": "Fluoroquinolones-R",
            "rpoB": "Rifampicin-R",
            "23S": "Linezolid-R", # added for linezolid res
        }

        self.min_id = min_id
        self.min_cov = min_cov
        
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

        self.gene_source: Dict[str, str] = {}
        self.load_ref_seqs()

    def load_ref_seqs(self) -> None:
        #for fpath in self.db_files.values():
        for cat_key, fpath in self.db_files.items():
            if not fpath.exists():
                continue
            for record in SeqIO.parse(fpath, "fasta"):
                gene_family = record.id.split("_")[0]
                self.gene_source[gene_family] = cat_key
                if "23S" in record.id: #23s rRNA exception
                    self.ref_prot_dict[record.id] = str(record.seq)
                    continue
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
            score = self.aligner_prot.score(ref, prot) if ref else 0 #fixed bug, it was self.aligner_dna
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
        if fam == "23S":
            aln = self.aligner_dna.align(ref, found)[0] #added to handle 23s rRNA
        else:
            aln = self.aligner_mut.align(ref, found)[0]
        r, f = aln[0], aln[1]
        pos = 0
        muts = []
        targets = self.known_mutations[fam]
        for rc, fc in zip(r, f):
            if rc != "-":
                pos += 1
            if pos in targets and fc != "-" and fc != rc:
                #ref_aa, allowed = targets[pos] #removed to handle 23s rRNA
                ref_base, allowed = targets[pos] #added to handle 23s rRNA
                if fc in allowed:
                    muts.append(f"{ref_base}{pos}{fc}") #modified to handle 23s rRNA
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

        # section for spurious hit
        df = df[
            (df["pident"] >= 80) &
            (df["coverage"] >= 40)
        ]
        
        #df = df[(df["pident"] >= self.min_id) & (df["coverage"] >= self.min_cov)].copy()
        # end of the section for spurious hit
        if df.empty:
            return df
            
        df["family"] = df["qseqid"].str.split("_").str[0]
        #return df.sort_values("bitscore", ascending=False)
        return df.sort_values(["pident", "bitscore"], ascending=[False, False]) #first id then bitscore

    def make_output(self) -> Dict[str, str]:
        keys = [
            "truncated_resistance_hits",
            "spurious_resistance_hits",
            "res_score",
            "res_gene_count",
            "res_class_count",
            "Amino_res",
            "Bla_res",
            "Flq_res",
            "Gly_res",
            "Mec_res",
            "MLSB_res",
            "Oxa_res",
            "Rif_res",
            "Tet_res",       
            ]
        return {k: "-" for k in keys}

    def run(self, assembly: Path) -> Dict[str, str]:
        out = self.make_output()
        if not self.check_db():
            return out

        seqs = SeqIO.to_dict(SeqIO.parse(assembly, "fasta"))
        df = self.filter_hits(self.run_blast(assembly))
        if df.empty:
            out["res_gene_count"] = "0"
            out["res_class_count"] = "0"
            return out

        best = df.drop_duplicates("family")

        strong, muts, trunc, spur = [], [], [], []
        cat_amino, cat_mec, cat_bla, cat_fq, cat_oxa, cat_mlsb,cat_tet, cat_gly, cat_rif = [], [], [], [], [], [], [], [], []
        mec_aa_found, mec_aa_ref = [], []
        #section mofified for spurious hit
        # for _, r in best.iterrows():
        #     gene_family = r["qseqid"].split("_")[0]
        #     pid = r["pident"]
        #     cov = r["coverage"]
        #     is_strong = (pid >= self.min_id) and (cov >= self.min_cov)

        #     if not is_strong:
        #         spur.append(gene_family)
        #         continue
        for _, r in best.iterrows():
            hit = GeneHit(
                qseqid=r.qseqid, sseqid=r.sseqid, pident=r.pident,
                length=r.length, slen=r.slen, qlen=r.qlen,
                sstart=r.sstart, send=r.send, bitscore=r.bitscore
            )
            dna = self.extract_gene(seqs, hit.sseqid, hit.sstart, hit.send)
            ref = self.ref_prot_dict.get(hit.qseqid, "")
            #found = self.trim_to_ref(self.best_translation(dna, ref), ref) #removed in 23s rRNA implementation
            if hit.family == "23S":
                found = str(dna) #added in 23s rRNA implementation
            else:
                found = self.trim_to_ref(self.best_translation(dna, ref), ref)

            display_str = hit.family
            #is_synonymous = False

            if hit.pident < 100:
                if ref and found == ref:
                    #is_synonymous = True
                    display_str += "^"
                else:
                    display_str += "*"
            
            if hit.coverage < 100:
                display_str += "?"

            is_strong = (hit.pident >= self.min_id) and (hit.coverage >= self.min_cov)
            if not is_strong:
                if hit.family in self.mutation_targets:
                    continue
                spur.append(display_str)
                continue

            if hit.family in self.mutation_targets and ref:
                mm = self.check_mutations(hit.family, found, ref)
                if mm:
                    mut_str = f"{hit.family} [{','.join(mm)}]"
                    muts.append(mut_str)
                    if hit.family in ["gyrA", "parC", "gyrB"]:
                        cat_fq.append(mut_str)
                    elif hit.family == "23S": # linezolid
                        cat_oxa.append(mut_str)
                    else:
                        cat_rif.append(mut_str) #rpoB
                continue

            if "*" in found:
                pct = int((found.find("*") / len(ref)) * 100) if ref else 0
                trunc.append(f"{hit.family}-{pct}%")
                continue

            strong.append(display_str)

            is_mec = hit.family in ["mecA", "mecC"]
            # is_perfect = (hit.pident == 100 and hit.coverage == 100)
            # is_valid = is_perfect or is_synonymous or is_mec
            # if not is_valid:
            #    spur.append(display_str)
            #    continue
            #end of the section modified for spurious hit

            if is_mec:
                #cat_mec.append(display_str)
                mec_aa_found.append(f"{hit.family}_Found:{found}")
                mec_aa_ref.append(f"{hit.family}_Ref:{ref}")
            # elif hit.family == "blaZ":
            #     cat_bla.append(display_str)
            # elif hit.family in ["AAC(6')-Ie-APH(2'')-Ia"]:
            #     cat_amino.append(display_str)    #aminoglycosides added
            # elif hit.family in ["cfrA"]:
            #     cat_lin.append(display_str) # lineozolid added
            # elif hit.family.startswith("tet"):
            #     cat_tet.append(display_str)
            # elif hit.family == "vanA":
            #     cat_van.append(display_str)
            # else:
            #     cat_rif.append(display_str)
            source = self.gene_source.get(hit.family, "other")
            if source == "mec_res": 
                 cat_mec.append(display_str)
            elif is_mec: 
                 cat_mec.append(display_str)
            elif source == "bla":
                 cat_bla.append(display_str)
            elif source == "amino":
                 cat_amino.append(display_str)
            elif source == "oxa":
                 cat_oxa.append(display_str)
            elif source == "tet":
                 cat_tet.append(display_str)
            elif source == "gly":
                 cat_gly.append(display_str)
            elif source == "mlsb":
                 cat_mlsb.append(display_str)
            elif source == "rif":
                 cat_rif.append(display_str)
            elif source == "flq":
                 cat_fq.append(display_str)


        classes_map = {
            "Amino_res": cat_amino,
            "Bla_res": cat_bla,
            "Flq_res": cat_fq,
            "Gly_res": cat_gly,            
            "Mec_res": cat_mec,
            "MLSB_res": cat_mlsb,
            "Oxa_res": cat_oxa,
            "Rif_res": cat_rif,
            "Tet_res": cat_tet,
        }
        
        active_classes = sum(1 for lst in classes_map.values() if len(lst) > 0)

        # classes_to_check = [cat_amino, cat_mec, cat_bla, cat_fq, cat_lin, cat_mlsb,cat_tet, cat_van, cat_rif]
        # class_count = sum(1 for c in classes_to_check if len(c) > 0)

        #out["res_gene_count"] = str(len(strong))
        #out["res_class_count"] = str(class_count)
        out["res_gene_count"] = str(len(strong))
        out["res_class_count"] = str(active_classes)
        out["truncated_resistance_hits"] = "; ".join(trunc) if trunc else "-"
        out["spurious_resistance_hits"] = "; ".join(spur) if spur else "-"
        out["Amino_res"] = "; ".join(cat_amino) if cat_amino else "-" #aminoglycosides
        out["Mec_res"] = "; ".join(cat_mec) if cat_mec else "-"
        out["Bla_res"] = "; ".join(cat_bla) if cat_bla else "-"
        out["Flq_res"] = "; ".join(cat_fq) if cat_fq else "-"
        out["Oxa_res"] = "; ".join(cat_oxa) if cat_oxa else "-" # linezolid
        out["MLSB_res"] = "; ".join(cat_mlsb) if cat_mlsb else "-" # mlsb
        out["Tet_res"] = "; ".join(cat_tet) if cat_tet else "-"
        out["Gly_res"] = "; ".join(cat_gly) if cat_gly else "-"
        out["Rif_res"] = "; ".join(cat_rif) if cat_rif else "-"        
        #out["Mec_AA_Found"] = " | ".join(mec_aa_found) if mec_aa_found else "-"
        #out["Mec_AA_Ref"] = " | ".join(mec_aa_ref) if mec_aa_ref else "-"

        score = 0
        is_van_positive = any("vanA" in x or "vanB" in x for x in cat_gly)
        if is_van_positive:
            score = 3
        elif cat_mec:
            score = 2
        elif cat_bla:
            score = 1
            
        out["res_score"] = str(score)

        return out