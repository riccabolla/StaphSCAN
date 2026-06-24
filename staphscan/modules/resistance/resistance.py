import pandas as pd
import subprocess
import io
import tempfile
import json
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


class Module:
    def __init__(self, min_id=90.0, min_cov=80.0):
        self.name = "resistance"
        self.module_dir = Path(__file__).resolve().parent
        self.data_dir = self.module_dir / "data"
        def get_db(prefix):
            try:
                return next(self.data_dir.glob(f"{prefix}*.fasta"))
            except StopIteration:
                raise FileNotFoundError(f"Database missing: No file matching '{prefix}*.fasta' found.")

        self.db_files = {
            "amino": get_db("amino_res"),
            "bla":   get_db("bla_res"),
            "flq":   get_db("flq_res"),
            "oxa":   get_db("oxa_res"), 
            "mlsb":  get_db("mlsb_res"), 
            "rif":   get_db("rif_res"), 
            "tet":   get_db("tet_res"),
            "gly":   get_db("gly_res"),
        }

        self.ref_prot_dict: Dict[str, str] = {}

        try:
            mutations_file = next(self.data_dir.glob("mutation*.json"))
        except StopIteration:
            raise FileNotFoundError(f"Database missing: No file matching 'mutation*.json' found.")
        
        with open(mutations_file, "r") as f:
            raw_mutations = json.load(f)
        self.known_mutations = {}
        for gene, mutations in raw_mutations.items():
            self.known_mutations[gene] = {
                int(pos): tuple(data) for pos, data in mutations.items()
            }
        self.mutation_targets = list(self.known_mutations.keys())
        #self.mutation_targets = ["gyrA", "parC", "gyrB", "23S", "rpoB"]
        
        #self.known_mutations = {
        #    "gyrA": {84: ('S', ['L']), 88: ('E', ['K', 'G'])}, 
        #    "parC": {80: ('S', ['F', 'Y']), 84: ('E', ['K', 'G', 'V'])},
        #    "rpoB": {481: ('H', ['Y', 'N']), 466: ('L', ['S']), 473: ('A', ['T']), 477: ('A', ['T'])}, 
        #    "gyrB": {451: ('T', ['S'])},
        #    "23S": {2576: ('G', ['T']), 2447: ('G', ['T']), 2500: ('T', ['A'])}, 
        #}

        self.min_id = min_id
        self.min_cov = min_cov
        
        # Aligner for 23S (DNA)
        self.aligner_dna = Align.PairwiseAligner()
        self.aligner_dna.mode = "global"
        self.aligner_dna.match_score = 1
        self.aligner_dna.mismatch_score = -1
        self.aligner_dna.open_gap_score = -2
        self.aligner_dna.extend_gap_score = -0.5

        # Aligner for Protein
        self.aligner_prot = Align.PairwiseAligner()
        self.aligner_prot.mode = "global"
        self.aligner_prot.match_score = 5
        self.aligner_prot.mismatch_score = -4
        self.aligner_prot.open_gap_score = -10 
        self.aligner_prot.extend_gap_score = -0.5

        self.gene_source: Dict[str, str] = {}
        self.load_ref_seqs()

    def load_ref_seqs(self) -> None:
        """Loads reference sequences, translating NT to AA exactly as provided."""
        for cat_key, fpath in self.db_files.items():
            if not fpath.exists():
                continue
            for record in SeqIO.parse(fpath, "fasta"):
                gene_family = record.id.split("_")[0]
                self.gene_source[gene_family] = cat_key
                
                if "23S" in record.id: 
                    self.ref_prot_dict[record.id] = str(record.seq)
                    continue
                
                seq = record.seq
                remainder = len(seq) % 3
                if remainder > 0:
                    seq = seq[:-remainder]
                
                prot = str(seq.translate(table=11)).strip("*")
                self.ref_prot_dict[record.id] = prot

    def check_db(self) -> bool:
        return all(f.exists() for f in self.db_files.values())

    def extract_gene(self, seqs: Dict[str, Seq], contig: str, start: int, end: int) -> Seq:
        if contig not in seqs:
            return Seq("")
        seq = seqs[contig].seq
        return seq[start - 1:end] if start < end else seq[end - 1:start].reverse_complement()

    def best_translation(self, dna: Seq, ref: str) -> str:
        """Finds the reading frame that best matches the reference protein."""
        best_score = -1e9
        best = ""
        for frame in range(3):
            frag = dna[frame:]
            frag = frag[:len(frag) - len(frag) % 3]
            if not frag:
                continue
            prot = str(frag.translate(table=11)).strip("*")
            score = self.aligner_prot.score(ref, prot) if ref else 0 
            if score > best_score:
                best_score = score
                best = prot
        return best

    def trim_to_ref(self, found: str, ref: str) -> str:
        """Trims leading gaps to align start of found with start of ref."""
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
        """
        Checks for specific point mutations using Alignment-Based Coordinates.
        Handles indels by counting biological reference positions.
        """
        if not found or not ref: # add safety check to prevent crushing
            return []
        
        if fam not in self.known_mutations:
            return []

        # Select Aligner
        if fam == "23S":
            aligner = self.aligner_dna
        else:
            aligner = self.aligner_prot

        # Create Alignment
        alignments = aligner.align(ref, found)
        if not alignments:
            return []
        
        aln = alignments[0]
        ref_aln = aln[0]
        found_aln = aln[1]

        muts = []
        targets = self.known_mutations[fam]
        
        # Iterate through alignment to track position
        ref_biological_pos = 0
        
        for r_char, f_char in zip(ref_aln, found_aln):
            # If r_char is a gap, it's an insertion in query; ref pos doesn't advance
            if r_char == "-":
                continue
            
            ref_biological_pos += 1
            
            if ref_biological_pos in targets:
                expected_ref_base, allowed_muts = targets[ref_biological_pos]

                if r_char != expected_ref_base:
                    # DB sequence doesn't match config. Skip to avoid false positive.
                    continue

                # If query has gap here, do not call mutation
                if f_char == "-":
                    continue
                
                # mutation check:
                if f_char != expected_ref_base and f_char in allowed_muts:
                    muts.append(f"{expected_ref_base}{ref_biological_pos}{f_char}")

        return muts

    def run_blast(self, assembly: Path) -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as td:
            db = Path(td) / "db"
            subprocess.run(
                ["makeblastdb", "-in", str(assembly), "-dbtype", "nucl", "-out", str(db)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
            )

            query = Path(td) / "query.fasta"
            with open(query, "w") as out:
                for f in self.db_files.values():
                    out.write(f.read_text() + "\n")

            cmd = [
                "blastn", "-query", str(query), "-db", str(db),
                "-outfmt", "6 qseqid sseqid pident length slen qlen sstart send bitscore",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout:
                return pd.DataFrame()

            return pd.read_csv(
                io.StringIO(res.stdout),
                sep="\t",
                names=["qseqid", "sseqid", "pident", "length", "slen", "qlen", "sstart", "send", "bitscore"],
                dtype={"qseqid": str, "sseqid": str} # add to fix contig name issues
            )

    def filter_hits(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df
        df = df.copy()
        df["coverage"] = (df["length"] / df["qlen"]) * 100

        # Spurious hit pre-filter
        df = df[(df["pident"] >= 80) & (df["coverage"] >= 40)]
        
        if df.empty:
            return df
            
        df["family"] = df["qseqid"].str.split("_").str[0]
        return df.sort_values(["pident", "bitscore"], ascending=[False, False])

    def make_output(self) -> Dict[str, str]:
        keys = ["truncated_resistance_hits", "spurious_resistance_hits", "res_score",
                "res_gene_count", "res_class_count", "Amino_res", "Bla_res", "Flq_res",
                "Gly_res", "Mec_res", "MLSB_res", "Oxa_res", "Rif_res", "Tet_res"]
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
        cat_amino, cat_mec, cat_bla, cat_fq, cat_oxa, cat_mlsb, cat_tet, cat_gly, cat_rif = [], [], [], [], [], [], [], [], []
        
        for _, r in best.iterrows():
            hit = GeneHit(
                qseqid=r.qseqid, sseqid=r.sseqid, pident=r.pident,
                length=r.length, slen=r.slen, qlen=r.qlen,
                sstart=r.sstart, send=r.send, bitscore=r.bitscore
            )
            dna = self.extract_gene(seqs, hit.sseqid, hit.sstart, hit.send)
            ref = self.ref_prot_dict.get(hit.qseqid, "")
            
            # translation
            if hit.family == "23S":
                found = str(dna)
            elif hit.family in self.mutation_targets:
                found = self.best_translation(dna, ref)
            else:
                found = self.trim_to_ref(self.best_translation(dna, ref), ref)

            display_str = hit.family

            if "*" in found:
                if hit.family in self.mutation_targets:
                    continue

                pct = int((found.find("*") / len(ref)) * 100) if ref else 0
                trunc.append(f"{hit.family}-{pct}%(Stop)")
                continue

            is_strong = (hit.pident >= self.min_id) and (hit.coverage >= self.min_cov)
            
            if not is_strong:
                if hit.family in self.mutation_targets:
                    continue 
                
                if hit.coverage < 100: display_str += "?"
                if hit.pident < 100: display_str += "*"
                spur.append(display_str)
                continue

            if hit.family in self.mutation_targets and ref:
                mm = self.check_mutations(hit.family, found, ref)
                
                if mm:
                    mut_str = f"{hit.family} [{','.join(mm)}]"
                    muts.append(mut_str)
                    if hit.family in ["gyrA", "parC", "gyrB"]: cat_fq.append(mut_str)
                    elif hit.family == "23S": cat_oxa.append(mut_str)
                    else: cat_rif.append(mut_str)
                continue

            if hit.pident < 100:
                display_str += "^" if (ref and found == ref) else "*"
            
            if hit.coverage < 100:
                display_str += "?"

            strong.append(display_str)
            
            is_mec = hit.family in ["mecA", "mecC"]
            source = self.gene_source.get(hit.family, "other")

            if source == "mec_res" or is_mec: cat_mec.append(display_str)
            elif source == "bla": cat_bla.append(display_str)
            elif source == "amino": cat_amino.append(display_str)
            elif source == "oxa": cat_oxa.append(display_str)
            elif source == "tet": cat_tet.append(display_str)
            elif source == "gly": cat_gly.append(display_str)
            elif source == "mlsb": cat_mlsb.append(display_str)
            elif source == "rif": cat_rif.append(display_str)
            elif source == "flq": cat_fq.append(display_str)

        classes_map = {
            "Amino_res": cat_amino, "Bla_res": cat_bla, "Flq_res": cat_fq,
            "Gly_res": cat_gly, "Mec_res": cat_mec, "MLSB_res": cat_mlsb,
            "Oxa_res": cat_oxa, "Rif_res": cat_rif, "Tet_res": cat_tet,
        }
        
        active_classes = sum(1 for lst in classes_map.values() if len(lst) > 0)

        out["res_gene_count"] = str(len(strong)) # do not count mutations
        out["res_class_count"] = str(active_classes)
        out["truncated_resistance_hits"] = "; ".join(trunc) if trunc else "-"
        out["spurious_resistance_hits"] = "; ".join(spur) if spur else "-"
        out["Amino_res"] = "; ".join(cat_amino) if cat_amino else "-"
        out["Mec_res"] = "; ".join(cat_mec) if cat_mec else "-"
        out["Bla_res"] = "; ".join(cat_bla) if cat_bla else "-"
        out["Flq_res"] = "; ".join(cat_fq) if cat_fq else "-"
        out["Oxa_res"] = "; ".join(cat_oxa) if cat_oxa else "-"
        out["MLSB_res"] = "; ".join(cat_mlsb) if cat_mlsb else "-"
        out["Tet_res"] = "; ".join(cat_tet) if cat_tet else "-"
        out["Gly_res"] = "; ".join(cat_gly) if cat_gly else "-"
        out["Rif_res"] = "; ".join(cat_rif) if cat_rif else "-"        

        score = 0
        is_van_positive = any("vanA" in x or "vanB" in x for x in cat_gly)
        if is_van_positive: score = 3
        elif cat_mec: score = 2
        elif cat_bla: score = 1
            
        out["res_score"] = str(score)

        return out