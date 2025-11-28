import pandas as pd
import subprocess
import io
from pathlib import Path
from Bio.Seq import Seq

class Module:
    def __init__(self):
        self.name = "virulence"
        self.module_dir = Path(__file__).parent
        self.db_fasta = self.module_dir / "data" / "targets.fasta"
        self.blast_db = self.module_dir / "data" / "vir_db"

    def check_db(self):
        if not self.db_fasta.exists(): return False
        if not (self.module_dir / "data" / "vir_db.nhr").exists():
            cmd = ["makeblastdb", "-in", str(self.db_fasta), "-dbtype", "nucl", "-out", str(self.blast_db)]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        return True

    def translate_dna(self, dna_seq, sstart, strand):
        """
        Translates DNA based on alignment frame and strand.
        Returns the protein sequence string (stopped at first *) and its length.
        """
        seq_obj = Seq(dna_seq)
        if strand == "minus":
            seq_obj = seq_obj.reverse_complement()
        
        frame = (sstart - 1) % 3
        trim_amt = (3 - frame) % 3
        if trim_amt == 3: trim_amt = 0
        
        clean_seq = str(seq_obj)[trim_amt:]
        clean_seq = clean_seq[:len(clean_seq)-(len(clean_seq)%3)]
        
        if not clean_seq: return "", 0

        prot_seq = str(Seq(clean_seq).translate())
        
        if "*" in prot_seq:
            stop_index = prot_seq.find("*")
            prot_seq = prot_seq[:stop_index]
            
        return prot_seq, len(prot_seq)

    def run(self, assembly_path):
        cmd = [
            "blastn", "-query", str(assembly_path), "-db", str(self.blast_db),
            "-outfmt", "6 sseqid pident length slen qseq sstart sstrand",
            "-max_target_seqs", "100"
        ]
        
        results = {"vir_pvl": "-", "vir_tsst": "-", "vir_et": "-", "vir_spurious": "-"}

        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: return results

            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", 
                             names=["sseqid", "pident", "length", "slen", "qseq", "sstart", "sstrand"])
            
            df['cov'] = (df['length'] / df['slen']) * 100
            
            df = df[(df['pident'] >= 80.0) & (df['cov'] >= 40.0)]
            
            strong_hits = []
            spurious_hits = []

            for gene in df['sseqid'].unique():
                gene_rows = df[df['sseqid'] == gene]
                best = gene_rows.sort_values('pident', ascending=False).iloc[0]
                
                pid = best['pident']
                cov = best['cov']
                
                is_weak_id = (pid >= 80.0 and pid < 90.0)
                is_weak_cov = (cov >= 40.0 and cov < 80.0)
                
                if is_weak_id or is_weak_cov:
                    spurious_hits.append(f"{gene} (Id:{int(pid)}% Cov:{int(cov)}%)")
                    continue 

                annotation = ""
                
                if pid < 100.0 or cov < 100.0:
                    annotation += "*"
                
                prot_seq, obs_aa_len = self.translate_dna(best['qseq'], int(best['sstart']), best['sstrand'])
                
                ref_aa_len = best['slen'] / 3
                ratio = (obs_aa_len / ref_aa_len) * 100
                
                if ratio < 90.0:
                    annotation += f"-{int(ratio)}%"
                
                strong_hits.append(f"{gene}{annotation}")

            def get_val(lookup):
                found = [x for x in strong_hits if lookup in x.lower()]
                return ", ".join(found) if found else "-"

            s = get_val("luks")
            f = get_val("lukf")
            
            if s != "-" and f != "-": 
                results["vir_pvl"] = "Positive"
            elif s != "-" or f != "-": 
                results["vir_pvl"] = f"Partial ({s if s!='-' else f})"

            results["vir_tsst"] = get_val("tst")
            
            et_list = [x for x in strong_hits if "eta" in x.lower() or "etb" in x.lower()]
            results["vir_et"] = ", ".join(et_list) if et_list else "-"
            
            results["vir_spurious"] = "; ".join(spurious_hits) if spurious_hits else "-"

            return results

        except Exception as e:
            return {"vir_pvl": "Error", "vir_tsst": str(e), "vir_et": "Error"}