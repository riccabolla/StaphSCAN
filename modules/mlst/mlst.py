import pandas as pd
import subprocess
import io
import sys
import re
from pathlib import Path
from Bio import SeqIO

class Module:
    def __init__(self):
        self.name = "mlst"
        self.module_dir = Path(__file__).parent
        self.data_dir = self.module_dir / "data"
        self.db_profiles = self.data_dir / "profiles.tsv"
        self.alleles_fasta = self.data_dir / "alleles.fasta"
        self.loci = ["arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"]     
        self.min_identity = 90.0
        self.min_coverage = 90.0 

    def check_db(self):
        if not self.db_profiles.exists(): return False

        if not self.alleles_fasta.exists():
            print("Building MLST Alleles FASTA (Reverse-BLAST Mode)...")
            combined_fasta = []
            
            for locus in self.loci:
                fpath = self.data_dir / f"{locus}.fas"
                if not fpath.exists():
                    print(f"Error: Missing {locus}.fas")
                    return False
                for record in SeqIO.parse(fpath, "fasta"):
                    raw_id = record.id
                    clean_num = raw_id.replace(f"{locus}_", "").replace(locus, "")
                    clean_num = clean_num.strip("-_")
                    if not clean_num: continue
                    
                    new_header = f">{locus}_{clean_num}"
                    combined_fasta.append(f"{new_header}\n{str(record.seq)}\n")
            
            try:
                with open(self.alleles_fasta, "w") as f:
                    f.write("".join(combined_fasta))
            except Exception as e:
                print(f"Error creating alleles FASTA: {e}")
                return False
                
        return True

    def run(self, assembly_path):
        cmd = [
            "blastn", "-task", "blastn", 
            "-query", str(self.alleles_fasta), 
            "-subject", str(assembly_path), 
            "-outfmt", "6 qseqid pident length qlen bitscore", 
            "-perc_identity", str(self.min_identity),
            "-qcov_hsp_perc", str(self.min_coverage),         
            "-dust", "no"
        ]
        
        result = {"ST": "Unknown"}
        for l in self.loci: result[l] = "-"

        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: 
                return result
            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["qseqid", "pident", "length", "qlen", "bitscore"])
            df['cov'] = (df['length'] / df['qlen']) * 100
            
            detected_profile = {} 
            
            for locus in self.loci:
                locus_hits = df[df['qseqid'].astype(str).str.startswith(f"{locus}_")].copy()
                
                if locus_hits.empty:
                    detected_profile[locus] = "-"
                    result[locus] = "-"
                    continue

                locus_hits = locus_hits.sort_values(['bitscore', 'pident', 'cov'], ascending=[False, False, False])
                best = locus_hits.iloc[0]
                qseqid = str(best['qseqid'])
                match = re.search(fr"{locus}_(\d+)", qseqid)
                allele_id = match.group(1) if match else "?"
                is_perfect = (best['pident'] >= 100.0) and (best['cov'] >= 100.0)                
                display_str = allele_id
                if not is_perfect: display_str += "*"
                
                detected_profile[locus] = allele_id
                result[locus] = display_str

            result["ST"] = self.resolve_st(detected_profile)
            return result

        except Exception as e:
            print(f"Error in MLST module: {e}", file=sys.stderr)
            result["ST"] = "Error"
            return result

    def resolve_st(self, observed_profile):
        try:
            dtype_map = {locus: 'object' for locus in self.loci}
            prof_df = pd.read_csv(self.db_profiles, sep="\t", dtype=dtype_map)
            
            best_st = "Undefined"
            min_mismatches = 7
            
            for _, row in prof_df.iterrows():
                mismatches = 0
                for locus in self.loci:
                    obs = observed_profile.get(locus, "-")
                    db = str(row[locus])
                    if obs != db: mismatches += 1
                
                if mismatches < min_mismatches:
                    min_mismatches = mismatches
                    best_st = str(row['ST'])
                if min_mismatches == 0: break
            
            if min_mismatches == 0: return f"ST{best_st}"
            elif min_mismatches == 1: return f"ST{best_st}-1LV"
            elif min_mismatches == 2: return f"ST{best_st}-2LV"
            return "Undefined"

        except Exception as e:
            return f"Error({e})"