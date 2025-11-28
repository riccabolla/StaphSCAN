import pandas as pd
import subprocess
import io
import sys
from pathlib import Path

class Module:
    def __init__(self):
        self.name = "mlst"
        self.module_dir = Path(__file__).parent
        self.data_dir = self.module_dir / "data"
        self.db_profiles = self.data_dir / "profiles.tsv"
        self.blast_db_prefix = self.data_dir / "mlst_db"
        self.loci = ["arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"]

    def check_db(self):
        if not self.db_profiles.exists(): return False
        if not (self.data_dir / "mlst_db.nhr").exists():
            missing = [l for l in self.loci if not (self.data_dir / f"{l}.fas").exists()]
            if missing: 
                print(f"Error: Missing MLST allele files: {missing}")
                return False
            
            print("Building MLST BLAST database...")
            combined = []
            for l in self.loci: 
                combined.append((self.data_dir/f"{l}.fas").read_text())
            cmd = [
                "makeblastdb", 
                "-in", "-", 
                "-dbtype", "nucl", 
                "-title", "StaphScan_MLST",
                "-out", str(self.blast_db_prefix)
            ]
            
            try:
                subprocess.run(cmd, input="\n".join(combined), text=True, check=True, stdout=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                print(f"Error creating BLAST DB: {e}")
                return False
                
        return True

    def run(self, assembly_path):
        cmd = [
            "blastn", "-query", str(assembly_path), "-db", str(self.blast_db_prefix),
            "-outfmt", "6 sseqid pident length slen bitscore", "-max_target_seqs", "2000"
        ]
        
        result = {"ST": "Unknown"}
        for l in self.loci: result[l] = "-"

        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: return result
            
            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["sseqid", "pident", "length", "slen", "bitscore"])
            df['cov'] = (df['length'] / df['slen']) * 100
            
            detected_profile = {} 
            
            for locus in self.loci:
                locus_hits = df[df['sseqid'].str.startswith(f"{locus}_")].copy()
                if locus_hits.empty:
                    detected_profile[locus] = "-"
                    result[locus] = "-"
                    continue
                
                locus_hits = locus_hits.sort_values('bitscore', ascending=False)
                best = locus_hits.iloc[0]
                
                parts = str(best['sseqid']).split('_')
                allele_id = parts[1] if len(parts) > 1 else "?"
                
                is_perfect = (best['pident'] == 100.0) and (best['cov'] >= 100.0)
                display_str = allele_id
                if not is_perfect: display_str += "*"
                
                detected_profile[locus] = allele_id
                result[locus] = display_str

            result["ST"] = self.resolve_st(detected_profile)
            return result

        except Exception as e:
            result["ST"] = "Error"
            return result

    def resolve_st(self, observed_profile):
        try:
            prof_df = pd.read_csv(self.db_profiles, sep="\t")
            cols = [c for c in prof_df.columns if c in self.loci]
            prof_df[cols] = prof_df[cols].astype(str)
            
            best_st = "Undefined"
            min_mismatches = 7
            
            for _, row in prof_df.iterrows():
                mismatches = 0
                for locus in self.loci:
                    if observed_profile.get(locus, "-") != row[locus]:
                        mismatches += 1
                
                if mismatches < min_mismatches:
                    min_mismatches = mismatches
                    best_st = str(row['ST'])
                if min_mismatches == 0: break
            
            if min_mismatches == 0: return f"ST{best_st}"
            elif min_mismatches == 1: return f"ST{best_st}-1LV"
            elif min_mismatches == 2: return f"ST{best_st}-2LV"
            return "Undefined"

        except: return "Error"