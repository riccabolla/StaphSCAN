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
        self.refs_fasta = self.data_dir / "refs.fasta" 
        self.loci = ["arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"]     
        self.min_identity = 95.0 
        self.min_coverage = 95.0       
        self.prof_df = None
        if self.db_profiles.exists():
            dtype_map = {locus: 'str' for locus in self.loci}
            dtype_map["ST"] = 'str'
            try:
                self.prof_df = pd.read_csv(self.db_profiles, sep="\t", dtype=dtype_map).fillna("-")
            except Exception as e:
                print(f"Warning: Could not load profiles.tsv: {e}", file=sys.stderr)

        self.allele_map = {l: {} for l in self.loci}
        if self.alleles_fasta.exists():
            self._load_allele_map()

    def _load_allele_map(self):
        """Loads all alleles into memory for exact string matching."""
        try:
            for record in SeqIO.parse(self.alleles_fasta, "fasta"):
                for locus in self.loci:
                    if record.id.startswith(f"{locus}_"):
                        seq_str = str(record.seq).upper()
                        clean_num = record.id.split('_')[-1]
                        self.allele_map[locus][seq_str] = clean_num
                        break
        except Exception as e:
            print(f"Error loading allele map: {e}", file=sys.stderr)

    def _get_closest_allele(self, locus, query_seq):
        """
        Fallback: If exact match fails, find the closest allele by SNP count.
        """
        best_id = "Novel"
        min_diffs = float('inf')

        for ref_seq, allele_id in self.allele_map[locus].items():
            if len(ref_seq) != len(query_seq):
                continue
                
            diffs = sum(1 for a, b in zip(ref_seq, query_seq) if a != b)
            
            if diffs < min_diffs:
                min_diffs = diffs
                best_id = allele_id
                
            if min_diffs == 1: 
                break
                
        if min_diffs == float('inf'):
            return "Novel*" 
            
        return f"{best_id}*"

    def check_db(self):
        if not self.db_profiles.exists(): return False

        if not self.alleles_fasta.exists() or not self.refs_fasta.exists():
            print("Building MLST Databases...")
            combined_fasta = []
            refs_fasta = []
            
            for locus in self.loci:
                fpath = self.data_dir / f"{locus}.fas"
                if not fpath.exists():
                    print(f"Error: Missing {locus}.fas")
                    return False
                
                first_record = True
                for record in SeqIO.parse(fpath, "fasta"):
                    clean_num = record.id.replace(f"{locus}_", "").replace(locus, "").strip("-_")
                    if not clean_num: continue
                    
                    header = f"{locus}_{clean_num}"
                    seq = str(record.seq).upper()
                    
                    combined_fasta.append(f">{header}\n{seq}\n")
                    
                    if first_record:
                        refs_fasta.append(f">{header}\n{seq}\n")
                        first_record = False
            
            try:
                with open(self.alleles_fasta, "w") as f:
                    f.write("".join(combined_fasta))
                with open(self.refs_fasta, "w") as f:
                    f.write("".join(refs_fasta))
                self._load_allele_map()
            except Exception as e:
                print(f"Error creating FASTA DBs: {e}")
                return False
                
        return True

    def run(self, assembly_path):
        cmd = [
            "blastn", "-task", "megablast", 
            "-query", str(self.refs_fasta), 
            "-subject", str(assembly_path), 
            "-outfmt", "6 qseqid sseq length qlen bitscore", 
            "-perc_identity", str(self.min_identity),
            "-qcov_hsp_perc", str(self.min_coverage),         
            "-dust", "no",
            "-max_target_seqs", "1" 
        ]
        
        result = {"ST": "Unknown"}
        for l in self.loci: result[l] = "-"

        try:
            if self.prof_df is None:
                result["ST"] = "DB_Error"
                return result

            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: 
                return result
            
            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["qseqid", "sseq", "length", "qlen", "bitscore"])
            
            detected_profile = {} 
            
            for locus in self.loci:
                mask = df['qseqid'].astype(str).str.startswith(f"{locus}_")
                locus_hits = df[mask]
                
                if locus_hits.empty:
                    detected_profile[locus] = "-"
                    result[locus] = "-"
                    continue

                best = locus_hits.sort_values('bitscore', ascending=False).iloc[0]
                
                genome_seq = best['sseq'].replace("-", "").upper()
                
                if best['length'] < best['qlen']:
                    result[locus] = "Partial"
                    detected_profile[locus] = "-"
                    continue

                if genome_seq in self.allele_map[locus]:
                    allele_id = self.allele_map[locus][genome_seq]
                    detected_profile[locus] = allele_id
                    result[locus] = allele_id
                else:
                    closest_str = self._get_closest_allele(locus, genome_seq)                  
                    detected_profile[locus] = "-" 
                    result[locus] = closest_str

            result["ST"] = self.resolve_st(detected_profile)
            return result

        except Exception as e:
            print(f"Error in MLST module: {e}", file=sys.stderr)
            result["ST"] = "Error"
            return result

    def resolve_st(self, observed_profile):
        """
        Vectorized ST resolution
        """
        if self.prof_df is None: return "DB_Error"
        
        try:
            obs_values = [str(observed_profile.get(locus, "-")) for locus in self.loci]
            
            matches_mask = (self.prof_df[self.loci] == obs_values)
            
            mismatch_counts = len(self.loci) - matches_mask.sum(axis=1)
            
            min_mismatches = mismatch_counts.min()
            
            if min_mismatches > 2:
                return "Undefined"
            
            best_match_idx = mismatch_counts.idxmin()
            best_st = str(self.prof_df.at[best_match_idx, 'ST'])
            
            if min_mismatches == 0: return f"ST{best_st}"
            elif min_mismatches == 1: return f"ST{best_st}-1LV"
            elif min_mismatches == 2: return f"ST{best_st}-2LV"
            return "Undefined"

        except Exception as e:
            return f"Error({e})"