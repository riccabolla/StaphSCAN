import pandas as pd
import subprocess
import io
from pathlib import Path

class Module:
    def __init__(self):
        self.name = "sccmec"
        self.module_dir = Path(__file__).parent
        self.data_dir = self.module_dir / "data"
        db_targets = list(self.data_dir.glob("targets*.fasta"))
        db_regions = list(self.data_dir.glob("regions*.fasta"))
        if not db_targets:
            raise FileNotFoundError(f"Missing database: No targets fasta file found in {self.data_dir}")
        if not db_regions:
            raise FileNotFoundError(f"Missing database: No regions fasta file found in {self.data_dir}")
        self.targets_fasta = db_targets[0]
        self.targets_db = self.data_dir / "sccmec_targets_db"
        self.regions_fasta = db_regions[0]
        self.regions_db = self.data_dir / "sccmec_regions_db"

    def check_db(self):
        
        if not self.targets_fasta.exists(): return False
        if not (self.data_dir / "sccmec_targets_db.nhr").exists():
            return False
            
        if not self.regions_fasta.exists(): return False
        if not (self.data_dir / "sccmec_regions_db.nhr").exists():
            return False
            
        return True
            
        return True

    def run(self, assembly_path):
        import re #added for type string parsing
        # Determine type
        target_res = self.run_blast(assembly_path, self.targets_db, min_ident=90.0, min_cov=80.0)
        
        sccmec_type = "Negative"
        sccmec_genes = "-"
        
        if not target_res.empty:
            found_genes = set()
            for raw_id in target_res['sseqid'].unique():
                lower_id = raw_id.lower()                
                if "ccra1" in lower_id: found_genes.add("ccrA1")
                elif "ccrb1" in lower_id: found_genes.add("ccrB1")
                elif "ccra2" in lower_id: found_genes.add("ccrA2")
                elif "ccrb2" in lower_id: found_genes.add("ccrB2")
                elif "ccra3" in lower_id: found_genes.add("ccrA3")
                elif "ccrb3" in lower_id: found_genes.add("ccrB3")
                elif "ccra4" in lower_id: found_genes.add("ccrA4")
                elif "ccrb4" in lower_id: found_genes.add("ccrB4")
                elif "ccrb6" in lower_id: found_genes.add("ccrB6")                
                # Differentiate ccrC1 vs ccrC2
                elif "ccrc1" in lower_id: found_genes.add("ccrC1")
                elif "ccrc2" in lower_id: found_genes.add("ccrC2")
                elif "ccrc" in lower_id: found_genes.add("ccrC1")                
                # Mec Complex Components
                elif "meci" in lower_id: found_genes.add("mecI")
                elif "mecr1" in lower_id: found_genes.add("mecR1")
                elif "is1272" in lower_id: found_genes.add("IS1272")
                elif "is431" in lower_id: found_genes.add("IS431")                
                # Specific Insertion Sequences
                elif "is12960d" in lower_id: found_genes.add("IS12960D")                
                # Resistance Genes
                elif "meca" in lower_id: found_genes.add("mecA")
                elif "mecc" in lower_id: found_genes.add("mecC")
                elif "blaz" in lower_id: found_genes.add("blaZ")

            mec_class = self.get_mec_class(found_genes)
            ccr_complexes = self.get_ccr_complex(found_genes)
            sccmec_type = self.assign_type(mec_class, ccr_complexes, found_genes)
            sccmec_genes = ";".join(sorted(list(found_genes)))

        base_type = None
        match = re.search(r"Type\s+([IVX]+)", sccmec_type)
        if match:
            base_type = match.group(1)    
        # Determine subtype
        sccmec_subtype = "-"
        if self.regions_fasta.exists():
            region_res = self.run_blast(assembly_path, self.regions_db, min_ident=90.0, min_cov=70.0)
            if not region_res.empty:
                region_res['subtype_name'] = region_res['sseqid'].apply(
                    lambda x: str(x).split("|")[-1] if "|" in str(x) else str(x)
                )
                if base_type:
                    pattern = f"^{base_type}[a-z]*$"
                    filtered_res = region_res[region_res['subtype_name'].str.match(pattern, na=False)]
                    if not filtered_res.empty:
                        region_res = filtered_res
                    else:
                        region_res = pd.DataFrame()            
            if not region_res.empty:
                region_res = region_res.sort_values(by=['cov', 'pident'], ascending=[False, False])
                sccmec_subtype = region_res.iloc[0]['subtype_name']

        return {
            "sccmec_type": sccmec_type,
            "sccmec_subtype": sccmec_subtype,
            "sccmec_genes": sccmec_genes
        }

    def run_blast(self, assembly_path, db_path, min_ident, min_cov):
        """Helper to run blast and calculate cumulative sequence coverage"""
        #added sstart and send to map fragmented contigs
        cmd = [
            "blastn", "-query", str(assembly_path), "-db", str(db_path),
            "-outfmt", "6 sseqid pident length slen sstart send",
            "-max_target_seqs", "100"
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout or not res.stdout.strip():
                return pd.DataFrame()

            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["sseqid", "pident", "length", "slen", "sstart", "send"])
            
            if df.empty:
                return pd.DataFrame()

            df = df[df['pident'] >= min_ident]
            if df.empty:
                return pd.DataFrame()
            
            records = []
            
            # it does not look to single blast-hit, instead it looks to cumulative coverage of all blast-hits 
            # to the same reference sequence, which is more robust to fragmented assemblies
            # added to fix subtyiping issues
            for sseqid, group in df.groupby('sseqid'):
                slen = group['slen'].iloc[0]                
                intervals = []
                for _, row in group.iterrows():
                    start, end = sorted([row['sstart'], row['send']])
                    intervals.append([start, end])
                
                intervals.sort()
                merged = []
                for interval in intervals:
                    if not merged:
                        merged.append(interval)
                    else:
                        prev = merged[-1]
                        if interval[0] <= prev[1] + 1:
                            prev[1] = max(prev[1], interval[1])
                        else:
                            merged.append(interval)
                
                # Calculate total cumulative coverage length and percentage
                cov_len = sum(interval[1] - interval[0] + 1 for interval in merged)
                cov = (cov_len / slen) * 100
                
                # Calculate weighted identity across all alignments
                weighted_ident = (group['pident'] * group['length']).sum() / group['length'].sum()
                
                records.append({
                    "sseqid": sseqid,
                    "pident": weighted_ident,
                    "length": cov_len,
                    "slen": slen,
                    "cov": cov
                })
                
            res_df = pd.DataFrame(records)
            if res_df.empty:
                return res_df
            
            res_df = res_df[res_df['cov'] >= min_cov]
            
            res_df = res_df.sort_values(by=['cov', 'pident'], ascending=[False, False])
            
            return res_df
        except Exception:
            return pd.DataFrame()

    def get_mec_class(self, genes):
        has_mec = "mecA" in genes or "mecC" in genes
        if not has_mec: return "None"
        
        # Class A: mecI + mecR1
        if "mecI" in genes and "mecR1" in genes:
            return "A"
        # Class B: IS1272 + mecA
        elif "IS1272" in genes:
            return "B"
        # Class C: IS431 + mecA (and no mecI/IS1272)
        elif "IS431" in genes and "mecI" not in genes and "IS1272" not in genes:
            return "C"
        # Class E: Contains mecC
        elif "mecC" in genes:
            return "E"

        # Fallbacks
        if "mecA" in genes: return "Unknown (mecA+)"
        if "mecC" in genes: return "Unknown (mecC+)"
        
        return "None"

    def get_ccr_complex(self, genes):
        complexes = []
        
        # Standard Pairs
        if "ccrA1" in genes and "ccrB1" in genes: complexes.append("1")
        if "ccrA2" in genes and "ccrB2" in genes: complexes.append("2")
        if "ccrA3" in genes and "ccrB3" in genes: complexes.append("3")
        if "ccrA4" in genes and "ccrB4" in genes: complexes.append("4")
        
        # C-type Complexes
        if "ccrC1" in genes: complexes.append("C1") 
        if "ccrC2" in genes: complexes.append("C2")
        
        # Mixed Pairs
        if "ccrA1" in genes and "ccrB6" in genes: complexes.append("A1B6") 
        if "ccrA1" in genes and "ccrB3" in genes: complexes.append("A1B3") 

        return complexes

    def assign_type(self, mec_class, ccr_list, all_genes):
        if mec_class == "None": return "Negative"
        if not ccr_list: return f"Orphan {mec_class} (No ccr)"
        
        types_found = []
        
        for ccr in ccr_list:
            # Standard Types (I - IV, VI, VIII)
            if ccr == "1" and mec_class == "B": types_found.append("Type I(1B)")
            elif ccr == "2" and mec_class == "A": types_found.append("Type II(2A)")
            elif ccr == "3" and mec_class == "A": types_found.append("Type III(3A)")
            elif ccr == "2" and mec_class == "B": types_found.append("Type IV(2B)")
            elif ccr == "4" and mec_class == "B": types_found.append("Type VI(4B)")
            elif ccr == "4" and mec_class == "A": types_found.append("Type VIII(4A)")
            
            # Type V / VII (ccrC1 based)
            elif ccr == "C1" and mec_class == "C":
                # Check for IS12960D to distinguish V from VII
                if "IS12960D" in all_genes:
                    types_found.append("Type VII (5C + IS12960D)")
                else:
                    types_found.append("Type V (5C)")
                    
            # Type XIV (ccrC1 based)
            elif ccr == "C1" and mec_class == "A":
                types_found.append("Type XIV (5A)")

            # Type IX (ccr1 + Class C)
            elif ccr == "1" and mec_class == "C":
                types_found.append("Type IX (1C)")

            # Type X / XV (ccrA1B6 based)
            elif ccr == "A1B6":
                if mec_class in ["C", "B"]: 
                    types_found.append("Type X (A1B6-C)")
                elif mec_class == "A":
                    types_found.append("Type XV (A1B6-A)")

            # Type XI (ccrA1B3 based)
            elif ccr == "A1B3":
                # Usually Class E/A
                types_found.append("Type XI (A1B3)")

            # Type XII / XIII (ccrC2 based)
            elif ccr == "C2" and mec_class == "C":
                types_found.append("Type XII (9C)")
            elif ccr == "C2" and mec_class == "A":
                types_found.append("Type XIII (9A)")
            
            else:
                types_found.append(f"Type ? ({ccr}-{mec_class})")

        if len(types_found) == 1:
            return types_found[0]
        else:
            return f"Composite ({' + '.join(types_found)})"