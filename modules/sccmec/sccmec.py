import pandas as pd
import subprocess
import io
from pathlib import Path

class Module:
    def __init__(self):
        self.name = "sccmec"
        self.module_dir = Path(__file__).parent
        self.db_fasta = self.module_dir / "data" / "targets.fasta"
        self.blast_db = self.module_dir / "data" / "sccmec_db"

    def check_db(self):
        if not self.db_fasta.exists(): return False
        if not (self.module_dir / "data" / "sccmec_db.nhr").exists():
            cmd = ["makeblastdb", "-in", str(self.db_fasta), "-dbtype", "nucl", "-out", str(self.blast_db)]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        return True

    def run(self, assembly_path):
        cmd = [
            "blastn", "-query", str(assembly_path), "-db", str(self.blast_db),
            "-outfmt", "6 sseqid pident length slen",
            "-max_target_seqs", "100"
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout:
                return {"sccmec_type": "Negative", "sccmec_genes": "-"}

            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["sseqid", "pident", "length", "slen"])
            
            df['cov'] = (df['length'] / df['slen']) * 100
            
            df = df[ (df['pident'] >= 90.0) & (df['cov'] >= 80.0) ]
            
            if df.empty:
                return {"sccmec_type": "Negative", "sccmec_genes": "-"}
            
            found_genes = set()
            for raw_id in df['sseqid'].unique():
                lower_id = raw_id.lower()
                
                # --- CCR Recombinases ---
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
                elif "ccrc" in lower_id: found_genes.add("ccrC1") # Fallback if C1 not specified
                
                # Mec Complex Components
                elif "meci" in lower_id: found_genes.add("mecI")
                elif "mecr1" in lower_id: found_genes.add("mecR1")
                elif "is1272" in lower_id: found_genes.add("IS1272")
                elif "is431" in lower_id: found_genes.add("IS431")
                
                # Specific Insertion Sequences
                elif "is12960d" in lower_id: found_genes.add("IS12960D") # Key for Type VII
                
                # Resistance Genes
                elif "meca" in lower_id: found_genes.add("mecA")
                elif "mecc" in lower_id: found_genes.add("mecC")
                elif "blaz" in lower_id: found_genes.add("blaZ") # Helper for Type XI

            mec_class = self.get_mec_class(found_genes)
            ccr_complexes = self.get_ccr_complex(found_genes)
            final_type = self.assign_type(mec_class, ccr_complexes, found_genes)
            genes_str = ";".join(sorted(list(found_genes)))

            return {
                "sccmec_type": final_type,
                "sccmec_genes": genes_str
            }

        except Exception as e:
            return {"sccmec_type": "Error", "sccmec_genes": str(e)}

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
        
        # Mixed Pairs for Higher Types
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
                if mec_class == "C" or mec_class == "B": 
                    types_found.append("Type X (A1B6-C)")
                elif mec_class == "A":
                    types_found.append("Type XV (A1B6-A)")

            # Type XI (ccrA1B3 based)
            elif ccr == "A1B3":
                # Usually Class A/E (mecA/mecC + mecI + mecR1)
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