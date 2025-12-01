import pandas as pd
import subprocess
import io
import re
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
                if "ccra1" in lower_id: found_genes.add("ccrA1")
                elif "ccrb1" in lower_id: found_genes.add("ccrB1")
                elif "ccra2" in lower_id: found_genes.add("ccrA2")
                elif "ccrb2" in lower_id: found_genes.add("ccrB2")
                elif "ccra3" in lower_id: found_genes.add("ccrA3")
                elif "ccrb3" in lower_id: found_genes.add("ccrB3")
                elif "ccra4" in lower_id: found_genes.add("ccrA4")
                elif "ccrb4" in lower_id: found_genes.add("ccrB4")
                elif "ccrc" in lower_id: found_genes.add("ccrC")
                elif "meci" in lower_id: found_genes.add("mecI")
                elif "mecr1" in lower_id: found_genes.add("mecR1")
                elif "is1272" in lower_id: found_genes.add("IS1272")
                elif "is431" in lower_id: found_genes.add("IS431")
                elif "meca" in lower_id: found_genes.add("mecA")
                elif "mecc" in lower_id: found_genes.add("mecC")

            mec_class = self.get_mec_class(found_genes)

            ccr_complexes = self.get_ccr_complex(found_genes)

            final_type = self.assign_type(mec_class, ccr_complexes)

            genes_str = ";".join(sorted(list(found_genes)))

            return {
                "sccmec_type": final_type,
                "sccmec_genes": genes_str
            }

        except Exception as e:
            return {"sccmec_type": "Error", "sccmec_genes": str(e)}

    def get_mec_class(self, genes):
        """
        Class A: mecI + mecR1 + mecA
        Class B: IS1272 + mecA (mecR1 truncated, no mecI)
        Class C: IS431 + mecA (mecR1 truncated, no mecI, no IS1272)
        """
        has_mec = "mecA" in genes or "mecC" in genes
        if not has_mec: return "None"
        
        if "mecI" in genes and "mecR1" in genes:
            return "A"
        elif "IS1272" in genes:
            return "B"
        elif "IS431" in genes and "mecI" not in genes and "IS1272" not in genes:
            return "C"

        if "mecA" in genes: return "Unknown (mecA+)"
        return "None"

    def get_ccr_complex(self, genes):
        """
        Returns a list of detected ccr complexes
        """
        complexes = []
        
        # Type 1: ccrA1 + ccrB1
        if "ccrA1" in genes and "ccrB1" in genes: complexes.append("1")
        
        # Type 2: ccrA2 + ccrB2
        if "ccrA2" in genes and "ccrB2" in genes: complexes.append("2")
        
        # Type 3: ccrA3 + ccrB3
        if "ccrA3" in genes and "ccrB3" in genes: complexes.append("3")
        
        # Type 4: ccrA4 + ccrB4
        if "ccrA4" in genes and "ccrB4" in genes: complexes.append("4")
        
        # Type 5: ccrC
        if "ccrC" in genes: complexes.append("5")
        
        return complexes

    def assign_type(self, mec_class, ccr_list):
        if mec_class == "None": return "Negative"
        if not ccr_list: return f"Orphan {mec_class} (No ccr)"
        
        types_found = []
        
        for ccr in ccr_list:
            
            if ccr == "1" and mec_class == "B": types_found.append("Type I(1B)")
            
            elif ccr == "2" and mec_class == "A": types_found.append("Type II(2A)")
            elif ccr == "2" and mec_class == "B": types_found.append("Type IV(2B)")
            
            elif ccr == "3" and mec_class == "A": types_found.append("Type III(3A)")
            
            elif ccr == "4" and mec_class == "B": types_found.append("Type VI(4B)")
            elif ccr == "4" and mec_class == "A": types_found.append("Type VIII(4A)")
            
            elif ccr == "5" and mec_class == "C": types_found.append("Type V(5C)")
            elif ccr == "5" and mec_class == "C2": types_found.append("Type V(5C2)")
            # Broad Type V catch
            elif ccr == "5" and "Unknown" in mec_class: types_found.append("Type V(5C)")
            elif ccr == "5" and mec_class == "B": types_found.append("Type V(5B) variant")

            else:
                types_found.append(f"Type ? ({ccr}{mec_class})")

        if len(types_found) == 1:
            return types_found[0]
        else:
            return f"Composite ({' + '.join(types_found)})"