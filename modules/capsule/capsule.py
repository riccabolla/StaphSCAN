import pandas as pd
import subprocess
import io
import sys
import re
from pathlib import Path

class Module:
    def __init__(self):
        self.name = "capsule"
        self.module_dir = Path(__file__).parent
        self.db_fasta = self.module_dir / "data" / "capsule_targets.fasta"
        
        # Genes specific to serotype determination (The central region)
        self.type_specific = ['H', 'I', 'J', 'K']
        
        # The full operon required for "Completeness"
        self.operon_genes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']

    def check_db(self):
        return self.db_fasta.exists()

    def run(self, assembly_path):
        # Default return values
        defaults = {
            "cap_type": "-", 
            "cap_completeness": "-", 
            "cap_genes": "-"
        }

        cmd = [
            "blastn", "-task", "blastn", 
            "-query", str(self.db_fasta), 
            "-subject", str(assembly_path),
            "-outfmt", "6 qseqid pident length qlen bitscore",
            "-perc_identity", "90",  
            "-qcov_hsp_perc", "80"   
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout:
                return defaults

            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["qseqid", "pident", "length", "qlen", "bitscore"])
            
            # Helper to extract gene name (e.g., cap5H from header)
            def get_gene_base(qseqid):
                match = re.search(r"(cap[58][A-P])", str(qseqid))
                if match:
                    return match.group(1)
                return str(qseqid).split('_')[0] 

            df['gene_full'] = df['qseqid'].apply(get_gene_base)
            
            # Create a set of all found genes
            found_set = set(df['gene_full'].unique())
            genes_str = ";".join(sorted(list(found_set)))
            
            # --- 1. Determine Serotype (Based on Specific Genes H-K) ---
            score_5 = 0
            score_8 = 0
            
            for locus in self.type_specific:
                if f"cap5{locus}" in found_set: score_5 += 1
                if f"cap8{locus}" in found_set: score_8 += 1
            
            detected_serotype = None
            display_type = "-"
            
            # Logic: If specific genes are present, assign type
            if score_5 > 0 and score_5 >= score_8:
                detected_serotype = "cap5"
                display_type = "Type 5"
            elif score_8 > 0 and score_8 > score_5:
                detected_serotype = "cap8"
                display_type = "Type 8"
            
            # --- 2. Determine Completeness (Based on full operon) ---
            completeness = "-"
            
            if detected_serotype:
                # Generate list of expected genes for this serotype (e.g., cap5A, cap5B...)
                expected_genes = [f"{detected_serotype}{locus}" for locus in self.operon_genes]
                
                # Count how many expected genes are actually found
                found_count = sum(1 for gene in expected_genes if gene in found_set)
                total_required = len(expected_genes)
                
                if found_count == total_required:
                    completeness = "Complete"
                else:
                    # Optional: You could add details like "Incomplete (12/16)" if desired
                    completeness = "Incomplete"

            return {
                "cap_type": display_type,
                "cap_completeness": completeness,
                "cap_genes": genes_str if genes_str else "-"
            }

        except Exception as e:
            return {"cap_type": "Error", "cap_completeness": "Error", "cap_genes": f"Error: {str(e)}"}