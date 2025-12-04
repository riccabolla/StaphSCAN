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
        self.type_specific = ['H', 'I', 'J', 'K']
        self.operon_genes = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']

    def check_db(self):
        return self.db_fasta.exists()

    def run(self, assembly_path):
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
            def get_gene_base(qseqid):
                match = re.search(r"(cap[58][A-P])", str(qseqid))
                if match:
                    return match.group(1)
                return str(qseqid).split('_')[0] 

            df['gene_full'] = df['qseqid'].apply(get_gene_base)
            found_set = set(df['gene_full'].unique())
            genes_str = ";".join(sorted(list(found_set)))
            score_5 = 0
            score_8 = 0
            
            for locus in self.type_specific:
                if f"cap5{locus}" in found_set: score_5 += 1
                if f"cap8{locus}" in found_set: score_8 += 1
            
            detected_serotype = None
            display_type = "-"
            if score_5 > 0 and score_5 >= score_8:
                detected_serotype = "cap5"
                display_type = "Type 5"
            elif score_8 > 0 and score_8 > score_5:
                detected_serotype = "cap8"
                display_type = "Type 8"
            completeness = "-"
            
            if detected_serotype:
                expected_genes = [f"{detected_serotype}{locus}" for locus in self.operon_genes]                
                found_count = sum(1 for gene in expected_genes if gene in found_set)
                total_required = len(expected_genes)
                
                if found_count == total_required:
                    completeness = "Complete"
                else:
                    completeness = "Incomplete"

            return {
                "cap_type": display_type,
                "cap_completeness": completeness,
                "cap_genes": genes_str if genes_str else "-"
            }

        except Exception as e:
            return {"cap_type": "Error", "cap_completeness": "Error", "cap_genes": f"Error: {str(e)}"}