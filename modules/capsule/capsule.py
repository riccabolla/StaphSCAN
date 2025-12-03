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
                return {"cap_type": "NT (Absent)", "cap_genes": "-"}

            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["qseqid", "pident", "length", "qlen", "bitscore"])
            
            def get_gene_base(qseqid):
                match = re.search(r"(cap[58][A-P])", str(qseqid))
                if match:
                    return match.group(1)
                return str(qseqid).split('_')[0] # Fallback

            df['gene_full'] = df['qseqid'].apply(get_gene_base)
            
            found_set = set(df['gene_full'].unique())
            genes_str = ";".join(sorted(list(found_set)))
            
            score_5 = 0
            score_8 = 0
            
            for locus in self.type_specific:
                if f"cap5{locus}" in found_set: score_5 += 1
                if f"cap8{locus}" in found_set: score_8 += 1
            
            target_type = None
            display_type = "NT"
            
            if score_5 > 0 and score_8 == 0:
                target_type = "cap5"
                display_type = "Type 5"
            elif score_8 > 0 and score_5 == 0:
                target_type = "cap8"
                display_type = "Type 8"
            elif score_5 > 0 and score_8 > 0:
                if score_5 > score_8: 
                    target_type = "cap5"
                    display_type = "Type 5 (Mixed)"
                elif score_8 > score_5: 
                    target_type = "cap8"
                    display_type = "Type 8 (Mixed)"
                else: 
                    display_type = "Conflict (5+8)"
            else:
                if genes_str: display_type = "NT (Only conserved genes found)"
                else: display_type = "-"

            if target_type:
                found_count = 0
                missing_genes = []
                
                for locus in self.operon_genes:
                    gene_name = f"{target_type}{locus}"
                    if gene_name in found_set:
                        found_count += 1
                    else:
                        missing_genes.append(locus)
                
                total_genes = len(self.operon_genes)
                
                if found_count == total_genes:
                    display_type += " (Complete)"
                else:
                    if len(missing_genes) > 4:
                        display_type += f" (Partial: {found_count}/{total_genes} genes)"
                    else:
                        display_type += f" (Missing: {','.join(missing_genes)})"

            return {
                "cap_type": display_type,
                "cap_genes": genes_str
            }

        except Exception as e:
            return {"cap_type": f"Error: {str(e)}", "cap_genes": "-"}