import pandas as pd
import subprocess
import io
from pathlib import Path

class Module:
    def __init__(self):
        self.name = "agr"
        self.module_dir = Path(__file__).parent
        self.db_fasta = self.module_dir / "data" / "targets.fasta"
        self.blast_db = self.module_dir / "data" / "agr_db"

    def check_db(self):
        if not self.db_fasta.exists():
            return False
        if not (self.module_dir / "data" / "agr_db.nhr").exists():
            cmd = ["makeblastdb", "-in", str(self.db_fasta), "-dbtype", "nucl", "-out", str(self.blast_db), "-parse_seqids"]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        return True

    def run(self, assembly_path):
        cmd = [
            "blastn", 
            "-task", "blastn-short",
            "-query", str(assembly_path), 
            "-db", str(self.blast_db),
            "-outfmt", "6 sseqid pident length slen bitscore",
            "-perc_identity", "90",  
            "-max_target_seqs", "5000" 
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout:
                return {"agr_type": "Negative", "agr_match_confidence": "0%"}
            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["sseqid", "pident", "length", "slen", "bitscore"])
            df['coverage'] = (df['length'] / df['slen']) * 100 #coverage percentage
            df = df.sort_values(by=['pident', 'coverage', 'bitscore'], ascending=[False, False, False]) # sorting
            best_hit = df.iloc[0]
            raw_id = str(best_hit['sseqid'])
            if '_' in raw_id:
                group_id = raw_id.split('_')[0]
            else:
                group_id = raw_id
            type_map = {'gp1': 'I', 'gp2': 'II', 'gp3': 'III', 'gp4': 'IV'}
            final_type = type_map.get(group_id, group_id) 
            confidence = f"{best_hit['pident']}%"
                    
            return {
                "agr_type": f"agr {final_type}",
                "agr_match_confidence": confidence
            }

        except Exception as e:
            return {"agr_type": "Error", "agr_match_confidence": "0%"}