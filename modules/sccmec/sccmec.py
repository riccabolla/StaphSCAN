import pandas as pd
import subprocess
import io
from pathlib import Path

class Module:
    def __init__(self):
        self.name = "sccmec"
        self.module_dir = Path(__file__).parent
        self.db_fasta = self.module_dir / "data" / "targets.fasta"
        self.db_tsv = self.module_dir / "data" / "targets.tsv"
        self.blast_db = self.module_dir / "data" / "sccmec_db"

    def check_db(self):
        if not self.db_fasta.exists():
            return False
        if not (self.module_dir / "data" / "sccmec_db.nhr").exists():
            cmd = ["makeblastdb", "-in", str(self.db_fasta), "-dbtype", "nucl", "-out", str(self.blast_db)]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        return True

    def run(self, assembly_path):
        cmd = [
            "blastn", "-query", str(assembly_path), "-db", str(self.blast_db),
            "-outfmt", "6 qseqid sseqid pident length slen stitle bitscore",
            "-max_target_seqs", "20"
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout:
                return {"sccmec_type": "Negative", "sccmec_genes": "-"}
                
            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", 
                             names=["qseqid", "sseqid", "pident", "length", "slen", "stitle", "bitscore"])

            df['cov'] = (df['length'] / df['slen']) * 100
            df = df[ (df['pident'] >= 90.0) & (df['cov'] >= 90.0) ]
            
            if df.empty:
                return {"sccmec_type": "Negative", "sccmec_genes": "-"}

            df = df.sort_values('bitscore', ascending=False)
            
            meta_df = pd.read_csv(self.db_tsv, sep="\t")

            def extract_accession(val):
                try:
                    parts = val.split(' ')
                    if len(parts) > 1:
                        raw_acc = parts[1]
                        if '|' in raw_acc: return raw_acc.split('|')[0]
                        return raw_acc
                    return val
                except: return val

            df['accession_key'] = df['stitle'].apply(extract_accession)
            merged = pd.merge(df, meta_df, left_on='accession_key', right_on='accession', how='left')
            merged['gene_family'] = merged['target'].astype(str).str[:4] 
            ccr_subset = merged[merged['target'].str.contains("ccr", case=False, na=False)].copy()
            other_subset = merged[~merged['target'].str.contains("ccr", case=False, na=False)].copy()
            ccr_subset = ccr_subset.sort_values('pident', ascending=False).drop_duplicates('gene_family')
            final_hits = pd.concat([ccr_subset, other_subset])
            all_genes_found = sorted(final_hits['target'].dropna().unique())
            genes_str = ";".join(all_genes_found)
            report_hits = final_hits[final_hits['target'].str.contains("ccr", case=False, na=False)]

            if report_hits.empty:
                 res_genes = [g for g in all_genes_found if "mec" in g or "bla" in g]
                 if res_genes:
                     return {"sccmec_type": "Orphan mecA (No ccr)", "sccmec_genes": genes_str}
                 else:
                     return {"sccmec_type": "Negative", "sccmec_genes": "-"}

            type_counts = report_hits['type'].value_counts()
            top_type = type_counts.idxmax()
            
            final_call = f"Type {top_type}"
            
            if len(type_counts) > 1:
                final_call += f" (Composite/Conflict: {list(type_counts.index)})"

            return {
                "sccmec_type": final_call,
                "sccmec_genes": genes_str
            }

        except Exception as e:
            return {"sccmec_type": "Error", "sccmec_genes": str(e)}