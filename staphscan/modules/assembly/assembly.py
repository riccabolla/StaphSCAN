import collections
import subprocess
import pandas as pd
import io
from pathlib import Path

class Module:
    def __init__(self):
        self.name = "assembly"
        self.module_dir = Path(__file__).parent
        self.db_sketch = self.module_dir / "data" / "staph_refs.msh"      
        self.min_size = 2600000 
        self.max_size = 3100000
        self.min_n50 = 10000

    def check_db(self):
        if not self.db_sketch.exists():
            print(f"Warning: Mash sketch '{self.db_sketch.name}' not found. Skipping assembly module.")
            return False
        return True

    def run(self, assembly_path):
        # if db is missing skip and report na results
        if not self.check_db():
            return {
                "Species": "Not Analyzed",
                "Mash_distance": "-",
                "Total_size": "-",
                "QC": "SKIPPED (Missing DB)",
                "contig_count": "-",
                "N50": "-",
                "largest_contig": "-",
                "ambiguous_bases": "-"
            }
        contig_count, n50, longest, total_size, ambig = self.get_contig_stats(assembly_path)
        
        species_call, mash_dist, is_mixed = self.run_mash(assembly_path)
        
        qc_failures = []
        if total_size < self.min_size: qc_failures.append("Too short")
        if total_size > self.max_size: qc_failures.append("Too long")
        if n50 < self.min_n50: qc_failures.append("Low_N50")
        if "yes" in ambig: qc_failures.append("Ambiguous_Bases")
        
        # Add contamination QC checks
        if is_mixed: qc_failures.append("Mixed")
        #if is_divergent: qc_failures.append("High_Mash_Distance")
        
        qc_status = "PASS"
        if qc_failures:
            qc_status = f"FAILED ({','.join(qc_failures)})"

        return {
            "Species": species_call,
            "Mash_distance": mash_dist,
            "Total_size": total_size,
            "QC": qc_status,
            "contig_count": contig_count,
            "N50": n50,
            "largest_contig": longest,
            "ambiguous_bases": ambig
        }

    def run_mash(self, assembly_path):
        try:
            cmd = ["mash", "dist", str(self.db_sketch), str(assembly_path)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: return "Unknown", "-", False, False
            
            rows = []
            for line in res.stdout.strip().split('\n'):
                parts = line.split('\t')
                ref_name = Path(parts[0]).name 
                dist = float(parts[2])
                rows.append((ref_name, dist))
                
            rows.sort(key=lambda x: x[1])
            best_match = rows[0]
            best_name = best_match[0]
            best_dist = best_match[1]

            # Contamination Logic: Check if multiple species are under 0.05 distance
            close_matches = [r for r in rows if r[1] <= 0.05]
            is_mixed = len(close_matches) > 1
            
            # Divergence Logic: Best match is acceptable, but unusually high for a pure genome
            #is_divergent = 0.02 < best_dist <= 0.04

            # Cleaner dictionary approach for species mapping
            species_map = {
                "S_aureus": "S. aureus",
                "S_epidermidis": "S. epidermidis",
                "S_lugdunensis": "S. lugdunensis",
                "S_haemolyticus": "S. haemolyticus",
                "S_argenteus": "S. argenteus",
                "S_capitis": "S. capitis",
                "S_schweitzeri": "S. schweitzeri"
            }
            
            display_name = next((v for k, v in species_map.items() if k in best_name), best_name)

            if best_dist <= 0.04:
                return display_name, str(best_dist), is_mixed
            else:
                return "No match found", str(best_dist), is_mixed

        except Exception as e:
            return f"Error ({e})", "-", False, False

    def get_contig_stats(self, fasta_path):
        # [Keep your existing get_contig_stats code here - no changes needed]
        lengths = []
        ambiguous_count = 0
        with open(fasta_path, 'r') as f:
            seq_buffer = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq_buffer:
                        s = "".join(seq_buffer)
                        lengths.append(len(s))
                        ambiguous_count += sum(1 for b in s if b.upper() not in "ATCG")
                    seq_buffer = []
                else:
                    seq_buffer.append(line)
            if seq_buffer:
                s = "".join(seq_buffer)
                lengths.append(len(s))
                ambiguous_count += sum(1 for b in s if b.upper() not in "ATCG")

        if not lengths: return 0, 0, 0, 0, "no"

        lengths.sort()
        longest = lengths[-1]
        total = sum(lengths)
        half = total / 2
        cum_sum = 0
        n50 = 0
        for l in reversed(lengths):
            cum_sum += l
            if cum_sum >= half:
                n50 = l
                break
        
        ambig_str = f"yes ({ambiguous_count})" if ambiguous_count > 0 else "no"
        return len(lengths), n50, longest, total, ambig_str