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
            print("Warning: Mash sketch 'staph_refs.msh' not found. Run setup_mash.py.")
            return False
        return True

    def run(self, assembly_path):
        contig_count, n50, longest, total_size, ambig = self.get_contig_stats(assembly_path)
        species_call, mash_dist = self.run_mash(assembly_path)
        qc_failures = []
        if total_size < self.min_size: qc_failures.append("Undersized")
        if total_size > self.max_size: qc_failures.append("Oversized")
        if n50 < self.min_n50: qc_failures.append("Low_N50")
        if "yes" in ambig: qc_failures.append("Ambiguous_Bases")
        
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
            if not res.stdout: return "Unknown", "-"
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
            display_name = "Unknown"
            if "S_aureus" in best_name: display_name = "S. aureus"
            elif "S_epidermidis" in best_name: display_name = "S. epidermidis"
            elif "S_lugdunensis" in best_name: display_name = "S. lugdunensis"
            elif "S_haemolyticus" in best_name: display_name = "S. haemolyticus"
            else: display_name = best_name
            if best_dist <= 0.02:
                if display_name == "S. aureus":
                    return "S. aureus (Strong match)", str(best_dist)
                else:
                    return display_name, str(best_dist)
            elif best_dist <= 0.04:
                if display_name == "S. aureus":
                    return "S. aureus (Weak match)", str(best_dist)
                else:
                    return display_name, str(best_dist)
            else:
                return "No match found", str(best_dist)

        except Exception as e:
            return f"Error ({e})", "-"

    def get_contig_stats(self, fasta_path):
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