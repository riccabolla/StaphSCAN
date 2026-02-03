import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO

class Module:
    def __init__(self):
        self.name = "spa"
        self.module_dir = Path(__file__).parent
        self.data_dir = self.module_dir / "data"
        self.repeats_file = self.data_dir / "repeats.fasta"
        self.types_file = self.data_dir / "spatypes.csv"

    def check_db(self):
        return self.repeats_file.exists() and self.types_file.exists()

    def load_repeats(self):
        """ 
        Load repeats into a dict: { 'SEQUENCE': '01' } 
        Strips 'r' prefix if present.
        """
        repeats = {}
        current_id = ""
        if not self.repeats_file.exists(): return {}

        with open(self.repeats_file, "r") as f:
            lines = f.read().splitlines()
        
        for line in lines:
            if line.startswith(">"):
                raw_id = line[1:]
                if raw_id.lower().startswith('r'):
                    current_id = raw_id[1:]
                else:
                    current_id = raw_id
            else:
                if current_id: repeats[line.upper()] = current_id
        return repeats

    def load_types(self):
        types = {}
        if not self.types_file.exists(): return {}

        with open(self.types_file, "r") as f:
            for line in f:
                if ";" in line: parts = line.strip().split(';')
                else: parts = line.strip().split(',')
                
                if len(parts) >= 2:
                    types[parts[1]] = parts[0]
        return types

    def run(self, assembly_path):
        try:
            sequences = {}
            for record in SeqIO.parse(assembly_path, "fasta"):
                sequences[record.id] = str(record.seq).upper()

            primer_sets = [
                ('TAAAGACGATCCTTCGGTGAG', 'CAGCAGTAGTGCCGTTTGCTT'),        # Set 1
                ('AGACGATCCTTCGGTGAGC', 'GCTTTTGCAATGTCATTTACTG'),         # Set 2
                ('ATAGCGTGATTTTGCGGTT', 'CTAAATATAAATAATGTTGTCACTTGGA'),   # Set 3
                ('CAACGCAATGGTTTCATCCA', 'GCTTTTGCAATGTCATTTACTG')         # Set 4
            ]
            
            roi = None
            
            for fwd, rev in primer_sets:
                for seq_id, seq in sequences.items():
                    amplicon = self.simulate_pcr(seq, fwd, rev)
                    if amplicon:
                        roi = amplicon
                        break
                    seq_rc = str(Seq(seq).reverse_complement())
                    amplicon = self.simulate_pcr(seq_rc, fwd, rev)
                    if amplicon:
                        roi = amplicon
                        break
                if roi: break

            if not roi:
                return {"spa_type": "Negative", "spa_repeats": "-"}

            repeats_map = self.load_repeats()
            found_repeats = []
            current_pos = 0

            while current_pos < len(roi):
                matched = False
                for size in [24, 21, 27, 30, 33]:
                    if current_pos + size > len(roi): continue
                    chunk = roi[current_pos : current_pos+size]
                    
                    if chunk in repeats_map:
                        found_repeats.append(repeats_map[chunk]) 
                        current_pos += size
                        matched = True
                        break
                
                if not matched:
                    current_pos += 1

            if not found_repeats:
                return {"spa_type": "Unknown (No repeats)", "spa_repeats": "-"}

            pattern_str = "-".join(found_repeats)
            types_map = self.load_types()
            final_type = types_map.get(pattern_str, "Novel")
            
            if final_type == "Novel":
                final_type = f"Novel ({pattern_str})"

            return {
                "spa_type": final_type,
                "spa_repeats": pattern_str
            }

        except Exception as e:
            return {"spa_type": "Error", "spa_repeats": str(e)}

    def simulate_pcr(self, seq, fwd, rev):
        rev_rc = str(Seq(rev).reverse_complement())
        
        start_idx = seq.find(fwd)
        if start_idx == -1: return None
        
        end_idx = seq.find(rev_rc, start_idx)
        if end_idx == -1: return None

        region_start = start_idx + len(fwd)
        region_end = end_idx
        
        dist = region_end - region_start
        if dist < 50 or dist > 5000: return None
        
        return seq[region_start:region_end]