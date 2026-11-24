import collections
import gzip
import io
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Union  # noqa: UP035

import pandas as pd


class Module:
    def __init__(self):
        self.name = "assembly"
        self.module_dir = Path(__file__).parent
        self.db_sketch = self.module_dir / "data" / "staph_refs.msh"      
        self.min_size = 2600000 
        self.min_size_fastq = 2500000
        self.max_size = 3100000
        self.max_size_fastq = 3300000
        self.min_n50 = 10000
        self.min_cov_fastq = 10.0

    def check_db(self):
        if not self.db_sketch.exists():
            print(f"Warning: Mash sketch '{self.db_sketch.name}' not found. Skipping assembly module.")
            return False
        return True

    def is_fastq(self, input_path: Union[Path, str, List[Path]]) -> bool:
        if isinstance(input_path, (list, tuple)):
            first_path = Path(input_path[0])
        else:
            first_path = Path(input_path)
        name_lower = first_path.name.lower()
        return any(name_lower.endswith(ext) for ext in (".fastq", ".fq", ".fastq.gz", ".fq.gz"))

    def run(self, assembly_path: Union[Path, str, List[Path]], is_fastq: bool = None, read_type: str = "sr"):
        if not self.check_db():
            return {
                "Species": "Not Analyzed",
                "Mash_distance": "-",
                "Total_size": "-",
                "QC": "SKIPPED (Missing DB)",
                "N_contig": "-",
                "N50": "-"
            }

        if is_fastq is None:
            is_fastq = self.is_fastq(assembly_path)

        # --- FASTQ WORKFLOW ---
        if is_fastq:
            read_files = [Path(p) for p in assembly_path] if isinstance(assembly_path, (list, tuple)) else [Path(assembly_path)]
            species_call, mash_dist, is_mixed, est_size, est_cov = self.run_mash_fastq(read_files, read_type=read_type)

            qc_failures = []
            if est_size > 0:
                if est_size < self.min_size_fastq:
                    qc_failures.append("Too short")
                if est_size > self.max_size_fastq:
                    qc_failures.append("Too long")
                if read_type == "sr" and est_size > self.max_size:
                    qc_failures.append("Too long")
            if 0 < est_cov < self.min_cov_fastq:
                qc_failures.append("Low_Coverage")
            if is_mixed:
                qc_failures.append("Mixed")

            qc_status = "PASS"
            if qc_failures:
                qc_status = f"FAILED ({','.join(qc_failures)})"

            return {
                "Species": species_call,
                "Mash_distance": mash_dist,
                "Total_size": est_size if est_size > 0 else "-",
                "QC": qc_status,
                "N_contig": "-",
                "N50": "-"
            }

        # --- EXISTING FASTA WORKFLOW ---
        fasta_path = Path(assembly_path[0]) if isinstance(assembly_path, (list, tuple)) else Path(assembly_path)
        contig_count, n50, total_size, ambig = self.get_contig_stats(fasta_path)
        species_call, mash_dist, is_mixed = self.run_mash(fasta_path)
        
        qc_failures = []
        if total_size < self.min_size: qc_failures.append("Too short")
        if total_size > self.max_size: qc_failures.append("Too long")
        if n50 < self.min_n50: qc_failures.append("Low_N50")
        if "yes" in ambig: qc_failures.append("Ambiguous_Bases")
        if is_mixed: qc_failures.append("Mixed")
        
        qc_status = "PASS"
        if qc_failures:
            qc_status = f"FAILED ({','.join(qc_failures)})"

        return {
            "Species": species_call,
            "Mash_distance": mash_dist,
            "Total_size": total_size,
            "QC": qc_status,
            "N_contig": contig_count,
            "N50": n50
        }

    def _parse_mash_stdout(self, stdout_str: str):
        if not stdout_str or not stdout_str.strip():
            return "Unknown", "-", False

        rows = []
        for line in stdout_str.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            ref_name = Path(parts[0]).name 
            dist = float(parts[2])
            rows.append((ref_name, dist))

        if not rows:
            return "Unknown", "-", False
            
        rows.sort(key=lambda x: x[1])
        best_match = rows[0]
        best_name = best_match[0]
        best_dist = best_match[1]

        close_matches = [r for r in rows if r[1] <= 0.04]
        is_mixed = len(close_matches) > 1

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

    def run_mash(self, assembly_path: Path):
        try:
            cmd = ["mash", "dist", str(self.db_sketch), str(assembly_path)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            return self._parse_mash_stdout(res.stdout)
        except Exception as e:
            return f"Error ({e})", "-", False

    def _sketch_reads(self, read_files: List[Path], sketch_prefix: Path, min_copies: int) -> Tuple[int, float]:
        sketch_cmd = [
            "mash", "sketch", "-r",
            "-m", str(min_copies),
            "-o", str(sketch_prefix),
            "-"
        ]
        proc = subprocess.Popen(
            sketch_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        try:
            for fpath in read_files:
                open_fn = gzip.open if fpath.name.lower().endswith(".gz") else open
                with open_fn(fpath, "rb") as fh:
                    shutil.copyfileobj(fh, proc.stdin, length=1024 * 1024)
        except BrokenPipeError:
            pass
        finally:
            if proc.stdin:
                try:
                    proc.stdin.close()
                except Exception:
                    pass
                proc.stdin = None

        _, stderr_bytes = proc.communicate()
        stderr_str = stderr_bytes.decode("utf-8", errors="replace")

        est_size = 0
        est_cov = 0.0
        if stderr_str:
            size_match = re.search(r"Estimated genome size:\s*([0-9.eE+-]+)", stderr_str)
            cov_match = re.search(r"Estimated coverage:\s*([0-9.eE+-]+)", stderr_str)
            if size_match:
                est_size = int(round(float(size_match.group(1))))
            if cov_match:
                est_cov = float(cov_match.group(1))

        return est_size, est_cov

    def run_mash_fastq(self, read_files: List[Path], read_type: str = "sr"):
        try:
            with tempfile.TemporaryDirectory() as td:
                sketch_prefix = Path(td) / "sample_reads"
                sketch_file = Path(td) / "sample_reads.msh"

                init_m = 4
                est_size, est_cov = self._sketch_reads(read_files, sketch_prefix, min_copies=init_m)

                adaptive_m = min(6, int(round(est_cov * 0.08)))
                if adaptive_m > init_m:
                    est_size, est_cov = self._sketch_reads(read_files, sketch_prefix, min_copies=adaptive_m)

                if not sketch_file.exists():
                    return "Unknown", "-", False, est_size, est_cov

                dist_cmd = ["mash", "dist", str(self.db_sketch), str(sketch_file)]
                dist_res = subprocess.run(dist_cmd, capture_output=True, text=True)
                species_call, mash_dist, is_mixed = self._parse_mash_stdout(dist_res.stdout)
                return species_call, mash_dist, is_mixed, est_size, est_cov

        except Exception as e:
            return f"Error ({e})", "-", False, 0, 0.0

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

        if not lengths: return 0, 0, 0, "no"

        lengths.sort()
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
        return len(lengths), n50, total, ambig_str