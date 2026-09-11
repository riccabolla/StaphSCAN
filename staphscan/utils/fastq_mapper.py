import collections
import pathlib
import shutil
import subprocess

import pysam
from Bio import SeqIO


def build_master_db(modules_dir: Path, out_fasta: Path) -> None:
    """
    Finds all reference FASTA files and concatenates them.
    Ensures all FASTA headers are strictly unique to prevent samtools crashes.
    """
    seq_count = 0
    seen_ids = set() # Track unique sequence IDs
    
    with open(out_fasta, "w") as out:
        for fasta_file in modules_dir.rglob("*.fasta"):
            if fasta_file.name == "alleles.fasta":
                continue 
            
            from Bio import SeqIO
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id not in seen_ids:
                    out.write(f">{record.id}\n{str(record.seq)}\n")
                    seen_ids.add(record.id)
                    seq_count += 1
                
    if seq_count == 0:
        raise RuntimeError(f"Master database is empty! Could not find reference FASTAs in {modules_dir}")

def align_reads(r1: Path, r2: Path | None, ref_db: Path, out_bam: Path, threads: int = 4) -> None:
    """
    Aligns reads to the master reference and generates a sorted BAM.
    Captures and reports standard error if the pipeline fails.
    """
    check_dependencies()
    
    # base minimap command
    cmd_minimap = ["minimap2", "-ax", "sr", "-t", str(threads), str(ref_db), str(r1)]
    if r2:
        cmd_minimap.append(str(r2))
        
    # Run the pipeline through bash
    pipe_cmd = f"{' '.join(cmd_minimap)} | samtools view -b -F 4 - | samtools sort -@ {threads} -o {out_bam} -"
    
    res = subprocess.run(pipe_cmd, shell=True, capture_output=True, text=True)
    
    # prints specific error when failing
    if res.returncode != 0 or not out_bam.exists():
        raise RuntimeError(f"Alignment pipeline failed!\nTerminal Error Log:\n{res.stderr}")

    pysam.index(str(out_bam))

def check_dependencies():
    for tool in ["minimap2", "samtools"]:
        if not shutil.which(tool):
            raise EnvironmentError(f"Required binary '{tool}' not found in PATH.")

from collections import Counter
from typing import Any, Dict  # noqa: UP035

import pysam


def extract_consensus_from_bam(bam_file: Path, min_depth: int = 5, min_breadth: float = 80.0) -> Dict[str, Any]:
    """
    Parses the BAM file, calculates coverage, and builds a 
    consensus DNA sequence for genes passing the threshold.
    """
    results = {}
    
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for ref_name, ref_len in zip(bam.references, bam.lengths):
            
            # Pre-fill the consensus array with Ns to match the exact reference length
            consensus_list = ["N"] * ref_len
            covered_bases = 0
            total_depth = 0
            
            # Iterate through the pileup and process each column immediately to save mem
            for pileupcolumn in bam.pileup(ref_name, truncate=True):
                pos = pileupcolumn.reference_pos
                depth = pileupcolumn.nsegments
                
                if depth >= min_depth:
                    covered_bases += 1
                    total_depth += depth
                    
                    bases = []
                    # Read the exact nucleotides aligned at this position
                    for pileupread in pileupcolumn.pileups:
                        # Skip deletions or reference skips
                        if not pileupread.is_del and not pileupread.is_refskip:
                            query_pos = pileupread.query_position
                            if query_pos is not None:
                                bases.append(pileupread.alignment.query_sequence[query_pos])
                    
                    if bases:
                        # Call the majority consensus base
                        most_common = Counter(bases).most_common(1)[0][0]
                        consensus_list[pos] = most_common
            
            breadth = (covered_bases / ref_len) * 100.0 if ref_len > 0 else 0.0
            mean_depth = (total_depth / ref_len) if ref_len > 0 else 0.0
            
            if breadth >= min_breadth:
                results[ref_name] = {
                    "dna": "".join(consensus_list),
                    "coverage": round(breadth, 2),
                    "mean_depth": round(mean_depth, 2)
                }
                
    return results