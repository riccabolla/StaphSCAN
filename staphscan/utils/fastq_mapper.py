import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict  # noqa: UP035

import pysam


def check_dependencies():
    for tool in ["minimap2", "samtools"]:
        if not shutil.which(tool):
            raise EnvironmentError(f"Required binary '{tool}' not found in PATH.")

def build_module_db(data_dir: Path, out_fasta: Path) -> bool:
    from Bio import SeqIO
    seq_count = 0
    with open(out_fasta, "w") as out:
        for fasta_file in data_dir.rglob("*.fasta"):
            if fasta_file.name == "alleles.fasta":
                continue 
            for record in SeqIO.parse(fasta_file, "fasta"):
                out.write(f">{record.id}\n{str(record.seq)}\n")
                seq_count += 1
    return seq_count > 0

def align_reads(r1: Path, r2: Path | None, ref_db: Path, out_bam: Path, threads: int = 4) -> None:
    check_dependencies()
    
    # keep up to 100 alignments down to 50% score
    cmd_minimap = ["minimap2", "-ax", "sr", "-t", str(threads), "-N", "100", "-p", "0.5", str(ref_db), str(r1)]
    if r2: 
        cmd_minimap.append(str(r2))
        
    pipe_cmd = f"{' '.join(cmd_minimap)} | samtools view -b - | samtools sort -@ {threads} -o {out_bam} -"
    res = subprocess.run(pipe_cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0 or not out_bam.exists():
        raise RuntimeError(f"Alignment pipeline failed!\n{res.stderr}")
    pysam.index(str(out_bam))

def extract_consensus_from_bam(bam_file: Path, min_depth: int = 1, min_breadth: float = 40.0) -> Dict[str, Any]:
    results = {}
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for ref_name, ref_len in zip(bam.references, bam.lengths):
            consensus_list = ["N"] * ref_len
            covered_bases = 0
            total_depth = 0
            
            for pileupcolumn in bam.pileup(ref_name, truncate=True, stepper="all", min_base_quality=0):
                pos = pileupcolumn.reference_pos
                depth = pileupcolumn.nsegments
                
                if depth >= min_depth:
                    covered_bases += 1
                    total_depth += depth
                    # extract bases
                    bases = [r.alignment.query_sequence[r.query_position] for r in pileupcolumn.pileups if r.query_position is not None]
                    if bases:
                        consensus_list[pos] = Counter(bases).most_common(1)[0][0]
            
            breadth = (covered_bases / ref_len) * 100.0 if ref_len > 0 else 0.0
            mean_depth = (total_depth / ref_len) if ref_len > 0 else 0.0
            
            if breadth >= min_breadth:
                results[ref_name] = {
                    "dna": "".join(consensus_list),
                    "coverage": round(breadth, 2),
                    "mean_depth": round(mean_depth, 2)
                }
    return results

def select_best_targets(bam_file: Path, min_reads: int = 3) -> Dict[str, str]:
    """
    From a BAM aligned against a multi-allele reference DB, pick the best-supported
    reference (allele) per gene family based on total primary-alignment mapping
    score (minimap2's AS tag). This lets the aligner's own scoring decide which allele a
    sample's reads actually match
    Only primary alignments are counted, since with
    -N 100 -p 0.5 a single read can align to many near-identical alleles

    """
    family_scores: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
    family_read_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            ref_name = read.reference_name
            if ref_name is None:
                continue
            family = ref_name.split("_")[0]
            try:
                score = read.get_tag("AS")  # minimap2 alignment score
            except KeyError:
                score = read.query_alignment_length or 0  # fallback if AS missing
            family_scores[family][ref_name] += score
            family_read_counts[family][ref_name] += 1

    best_targets: Dict[str, str] = {}
    for family, ref_scores in family_scores.items():
        # Require a minimum read count
        candidates = {r: s for r, s in ref_scores.items()
                      if family_read_counts[family][r] >= min_reads}
        if not candidates:
            candidates = ref_scores
        best_targets[family] = max(candidates, key=candidates.get)

    return best_targets

def build_focused_db(data_dir: Path, best_targets: Dict[str, str], out_fasta: Path) -> bool:
    """
    Writes a reference fasta containing only the best targets
    """
    from Bio import SeqIO
    wanted_ids = set(best_targets.values())
    seq_count = 0
    with open(out_fasta, "w") as out:
        for fasta_file in data_dir.rglob("*.fasta"):
            if fasta_file.name == "alleles.fasta":
                continue
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in wanted_ids:
                    out.write(f">{record.id}\n{str(record.seq)}\n")
                    seq_count += 1
    return seq_count > 0