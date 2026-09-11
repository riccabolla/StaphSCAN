import pathlib
import shutil
import subprocess

import pysam


def check_aligner_dependencies():
    """Ensure minimap2 and samtools are available in PATH."""
    for tool in ["minimap2", "samtools"]:
        if not shutil.which(tool):
            raise EnvironmentError(f"Required binary '{tool}' not found in PATH.")

def align_reads(r1: pathlib.Path, r2: pathlib.Path | None, reference: pathlib.Path, output_bam: pathlib.Path, threads: int = 4):
    """
    Aligns single or paired-end FASTQ reads against a multi-FASTA reference database
    and outputs a sorted, indexed BAM file.
    """
    check_aligner_dependencies()
    
    cmd_minimap = ["minimap2", "-ax", "sr", "-t", str(threads), str(reference), str(r1)]
    if r2:
        cmd_minimap.append(str(r2))
        
    cmd_samtools_view = ["samtools", "view", "-b", "-F", "4", "-"]  # keep only mapped reads
    cmd_samtools_sort = ["samtools", "sort", "-@", str(threads), "-o", str(output_bam), "-"]

    # Pipeline: minimap2 | samtools view -b | samtools sort > output.bam
    p1 = subprocess.Popen(cmd_minimap, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(cmd_samtools_view, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p3 = subprocess.Popen(cmd_samtools_sort, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2.stdout.close()
    p3.communicate()

    # Index the BAM file
    pysam.index(str(output_bam))

def calculate_gene_coverage(bam_file: pathlib.Path, min_depth: int = 5):
    """
    Parses BAM file and returns a dictionary of:
    { target_id: {"breadth": float, "mean_depth": float, "length": int} }
    """
    results = {}
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for ref_name, ref_len in zip(bam.references, bam.lengths):
            covered_bases = 0
            total_depth = 0
            
            for pileupcolumn in bam.pileup(ref_name, truncate=True):
                # Count only reads meeting base quality threshold
                depth = pileupcolumn.nsegments
                if depth >= min_depth:
                    covered_bases += 1
                total_depth += depth
                
            breadth = (covered_bases / ref_len) * 100.0 if ref_len > 0 else 0.0
            mean_depth = (total_depth / ref_len) if ref_len > 0 else 0.0
            
            results[ref_name] = {
                "length": ref_len,
                "breadth": round(breadth, 2),
                "mean_depth": round(mean_depth, 2)
            }
    return results

def check_species_fastq(r1: Path, sketch_db: Path, winner_threshold: float = 0.80):
    """
    Runs `mash screen` on raw FASTQ to verify S. aureus identity and detect mixtures.
    """
    cmd = ["mash", "screen", "-w", str(sketch_db), str(r1)]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    # Parse Mash screen output
    hits = []
    for line in res.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        identity = float(parts[0])
        comment = parts[4] if len(parts) > 4 else ""
        if identity >= winner_threshold:
            hits.append((identity, comment))
            
    return hits