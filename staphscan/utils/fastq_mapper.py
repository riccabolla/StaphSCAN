from __future__ import annotations

import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import pysam
from Bio import SeqIO

REFERENCE_EXTENSIONS = {".fasta", ".fas", ".fa"}

# Files that are not intended to be generic FASTQ mapping targets.
EXCLUDED_REFERENCE_NAMES = {
    "alleles.fasta",
    "refs.fasta",
}


def check_dependencies() -> None:
    """Check external programs required by the FASTQ workflow."""
    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            raise EnvironmentError(
                f"Required binary '{tool}' not found in PATH."
            )


def infer_sample_name(r1: Path) -> str:
    """
    Infer a sample name from common FASTQ naming conventions.

    Examples:
        sample_R1.fastq.gz      -> sample
        sample_1.fastq.gz       -> sample
        sample_R1_001.fastq.gz  -> sample_001
    """
    name = r1.name

    for suffix in (
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq",
    ):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    for token in ("_R1", "_r1", "_1"):
        if name.endswith(token):
            name = name[: -len(token)]
            break

    return name


def _iter_reference_files(
    data_dir: Path,
    reference_patterns: list[str] | tuple[str, ...] | None = None,
):
    """
    Yield reference files explicitly selected by the module, or otherwise
    all supported FASTA extensions under data_dir.

    The fallback includes .fas because the MLST allele files use .fas.
    """
    data_dir = Path(data_dir)

    if reference_patterns:
        files = []

        for pattern in reference_patterns:
            files.extend(sorted(data_dir.glob(pattern)))

        seen = set()

        for path in files:
            path = path.resolve()

            if path in seen or not path.is_file():
                continue

            seen.add(path)
            yield path

        return

    for path in sorted(data_dir.rglob("*")):
        if not path.is_file():
            continue

        if path.suffix.lower() not in REFERENCE_EXTENSIONS:
            continue

        if path.name.lower() in EXCLUDED_REFERENCE_NAMES:
            continue

        yield path


def build_module_db(
    data_dir: Path,
    out_fasta: Path,
    reference_patterns: list[str] | tuple[str, ...] | None = None,
) -> bool:
    """
    Build a FASTQ mapping database from module reference sequences.

    The fallback intentionally supports .fas as well as .fasta/.fa so that
    PubMLST locus files are included.
    """
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    seq_count = 0
    seen_ids: set[str] = set()

    with out_fasta.open("w") as out:
        for fasta_file in _iter_reference_files(
            Path(data_dir),
            reference_patterns=reference_patterns,
        ):
            for record in SeqIO.parse(str(fasta_file), "fasta"):
                if record.id in seen_ids:
                    raise ValueError(
                        "Duplicate reference ID detected while building "
                        f"FASTQ database: {record.id}"
                    )

                seen_ids.add(record.id)
                out.write(f">{record.id}\n{str(record.seq)}\n")
                seq_count += 1

    return seq_count > 0


def _run_pipeline(
    minimap_cmd: list[str],
    output_bam: Path,
    threads: int,
) -> None:
    """
    Run:

        minimap2 -> samtools view -> samtools sort
    """
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    p1 = subprocess.Popen(
        minimap_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    p2 = subprocess.Popen(
        ["samtools", "view", "-b", "-"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if p1.stdout is not None:
        p1.stdout.close()

    p3 = subprocess.Popen(
        [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(output_bam),
            "-",
        ],
        stdin=p2.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if p2.stdout is not None:
        p2.stdout.close()

    _, sort_stderr = p3.communicate()
    _, view_stderr = p2.communicate()
    _, minimap_stderr = p1.communicate()

    if p1.returncode != 0:
        raise RuntimeError(
            "minimap2 failed:\n"
            + minimap_stderr.decode(errors="replace")
        )

    if p2.returncode != 0:
        raise RuntimeError(
            "samtools view failed:\n"
            + view_stderr.decode(errors="replace")
        )

    if p3.returncode != 0:
        raise RuntimeError(
            "samtools sort failed:\n"
            + sort_stderr.decode(errors="replace")
        )

    if not output_bam.exists() or output_bam.stat().st_size == 0:
        raise RuntimeError(
            f"Alignment pipeline completed but BAM was not created: "
            f"{output_bam}"
        )


def align_reads(
    r1: Path,
    r2: Path | None,
    ref_db: Path,
    out_bam: Path,
    threads: int = 4,
    secondary: bool = True,
) -> None:
    """
    Align single- or paired-end short reads against a reference FASTA.

    secondary=True retains secondary mappings so that highly similar
    candidate references can be detected (pass 1).

    secondary=False suppresses them (used for per-candidate remapping).
    """
    check_dependencies()

    if threads < 1:
        raise ValueError("threads must be >= 1")

    r1 = Path(r1)
    ref_db = Path(ref_db)
    out_bam = Path(out_bam)

    if not r1.exists():
        raise FileNotFoundError(f"FASTQ file not found: {r1}")

    if r2 is not None and not Path(r2).exists():
        raise FileNotFoundError(f"FASTQ file not found: {r2}")

    if not ref_db.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_db}")

    cmd = [
        "minimap2",
        "-ax",
        "sr",
        "-t",
        str(threads),
    ]

    if secondary:
        cmd.extend(["-p", "0.8", "-N", "20"])
    else:
        cmd.append("--secondary=no")

    cmd.append(str(ref_db))
    cmd.append(str(r1))

    if r2 is not None:
        cmd.append(str(r2))

    _run_pipeline(cmd, out_bam, threads)

    pysam.index(str(out_bam))


def _family_from_reference(ref_name: str) -> str:
    """Extract a gene family from a reference ID."""
    return ref_name.split("_", 1)[0]


def select_candidate_targets(
    bam_file: Path,
    max_candidates_per_family: int = 3,
    min_reads: int = 3,
) -> dict[str, list[str]]:
    """
    Select up to N candidate references per gene family from pass-1 BAM.

    Candidates are ranked by aggregate primary-alignment AS score.

    Final selection is based on consensus identity + breadth in the
    module's FASTQ section.
    """
    if max_candidates_per_family < 1:
        raise ValueError("max_candidates_per_family must be >= 1")

    if min_reads < 1:
        raise ValueError("min_reads must be >= 1")

    family_scores: dict[str, dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    family_read_counts: dict[str, dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )

    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            if read.is_secondary or read.is_supplementary:
                continue

            if read.reference_name is None:
                continue

            ref_name = read.reference_name
            family = _family_from_reference(ref_name)

            try:
                score = float(read.get_tag("AS"))
            except KeyError:
                score = float(read.query_alignment_length or 0)

            family_scores[family][ref_name] += score
            family_read_counts[family][ref_name] += 1

    candidates: dict[str, list[str]] = {}

    for family, refs in family_scores.items():
        eligible = [
            ref
            for ref in refs
            if family_read_counts[family][ref] >= min_reads
        ]

        if not eligible:
            eligible = list(refs)

        eligible.sort(
            key=lambda ref: (
                refs[ref],
                family_read_counts[family][ref],
            ),
            reverse=True,
        )

        candidates[family] = eligible[:max_candidates_per_family]

    return candidates


def build_focused_db(
    data_dir: Path,
    candidate_targets: dict[str, list[str]],
    out_fasta: Path,
    reference_patterns: list[str] | tuple[str, ...] | None = None,
) -> bool:
    """Build a FASTA containing all selected candidate references."""
    wanted_ids = {
        ref_id
        for ref_ids in candidate_targets.values()
        for ref_id in ref_ids
    }

    if not wanted_ids:
        return False

    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    seq_count = 0
    written_ids: set[str] = set()

    with out_fasta.open("w") as out:
        for fasta_file in _iter_reference_files(
            Path(data_dir),
            reference_patterns=reference_patterns,
        ):
            for record in SeqIO.parse(str(fasta_file), "fasta"):
                if record.id not in wanted_ids:
                    continue

                if record.id in written_ids:
                    continue

                out.write(f">{record.id}\n{str(record.seq)}\n")
                written_ids.add(record.id)
                seq_count += 1

    missing = wanted_ids - written_ids

    if missing:
        raise ValueError(
            "Selected FASTQ reference IDs were not found in the module "
            f"reference files: {', '.join(sorted(missing))}"
        )

    return seq_count > 0


def _usable_bases_at_column(
    pileupcolumn,
    min_mapq: int,
    min_baseq: int,
) -> list[str]:
    """
    Return high-quality nucleotide observations at one reference position.

    Secondary, supplementary, duplicate, deleted and reference-skipped
    observations are excluded.
    """
    bases = []

    for pileup_read in pileupcolumn.pileups:
        alignment = pileup_read.alignment

        if (
            alignment.is_unmapped
            or alignment.is_secondary
            or alignment.is_supplementary
            or alignment.is_duplicate
        ):
            continue

        if alignment.mapping_quality < min_mapq:
            continue

        if pileup_read.is_del or pileup_read.is_refskip:
            continue

        query_pos = pileup_read.query_position

        if query_pos is None:
            continue

        if (
            alignment.query_qualities is not None
            and alignment.query_qualities[query_pos] < min_baseq
        ):
            continue

        base = alignment.query_sequence[query_pos].upper()

        if base in {"A", "C", "G", "T"}:
            bases.append(base)

    return bases


def extract_consensus_from_bam(
    bam_file: Path,
    reference_fasta: Path,
    min_depth: int = 3,
    min_mapq: int = 0,
    min_baseq: int = 20,
    min_breadth: float = 0.0,
) -> dict[str, dict[str, Any]]:
    """
    Build a reference-oriented consensus for each reference in a BAM.

    Returned metrics:

        dna
            Reference-oriented consensus sequence. Positions below
            min_depth remain N.

        coverage
            Percentage of reference positions with >= min_depth usable bases.

        mean_depth
            Mean usable depth across the entire reference.

        pident
            Percentage identity between consensus and reference among
            covered positions only.

        matches / mismatches
            Covered positions matching / differing from the reference.

    Identity is measured only where sufficient evidence exists, so it is
    not inflated by uncovered positions. Callers must therefore consider
    identity together with coverage.

    Pileup uses ignore_orphans=False and ignore_overlaps=False so that reads
    whose mate falls outside a short gene, and overlapping mates, are all
    counted. Otherwise gene ends lose depth and breadth is underestimated.
    """
    if min_depth < 1:
        raise ValueError("min_depth must be >= 1")

    if min_mapq < 0:
        raise ValueError("min_mapq must be >= 0")

    if min_baseq < 0:
        raise ValueError("min_baseq must be >= 0")

    reference_sequences = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(str(reference_fasta), "fasta")
    }

    results: dict[str, dict[str, Any]] = {}

    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        for ref_name, ref_len in zip(bam.references, bam.lengths):
            reference = reference_sequences.get(ref_name)

            if reference is None:
                raise ValueError(
                    f"Reference '{ref_name}' is present in BAM but absent "
                    f"from {reference_fasta}"
                )

            if len(reference) != ref_len:
                raise ValueError(
                    f"Reference length mismatch for '{ref_name}': "
                    f"BAM={ref_len}, FASTA={len(reference)}"
                )

            consensus = ["N"] * ref_len
            covered_bases = 0
            total_depth = 0
            matches = 0
            mismatches = 0

            for pileupcolumn in bam.pileup(
                ref_name,
                truncate=True,
                stepper="all",
                min_base_quality=0,
                ignore_orphans=False,
                ignore_overlaps=False,
                max_depth=100000,
            ):
                pos = pileupcolumn.reference_pos

                if pos < 0 or pos >= ref_len:
                    continue

                bases = _usable_bases_at_column(
                    pileupcolumn,
                    min_mapq=min_mapq,
                    min_baseq=min_baseq,
                )

                depth = len(bases)
                total_depth += depth

                if depth < min_depth:
                    continue

                covered_bases += 1

                counts = defaultdict(int)
                for base in bases:
                    counts[base] += 1

                consensus_base = max(
                    counts,
                    key=lambda base: (counts[base], base),
                )

                consensus[pos] = consensus_base

                if consensus_base == reference[pos]:
                    matches += 1
                else:
                    mismatches += 1

            breadth = covered_bases / ref_len * 100.0 if ref_len else 0.0
            mean_depth = total_depth / ref_len if ref_len else 0.0
            identity = matches / covered_bases * 100.0 if covered_bases else 0.0

            if breadth < min_breadth:
                continue

            results[ref_name] = {
                "dna": "".join(consensus),
                "coverage": round(breadth, 2),
                "mean_depth": round(mean_depth, 2),
                "pident": round(identity, 2),
                "matches": matches,
                "mismatches": mismatches,
                "length": ref_len,
            }

    return results


def consensus_per_candidate(
    r1: Path,
    r2: Path | None,
    focused_db: Path,
    work_dir: Path,
    threads: int = 4,
    min_depth: int = 3,
    min_mapq: int = 0,
    min_baseq: int = 20,
) -> dict[str, dict[str, Any]]:
    """
    Map the reads to each candidate reference separately and build a
    consensus for each.

    Mapping all candidates in one run lets near-identical alleles compete:
    reads in shared regions get MAPQ 0 and are discarded, so even the true
    allele loses breadth. With one reference per run there is no ambiguity,
    so MAPQ is not informative and min_mapq should normally be 0. Paralog
    cross-mapping is still rejected downstream through consensus identity.

    Cost: one minimap2 run per candidate, each reading the FASTQ files.
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, dict[str, Any]] = {}

    for record in SeqIO.parse(str(focused_db), "fasta"):
        safe = re.sub(r"[^\w.-]", "_", record.id)
        ref = work_dir / f"{safe}.fasta"
        bam = work_dir / f"{safe}.bam"

        SeqIO.write(record, str(ref), "fasta")

        align_reads(r1, r2, ref, bam, threads=threads, secondary=False)

        results.update(
            extract_consensus_from_bam(
                bam,
                ref,
                min_depth=min_depth,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
            )
        )

    return results


def select_best_consensus_per_family(
    consensus_dict: dict[str, dict[str, Any]],
    min_id: float = 90.0,
    min_cov: float = 80.0,
) -> dict[str, dict[str, Any]]:
    """
    Select the best-supported reference per family.

    Ranking:
      1. candidates passing min_id and min_cov first
      2. breadth
      3. identity
      4. mean depth

    Modules may still apply their own thresholds afterwards.
    """

    def _key(stats: dict[str, Any]) -> tuple:
        cov = stats.get("coverage", 0.0)
        pid = stats.get("pident", 0.0)
        return (
            pid >= min_id and cov >= min_cov,
            cov,
            pid,
            stats.get("mean_depth", 0.0),
        )

    best: dict[str, dict[str, Any]] = {}

    for ref_name, stats in consensus_dict.items():
        family = _family_from_reference(ref_name)
        current = best.get(family)

        if current is None or _key(stats) > _key(current):
            best[family] = {**stats, "reference": ref_name}

    return best