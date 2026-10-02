from dataclasses import dataclass
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Seq import Seq


@dataclass
class GeneHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    slen: int
    qlen: int
    sstart: int
    send: int
    bitscore: float

    @property
    def coverage(self) -> float:
        return (self.length / self.qlen) * 100.0 if self.qlen else 0.0

    @property
    def family(self) -> str:
        return self.qseqid.split("_")[0]


class Module:
    """
    StaphSCAN virulence module.

    Supports:
      - FASTA/assembly analysis through BLASTn
      - FASTQ analysis through a reference-oriented consensus generated
        by the FASTQ mapping layer

    The biological interpretation/scoring is shared between both input modes.
    FASTA mode is the reference behaviour; FASTQ mode is tuned to reproduce it.
    """

    name = "virulence"

    # FASTQ noise floor. Mirrors the BLAST pre-filter used in FASTA mode
    # (-perc_identity 80, -qcov_hsp_perc 40): candidates below it are treated
    # as absent rather than reported as spurious.
    NOISE_MIN_ID = 80.0
    NOISE_MIN_COV = 40.0

    def __init__(self, min_id: float = 90.0, min_cov: float = 80.0):
        self.min_id = float(min_id)
        self.min_cov = float(min_cov)

        self.data_dir = Path(__file__).parent / "data"
        self.target_db = self._find_target_db()

        if self.target_db is None:
            raise FileNotFoundError(
                f"No virulence target FASTA found in {self.data_dir}"
            )

        # Restrict the FASTQ mapper to the targets file only, so other FASTA
        # files placed in data/ are never mapped.
        self.fastq_reference_patterns = [self.target_db.name]

        self.targets = list(SeqIO.parse(self.target_db, "fasta"))
        if not self.targets:
            raise ValueError(f"Virulence target database is empty: {self.target_db}")

        self.target_by_id = {
            record.id: str(record.seq).upper() for record in self.targets
        }

        # Translate nucleotide reference targets once. These proteins are used
        # for the protein confirmation step.
        self.target_proteins = {}
        for record in self.targets:
            seq = str(record.seq).upper().replace("-", "")
            try:
                protein = str(Seq(seq).translate(table=11, to_stop=False))
            except Exception:
                protein = ""
            self.target_proteins[record.id] = protein

        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = "global"
        self.aligner.match_score = 1
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5

    # database

    def _find_target_db(self) -> Path | None:
        """
        Locate the virulence target database.

        Prefer the conventional targets*.fasta naming used by the module.
        Fall back to a single FASTA/FA file if no targets*.fasta exists.
        """
        preferred = sorted(
            p for p in self.data_dir.glob("targets*.fasta")
            if p.is_file()
        )
        if preferred:
            return preferred[0]

        preferred = sorted(
            p for p in self.data_dir.glob("targets*.fas")
            if p.is_file()
        )
        if preferred:
            return preferred[0]

        candidates = sorted(
            p for p in self.data_dir.iterdir()
            if p.is_file() and p.suffix.lower() in {".fasta", ".fas", ".fa"}
        )
        return candidates[0] if len(candidates) == 1 else None

    def check_db(self) -> bool:
        return (
            self.target_db is not None
            and self.target_db.exists()
            and len(self.targets) > 0
        )

    def get_fastq_reference(self) -> Path:
        """
        Return the explicit reference database used by the FASTQ mapper.
        """
        if not self.check_db():
            raise FileNotFoundError(
                f"Virulence database is unavailable: {self.data_dir}"
            )
        return self.target_db

    # biological interpretation

    @staticmethod
    def _is_pair_present(families: set[str], gene_a: str, gene_b: str) -> bool:
        return gene_a in families and gene_b in families

    @staticmethod
    def _classify_families(families: dict[str, str]) -> dict:
        """
        Apply the existing StaphSCAN virulence scoring scheme.

        Score:
          PVL (lukS + lukF) or TSST1         -> 3
          Exfoliative toxins eta/etb/etd/ete -> 2
          Enterotoxins / LukED               -> 1
        """
        pvl = "lukS" in families and "lukF" in families
        tsst = "tst" in families or "tsst1" in families

        exfoliative = sorted(
            families.intersection({"eta", "etb", "etd", "ete"})
        )

        enterotoxins = sorted(
            families.intersection({"sea", "sec", "seh", "selk", "sell", "selq"})
        )

        luked = "lukE" in families and "lukD" in families

        score = 0

        if pvl or tsst:
            score += 3

        if exfoliative:
            score += 2

        if enterotoxins or luked:
            score += 1

        return {
            "vir_score": score,
            "vir_pvl": "Positive" if pvl else "-",
            "vir_tsst": "Positive" if tsst else "-",
            "vir_et": ",".join(exfoliative) if exfoliative else "-",
            "vir_lukED": "Positive" if luked else "-",
            "vir_se": ",".join(enterotoxins) if enterotoxins else "-",
        }

    @staticmethod
    def _empty_result() -> dict:
        return {
            "vir_score": 0,
            "vir_pvl": "-",
            "vir_tsst": "-",
            "vir_et": "-",
            "vir_lukED": "-",
            "vir_se": "-",
            "spurious_virulence_hits": "-",
            "truncated_virulence_hits": "-",
        }

    # fasta mode (unchanged: this is the reference behaviour)

    @staticmethod
    def _run_blast(query: Path, subject: Path) -> list[GeneHit]:
        import shutil
        import subprocess

        if shutil.which("blastn") is None:
            raise EnvironmentError("Required binary 'blastn' not found in PATH.")

        outfmt = (
            "6 qseqid sseqid pident length slen qlen "
            "sstart send bitscore"
        )

        cmd = [
            "blastn",
            "-task",
            "blastn",
            "-query",
            str(query),
            "-subject",
            str(subject),
            "-outfmt",
            outfmt,
            "-perc_identity",
            "80",
            "-qcov_hsp_perc",
            "40",
            "-dust",
            "no",
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"blastn failed with exit code {result.returncode}: "
                f"{result.stderr.strip()}"
            )

        hits = []
        for line in result.stdout.splitlines():
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            try:
                hits.append(
                    GeneHit(
                        qseqid=fields[0],
                        sseqid=fields[1],
                        pident=float(fields[2]),
                        length=int(fields[3]),
                        slen=int(fields[4]),
                        qlen=int(fields[5]),
                        sstart=int(fields[6]),
                        send=int(fields[7]),
                        bitscore=float(fields[8]),
                    )
                )
            except ValueError:
                continue

        return hits

    @staticmethod
    def _extract_gene(
        assembly_seq: str,
        start: int,
        end: int,
    ) -> str:
        """
        Extract a BLAST subject interval.

        BLAST coordinates are 1-based inclusive. Reverse-strand hits are
        reverse-complemented so the returned sequence is in reference
        orientation.
        """
        left = min(start, end) - 1
        right = max(start, end)

        seq = assembly_seq[left:right]

        if start > end:
            seq = str(Seq(seq).reverse_complement())

        return seq

    @staticmethod
    def _translate_best_frame(dna: str, reference_protein: str) -> tuple[str, float]:
        """
        Translate the three forward frames and select the one with the best
        global protein alignment against the reference protein.

        The DNA sequence is expected to already be in reference orientation.
        """
        dna = dna.upper()
        best_protein = ""
        best_score = float("-inf")

        for frame in range(3):
            usable = dna[frame:]
            usable = usable[: len(usable) - (len(usable) % 3)]

            if len(usable) < 3:
                continue

            try:
                protein = str(
                    Seq(usable).translate(table=11, to_stop=False)
                )
            except Exception:
                continue

            if not reference_protein:
                score = 0.0
            else:
                score = float(
                    Module._protein_alignment_score(
                        protein,
                        reference_protein,
                    )
                )

            if score > best_score:
                best_score = score
                best_protein = protein

        return best_protein, best_score

    @staticmethod
    def _protein_alignment_score(
        query: str,
        reference: str,
    ) -> float:
        if not query or not reference:
            return float("-inf")

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        return float(aligner.score(query, reference))

    @staticmethod
    def _protein_identity(
        query: str,
        reference: str,
    ) -> float:
        if not query or not reference:
            return 0.0

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        alignment = aligner.align(query, reference)[0]

        # Identity is calculated from aligned coordinates, since gapped
        # strings are not exposed in all Biopython versions.
        matches = 0
        compared = 0

        for (q_start, q_end), (r_start, r_end) in zip(
            alignment.aligned[0],
            alignment.aligned[1],
        ):
            length = min(q_end - q_start, r_end - r_start)
            for i in range(length):
                compared += 1
                if query[q_start + i] == reference[r_start + i]:
                    matches += 1

        if compared == 0:
            return 0.0

        return (matches / compared) * 100.0

    @staticmethod
    def _best_hits_by_family(hits: list[GeneHit]) -> dict[str, GeneHit]:
        best = {}

        for hit in hits:
            current = best.get(hit.family)

            if current is None:
                best[hit.family] = hit
                continue

            current_key = (current.pident, current.coverage, current.bitscore)
            hit_key = (hit.pident, hit.coverage, hit.bitscore)

            if hit_key > current_key:
                best[hit.family] = hit

        return best

    def run(self, fasta_path: Path) -> dict:
        result = self._empty_result()

        assembly_records = list(SeqIO.parse(fasta_path, "fasta"))
        if not assembly_records:
            raise ValueError(f"No FASTA records found in {fasta_path}")

        assembly_seq = "".join(str(record.seq).upper() for record in assembly_records)

        hits = self._run_blast(self.target_db, fasta_path)

        if not hits:
            return result

        best_hits = self._best_hits_by_family(hits)

        strong_families = set()
        spurious = set()
        truncated = set()

        for family, hit in best_hits.items():
            if hit.pident < self.min_id or hit.coverage < self.min_cov:
                spurious.add(family)
                continue

            reference_protein = self.target_proteins.get(hit.qseqid, "")
            dna = self._extract_gene(
                assembly_seq,
                hit.sstart,
                hit.send,
            )

            protein, _ = self._translate_best_frame(
                dna,
                reference_protein,
            )

            protein_id = self._protein_identity(
                protein,
                reference_protein,
            )

            display_str = family
            if hit.pident < 100.0:
                display_str+= "^"
            else: 
                display_str+= "*"
            if hit.coverage < 100.0:
                display_str+= "?"

            is_strong = hit.pident >= self.min_id and hit.coverage >= self.min_cov
            is_truncated = False

            # A nucleotide hit can satisfy the BLAST threshold while still
            # containing a disruptive protein-level change.
            if reference_protein and protein_id < self.min_id:
                spurious.add(family)
                continue

            # Detect clear premature termination/truncation.
            ref_len = len(reference_protein)
            if ref_len:
                usable_stop = protein.find("*")
                if usable_stop >= 0:
                    retained_pct = (usable_stop / ref_len) * 100.0

                    if retained_pct < self.min_cov:
                        truncated.append(f"{family}-{int(retained_pct)}%(Stop)")
                        is_truncated = True
            if is_truncated:
                continue
            elif is_strong:
                strong_families[family] = display_str
            else:
                spurious.append(display_str)

            #strong_families.add(family)

        result.update(self._classify_families(strong_families))

        result["spurious_virulence_hits"] = (
            ",".join(sorted(spurious)) if spurious else "-"
        )
        result["truncated_virulence_hits"] = (
            ",".join(sorted(truncated)) if truncated else "-"
        )

        return result

    # fastq mode

    @staticmethod
    def _normalise_family(name: str) -> str:
        """Convert a reference/consensus ID to the biological family name."""
        return name.split("_")[0]

    @staticmethod
    def _trim_consensus_to_reference(
        dna: str,
        reference_dna: str,
    ) -> str:
        """
        Force the consensus to the reference length.

        A short consensus is right-padded with Ns and a long one is cut.
        Internal Ns are retained because they represent unresolved bases and
        are not silently converted into matches.
        """
        if not dna:
            return ""

        dna = dna.upper()
        reference_dna = reference_dna.upper()

        if len(dna) < len(reference_dna):
            dna = dna + ("N" * (len(reference_dna) - len(dna)))
        elif len(dna) > len(reference_dna):
            dna = dna[: len(reference_dna)]

        return dna

    @staticmethod
    def _fill_ns_from_reference(dna: str, reference: str) -> str:
        """
        Replace unresolved positions with the reference base.

        Used only before translation, so that a few masked bases do not turn
        codons into 'X' and count as protein mismatches. The nucleotide
        identity already ignores Ns, and the breadth / N-block checks still
        limit how many Ns can reach this point.
        """
        return "".join(
            reference[i] if (base == "N" and i < len(reference)) else base
            for i, base in enumerate(dna)
        )

    @staticmethod
    def _dna_identity(
        consensus: str,
        reference: str,
        min_depth_mask: list[bool] | None = None,
    ) -> float:
        """
        Calculate nucleotide identity over usable consensus positions.

        Ns are excluded rather than counted as mismatches.
        """
        if not consensus or not reference:
            return 0.0

        length = min(len(consensus), len(reference))

        matches = 0
        compared = 0

        for i in range(length):
            if min_depth_mask is not None and i < len(min_depth_mask):
                if not min_depth_mask[i]:
                    continue

            base = consensus[i]
            if base == "N":
                continue

            compared += 1
            if base == reference[i]:
                matches += 1

        return (matches / compared) * 100.0 if compared else 0.0

    def _passes_thresholds(self, pident: float, coverage: float) -> bool:
        return pident >= self.min_id and coverage >= self.min_cov

    def _select_fastq_hits(
        self,
        consensus_dict: dict,
    ) -> dict[str, tuple[str, dict]]:
        """
        Select the best FASTQ-supported target for each virulence family.

        Candidate ranking:
          1. candidates passing both min_id and min_cov first
          2. reference breadth
          3. identity
          4. mean depth

        Passing candidates are ranked first so that a short, high-identity
        fragment cannot hide a full-length allele (FASTA mode avoids this
        because BLAST HSPs already span most of the gene).
        """
        candidates = {}

        for target_id, stats in consensus_dict.items():
            if not isinstance(stats, dict):
                continue

            reference = self.target_by_id.get(target_id)
            if reference is None:
                continue

            dna = self._trim_consensus_to_reference(
                str(stats.get("dna", "")).upper(),
                reference,
            )
            coverage = float(stats.get("coverage", 0.0))
            mean_depth = float(stats.get("mean_depth", 0.0))

            pident = stats.get("pident")
            if pident is None:
                pident = self._dna_identity(dna, reference)
            pident = float(pident)

            key = (
                self._passes_thresholds(pident, coverage),
                coverage,
                pident,
                mean_depth,
            )

            candidate = {
                **stats,
                "dna": dna,
                "pident": pident,
                "coverage": coverage,
                "mean_depth": mean_depth,
            }

            family = self._normalise_family(target_id)
            current = candidates.get(family)

            if current is None or key > current[2]:
                candidates[family] = (target_id, candidate, key)

        return {
            family: (target_id, candidate)
            for family, (target_id, candidate, _) in candidates.items()
        }

    def run_fastq(self, consensus_dict: dict) -> dict:
        """
        Interpret reference-oriented FASTQ consensus sequences.

        Decision order per family (best candidate only):
          1. below the noise floor (cov < 40 or id < 80) -> absent
          2. below min_cov / min_id                       -> spurious
          3. terminal N-block / low retained span         -> truncated
          4. protein identity < min_id                    -> spurious
          5. premature stop                               -> truncated
          6. otherwise                                    -> strong hit
        """
        result = self._empty_result()

        candidates = self._select_fastq_hits(consensus_dict)

        strong_families = set()
        spurious = set()
        truncated = set()

        for family, (target_id, stats) in candidates.items():
            coverage = float(stats.get("coverage", 0.0))
            pident = float(stats.get("pident", 0.0))
            dna = str(stats.get("dna", "")).upper()
            reference = self.target_by_id.get(target_id, "")

            if not reference:
                continue

            # Noise floor: equivalent to the BLAST pre-filter in FASTA mode.
            if coverage < self.NOISE_MIN_COV or pident < self.NOISE_MIN_ID:
                continue

            if not self._passes_thresholds(pident, coverage):
                spurious.add(family)
                continue

            non_n_positions = [i for i, base in enumerate(dna) if base != "N"]

            if not non_n_positions:
                continue

            first = min(non_n_positions)
            last = max(non_n_positions)

            retained_pct = ((last - first + 1) / len(reference)) * 100.0

            if retained_pct < self.min_cov:
                truncated.add(family)
                continue

            # Protein-level sanity check. Unresolved positions are filled
            # from the reference so they do not become 'X' mismatches.
            reference_protein = self.target_proteins.get(target_id, "")

            if reference_protein:
                filled = self._fill_ns_from_reference(dna, reference)
                usable = filled[: len(reference)]
                usable = usable[: len(usable) - (len(usable) % 3)]

                protein, _ = self._translate_best_frame(
                    usable,
                    reference_protein,
                )

                if protein:
                    protein_identity = self._protein_identity(
                        protein,
                        reference_protein,
                    )

                    if protein_identity < self.min_id:
                        spurious.add(family)
                        continue

                    first_stop = protein.find("*")

                    if first_stop >= 0:
                        protein_retained_pct = (
                            first_stop / len(reference_protein)
                        ) * 100.0

                        if protein_retained_pct < self.min_cov:
                            truncated.add(family)
                            continue

            strong_families.add(family)

        result.update(self._classify_families(strong_families))

        result["spurious_virulence_hits"] = (
            ",".join(sorted(spurious)) if spurious else "-"
        )
        result["truncated_virulence_hits"] = (
            ",".join(sorted(truncated)) if truncated else "-"
        )

        return result