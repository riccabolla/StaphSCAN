import io
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


class Module:
    def __init__(self, min_id=100.0, min_cov=100.0, db_dir=None):
        self.name = "mlst"

        self.module_dir = Path(__file__).parent

        if db_dir:
            self.data_dir = Path(db_dir) / "mlst"
            self.data_dir.mkdir(parents=True, exist_ok=True)
        else:
            self.data_dir = self.module_dir / "data"

        self.db_profiles = self.data_dir / "profiles.tsv"
        self.alleles_fasta = self.data_dir / "alleles.fasta"

        # Kept for backward compatibility with existing installations.
        # It is no longer used for allele detection.
        self.refs_fasta = self.data_dir / "refs.fasta"

        self.loci = [
            "arcC",
            "aroE",
            "glpF",
            "gmk",
            "pta",
            "tpi",
            "yqiL",
        ]

        self.min_identity = float(min_id)
        self.min_coverage = float(min_cov)

        self.prof_df = None
        self.allele_map = {locus: {} for locus in self.loci}

        self._load_profiles()
        self._load_allele_map()

    def _load_profiles(self):
        """Load PubMLST ST profiles."""

        if not self.db_profiles.exists():
            return

        dtype_map = {locus: "str" for locus in self.loci}
        dtype_map["ST"] = "str"

        try:
            self.prof_df = pd.read_csv(
                self.db_profiles,
                sep="\t",
                dtype=dtype_map,
            ).fillna("-")

        except Exception as e:
            print(
                f"Warning: Could not load profiles.tsv: {e}",
                file=sys.stderr,
            )

    def _load_allele_map(self):
        """
        Load all allele sequences into memory.
        """

        if not self.alleles_fasta.exists():
            return

        try:
            for record in SeqIO.parse(self.alleles_fasta, "fasta"):
                record_id = record.id

                for locus in self.loci:
                    prefix = f"{locus}_"

                    if record_id.startswith(prefix):
                        allele_id = record_id[len(prefix):]
                        sequence = str(record.seq).upper()

                        self.allele_map[locus][allele_id] = sequence
                        break

        except Exception as e:
            print(
                f"Error loading allele map: {e}",
                file=sys.stderr,
            )


    def _source_fastas(self):
        """Return the downloaded PubMLST locus FASTA files."""

        return [
            self.data_dir / f"{locus}.fas"
            for locus in self.loci
        ]

    def _alleles_fasta_is_current(self):
        """
        Check whether the combined allele FASTA is newer than all
        downloaded PubMLST locus FASTAs.
        """

        if not self.alleles_fasta.exists():
            return False

        combined_mtime = self.alleles_fasta.stat().st_mtime

        for fasta in self._source_fastas():
            if not fasta.exists():
                return False

            if fasta.stat().st_mtime > combined_mtime:
                return False

        return True

    def _build_alleles_fasta(self):
        """
        Build a single FASTA containing every PubMLST allele.
        """

        combined_records = []

        for locus in self.loci:
            fasta_path = self.data_dir / f"{locus}.fas"

            if not fasta_path.exists():
                print(
                    f"Error: Missing {fasta_path}",
                    file=sys.stderr,
                )
                return False

            try:
                for record in SeqIO.parse(fasta_path, "fasta"):
                    # Extract the numeric/string allele identifier.
                    raw_id = record.id

                    if raw_id.startswith(f"{locus}_"):
                        allele_id = raw_id[len(locus) + 1:]
                    elif raw_id.startswith(locus):
                        allele_id = raw_id[len(locus):].lstrip("_-")
                    else:
                        allele_id = raw_id

                    if not allele_id:
                        continue

                    sequence = str(record.seq).upper()

                    header = f"{locus}_{allele_id}"

                    combined_records.append(
                        f">{header}\n{sequence}\n"
                    )

            except Exception as e:
                print(
                    f"Error reading {fasta_path}: {e}",
                    file=sys.stderr,
                )
                return False

        try:
            with open(self.alleles_fasta, "w") as handle:
                handle.write("".join(combined_records))

        except Exception as e:
            print(
                f"Error writing {self.alleles_fasta}: {e}",
                file=sys.stderr,
            )
            return False

        return True

    def check_db(self):
        """
        Check that the MLST database is available and rebuild the
        combined allele FASTA if necessary.
        """

        if not self.db_profiles.exists():
            print(
                f"Error: Missing MLST profile database: "
                f"{self.db_profiles}",
                file=sys.stderr,
            )
            return False

        for fasta in self._source_fastas():
            if not fasta.exists():
                print(
                    f"Error: Missing MLST allele file: {fasta}",
                    file=sys.stderr,
                )
                return False

        if not self._alleles_fasta_is_current():
            print("Building MLST allele database...")

            if not self._build_alleles_fasta():
                return False

            self.allele_map = {locus: {} for locus in self.loci}
            self._load_allele_map()

        return True

    def _run_blast(self, assembly_path):
        """
        BLAST every PubMLST allele against the assembly.
        """

        cmd = [
            "blastn",
            "-task",
            "megablast",

            "-query",
            str(self.alleles_fasta),

            "-subject",
            str(assembly_path),

            "-outfmt",
            (
                "6 "
                "qseqid "
                "sseqid "
                "pident "
                "length "
                "qlen "
                "qcovhsp "
                "bitscore "
                "qstart "
                "qend "
                "sstart "
                "send"
            ),

            "-perc_identity",
            str(self.min_identity),

            "-qcov_hsp_perc",
            str(self.min_coverage),

            "-dust",
            "no",
            "-max_target_seqs",
            "1",
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )

        except FileNotFoundError:
            raise RuntimeError(
                "blastn was not found. "
                "Please make sure BLAST+ is installed and available "
                "in PATH."
            )

        if result.returncode != 0:
            raise RuntimeError(
                f"BLAST failed with exit code {result.returncode}: "
                f"{result.stderr.strip()}"
            )

        if not result.stdout.strip():
            return pd.DataFrame(
                columns=[
                    "qseqid",
                    "sseqid",
                    "pident",
                    "length",
                    "qlen",
                    "qcovhsp",
                    "bitscore",
                    "qstart",
                    "qend",
                    "sstart",
                    "send",
                ]
            )

        columns = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "qcovhsp",
            "bitscore",
            "qstart",
            "qend",
            "sstart",
            "send",
        ]

        return pd.read_csv(
            io.StringIO(result.stdout),
            sep="\t",
            names=columns,
        )

    def _extract_locus_and_allele(self, qseqid):
        """
        Convert a BLAST query ID 
        """

        qseqid = str(qseqid)

        for locus in self.loci:
            prefix = f"{locus}_"

            if qseqid.startswith(prefix):
                return locus, qseqid[len(prefix):]

        return None, None

    def _select_allele(self, locus, hits):
        """
        Select the best allele hit for a locus.
        """

        if hits.empty:
            return "-", "-"

        hits = hits.copy()

        # Convert numeric BLAST columns explicitly.
        for column in [
            "pident",
            "length",
            "qlen",
            "qcovhsp",
            "bitscore",
        ]:
            hits[column] = pd.to_numeric(
                hits[column],
                errors="coerce",
            )

        # Full-length means that the complete PubMLST allele was aligned.
        hits["full_length"] = (
            (hits["length"] >= hits["qlen"]) &
            (hits["qcovhsp"] >= 99.999)
        )

        hits = hits.sort_values(
            by=[
                "full_length",
                "pident",
                "qcovhsp",
                "length",
                "bitscore",
            ],
            ascending=[
                False,
                False,
                False,
                False,
                False,
            ],
        )

        best = hits.iloc[0]

        qseqid = str(best["qseqid"])
        _, allele_id = self._extract_locus_and_allele(qseqid)

        if allele_id is None:
            return "-", "-"

        # Full-length hit satisfying the requested identity threshold.
        if (
            bool(best["full_length"])
            and best["pident"] >= self.min_identity
        ):
            return allele_id, allele_id

        # Partial hit.
        if not bool(best["full_length"]):
            return "Partial", "-"

        # Full-length but below the requested identity threshold.
        # This branch is mainly relevant when min_id < 100.
        return f"{allele_id}*", "-"


    def run(self, assembly_path):

        result = {"ST": "Unknown"}

        for locus in self.loci:
            result[locus] = "-"

        try:
            if not self.check_db():
                result["ST"] = "DB_Error"
                return result

            blast_df = self._run_blast(assembly_path)

            if blast_df.empty:
                return result

            detected_profile = {}

            for locus in self.loci:

                # Allele query IDs are locus-specific.
                locus_mask = (
                    blast_df["qseqid"]
                    .astype(str)
                    .str.startswith(f"{locus}_")
                )

                locus_hits = blast_df[locus_mask]

                allele_call, profile_allele = self._select_allele(
                    locus,
                    locus_hits,
                )

                result[locus] = allele_call
                detected_profile[locus] = profile_allele

            result["ST"] = self.resolve_st(detected_profile)

            return result

        except Exception as e:
            print(
                f"Error in MLST module: {e}",
                file=sys.stderr,
            )
            result["ST"] = "Error"
            return result

    def resolve_st(self, observed_profile):
        """
        Resolve the detected 7-locus allele profile against profiles.tsv.
        """

        if self.prof_df is None:
            return "DB_Error"

        try:
            obs_values = [
                str(observed_profile.get(locus, "-"))
                for locus in self.loci
            ]

            # Compare each observed allele against every PubMLST profile.
            matches_mask = (
                self.prof_df[self.loci] == obs_values
            )

            mismatch_counts = (
                len(self.loci)
                - matches_mask.sum(axis=1)
            )

            min_mismatches = int(mismatch_counts.min())

            if min_mismatches > 2:
                return "Undefined"

            best_match_idx = mismatch_counts.idxmin()

            best_st = str(
                self.prof_df.at[best_match_idx, "ST"]
            )

            if min_mismatches == 0:
                return f"ST{best_st}"

            elif min_mismatches == 1:
                return f"ST{best_st}-1LV"

            elif min_mismatches == 2:
                return f"ST{best_st}-2LV"

            return "Undefined"

        except Exception as e:
            return f"Error({e})"