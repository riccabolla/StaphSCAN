import argparse
import importlib
import sys
import tempfile
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pandas as pd


if sys.version_info < (3, 10):
    sys.exit("StaphScan requires Python 3.10+")


def get_version():
    try:
        return version("staphscan")
    except PackageNotFoundError:
        return "unknown"


def get_available_modules():
    modules_dir = Path(__file__).parent / "modules"

    if not modules_dir.exists():
        sys.exit(f"Error: 'modules' directory not found at {modules_dir}")

    return sorted(
        d.name
        for d in modules_dir.iterdir()
        if d.is_dir() and (d / f"{d.name}.py").exists()
    )


def load_module(module_name, args):
    try:
        mod_pkg = importlib.import_module(
            f"staphscan.modules.{module_name}.{module_name}"
        )

        if module_name == "mlst":
            return mod_pkg.Module(
                min_id=args.min_id_mlst,
                min_cov=args.min_cov_mlst,
                db_dir=args.db_mlst,
            )
        if module_name == "capsule":
            return mod_pkg.Module(
                min_id=args.min_id_capsule,
                min_cov=args.min_cov_capsule,
            )
        if module_name == "virulence":
            return mod_pkg.Module(
                min_id=args.min_id_vir,
                min_cov=args.min_cov_vir,
            )
        if module_name == "resistance":
            return mod_pkg.Module(
                min_id=args.min_id_res,
                min_cov=args.min_cov_res,
            )
        if module_name == "biofilm":
            return mod_pkg.Module(
                min_id=args.min_id_biofilm,
                min_cov=args.min_cov_biofilm,
            )

        return mod_pkg.Module()

    except Exception as e:
        sys.exit(f"Failed to import module '{module_name}': {e}")


def parse_arguments(available_modules):
    parser = argparse.ArgumentParser(
        description=(
            "Staphylococcus aureus Surveillance through Comprehensive "
            f"Analysis and staNdardized reporting (v{get_version()})"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"staphscan {get_version()}",
    )
    parser.add_argument(
        "--list-modules",
        action="store_true",
        help="List available modules and exit",
    )
    parser.add_argument(
        "--mlst_update",
        action="store_true",
        help="Authenticate and update the local PubMLST database",
    )
    parser.add_argument(
        "--db_mlst",
        type=str,
        default=None,
        help="Path to custom db folder for MLST",
    )

    io_group = parser.add_argument_group("Input/Output")
    input_group = io_group.add_mutually_exclusive_group()

    input_group.add_argument(
        "-i",
        "--input",
        nargs="+",
        help="Input FASTA files",
    )
    input_group.add_argument(
        "--r1",
        type=str,
        help="Input FASTQ reads 1",
    )
    io_group.add_argument(
        "--r2",
        type=str,
        help="Input FASTQ reads 2 (paired-end only)",
    )
    io_group.add_argument(
        "-o",
        "--outdir",
        default=None,
        help="Output directory for results",
    )

    io_group.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads used for FASTQ alignment",
    )
    io_group.add_argument(
        "--fastq_min_depth",
        type=int,
        default=3,
        help=(
            "Minimum read depth required for a reference position to be "
            "considered covered in FASTQ mode"
        ),
    )
    io_group.add_argument(
        "--fastq_min_mapq",
        type=int,
        default=0,
        help=(
            "Minimum mapping quality for FASTQ consensus generation. "
            "Each candidate is mapped on its own, so MAPQ does not "
            "reflect allele ambiguity and 0 is appropriate"
        ),
    )
    io_group.add_argument(
        "--fastq_min_baseq",
        type=int,
        default=20,
        help="Minimum base quality for FASTQ consensus generation",
    )

    mod_group = parser.add_argument_group("Modules")
    mod_group.add_argument(
        "-m",
        "--modules",
        default="all",
        help=(
            "Comma-separated list of modules to run. "
            f"Available: {', '.join(available_modules)}"
        ),
    )

    thresh_group = parser.add_argument_group("Thresholds")
    thresh_group.add_argument(
        "--min_id_mlst",
        type=float,
        default=100.0,
        help="Min identity for MLST",
    )
    thresh_group.add_argument(
        "--min_cov_mlst",
        type=float,
        default=100.0,
        help="Min coverage for MLST",
    )
    thresh_group.add_argument(
        "--min_id_capsule",
        type=float,
        default=90.0,
        help="Min identity for Capsule",
    )
    thresh_group.add_argument(
        "--min_cov_capsule",
        type=float,
        default=80.0,
        help="Min coverage for Capsule",
    )
    thresh_group.add_argument(
        "--min_id_vir",
        type=float,
        default=90.0,
        help="Min identity for Virulence",
    )
    thresh_group.add_argument(
        "--min_cov_vir",
        type=float,
        default=80.0,
        help="Min coverage for Virulence",
    )
    thresh_group.add_argument(
        "--min_id_res",
        type=float,
        default=90.0,
        help="Min identity for Resistance",
    )
    thresh_group.add_argument(
        "--min_cov_res",
        type=float,
        default=80.0,
        help="Min coverage for Resistance",
    )
    thresh_group.add_argument(
        "--min_id_biofilm",
        type=float,
        default=90.0,
        help="Min identity for Biofilm",
    )
    thresh_group.add_argument(
        "--min_cov_biofilm",
        type=float,
        default=80.0,
        help="Min coverage for Biofilm",
    )

    rep_group = parser.add_argument_group("Reporting")
    rep_group.add_argument(
        "--report",
        type=str,
        default=None,
        help="Custom filename for the report output",
    )

    args = parser.parse_args()

    if args.threads < 1:
        parser.error("--threads must be >= 1")

    if args.fastq_min_depth < 1:
        parser.error("--fastq_min_depth must be >= 1")

    if not 0 <= args.fastq_min_mapq:
        parser.error("--fastq_min_mapq must be >= 0")

    if not 0 <= args.fastq_min_baseq:
        parser.error("--fastq_min_baseq must be >= 0")

    if args.r2 and not args.r1:
        parser.error("--r2 requires --r1")

    if not (args.list_modules or args.mlst_update):
        if not args.outdir:
            parser.error("The following argument is required: -o/--outdir")
        if not args.input and not args.r1:
            parser.error(
                "One of the following arguments is required: "
                "-i/--input, --r1"
            )

    return args


def run_fasta_inputs(args, loaded_modules):
    all_results = []

    print(f"Inputs : {len(args.input)} file(s)")

    for fasta_file in args.input:
        fpath = Path(fasta_file)

        if not fpath.exists():
            print(f"Warning: File not found {fpath}")
            continue

        print(f"Processing: {fpath.stem}...")
        record = {"Sample": fpath.stem}
        perform_downstream = True

        if "assembly" in loaded_modules:
            try:
                asm_res = loaded_modules["assembly"].run(fpath)
                record.update(asm_res)

                species = asm_res.get("Species", "Unknown")
                if species != "S. aureus":
                    print(
                        f"Species identified as '{species}'. "
                        "Skipping downstream analyses."
                    )
                    perform_downstream = False

            except Exception as e:
                print(f"Error running assembly: {e}")
                record["assembly_error"] = "Fail"
                perform_downstream = False

        for name, mod in loaded_modules.items():
            if name == "assembly" or not perform_downstream:
                continue

            try:
                record.update(mod.run(fpath))
            except Exception as e:
                print(f"Error running {name}: {e}")
                record[f"{name}_error"] = "Fail"

        all_results.append(record)

    return all_results


def run_fastq_input(args, loaded_modules):
    """
    FASTQ workflow

    Each module is analysed independently:
      1. build a module-specific reference DB;
      2. map reads to the complete DB (pass 1, secondary alignments kept);
      3. identify candidate targets per family;
      4. map reads to EACH candidate separately and build a consensus
         (pass 2), so near-identical alleles cannot compete for reads;
      5. pass consensus + identity/coverage/depth metrics to run_fastq().
    """
    import staphscan.utils.fastq_mapper as fastq_mapper

    r1_path = Path(args.r1)
    r2_path = Path(args.r2) if args.r2 else None

    if not r1_path.exists():
        sys.exit(f"Input FASTQ not found: {r1_path}")

    if r2_path and not r2_path.exists():
        sys.exit(f"Input FASTQ not found: {r2_path}")

    print(f"Processing Raw Reads: {r1_path.name}...")

    sample_name = fastq_mapper.infer_sample_name(r1_path)
    record = {
        "Sample": sample_name,
        # Temporary placeholder until FASTQ species identification is
        # implemented in the complete workflow.
        "Species": "S. aureus",
    }

    with tempfile.TemporaryDirectory(prefix="staphscan_fastq_") as tmpdir:
        tmpdir = Path(tmpdir)

        for name, mod in loaded_modules.items():
            if name == "assembly":
                continue

            if not hasattr(mod, "run_fastq"):
                print(
                    f" -> Warning: Module '{name}' does not yet "
                    "support FASTQ reads."
                )
                continue

            try:
                print(f" -> FASTQ analysis: {name.upper()}")

                mod_db = tmpdir / f"{name}_refs_all.fasta"

                reference_patterns = getattr(
                    mod,
                    "fastq_reference_patterns",
                    None,
                )

                if not fastq_mapper.build_module_db(
                    mod.data_dir,
                    mod_db,
                    reference_patterns=reference_patterns,
                ):
                    print(
                        f"    No FASTQ reference sequences found for "
                        f"module '{name}'."
                    )
                    continue

                # Pass 1: broad mapping. Secondary alignments are retained
                # because they are useful for identifying near-identical
                # candidate alleles.
                bam_pass1 = tmpdir / f"{name}_pass1.bam"

                fastq_mapper.align_reads(
                    r1_path,
                    r2_path,
                    mod_db,
                    bam_pass1,
                    threads=args.threads,
                    secondary=True,
                )

                candidate_targets = fastq_mapper.select_candidate_targets(
                    bam_pass1,
                    max_candidates_per_family=3,
                    min_reads=3,
                )

                if not candidate_targets:
                    print(f"    No candidate targets detected for {name}.")
                    continue

                focused_db = tmpdir / f"{name}_refs_focused.fasta"

                if not fastq_mapper.build_focused_db(
                    mod.data_dir,
                    candidate_targets,
                    focused_db,
                    reference_patterns=reference_patterns,
                ):
                    print(
                        f"    Failed to build focused reference DB for {name}."
                    )
                    continue

                # Pass 2: map reads to each candidate separately so that
                # near-identical alleles cannot compete (MAPQ 0 discards).
                consensus_dict = fastq_mapper.consensus_per_candidate(
                    r1_path,
                    r2_path,
                    focused_db,
                    tmpdir / f"{name}_candidates",
                    threads=args.threads,
                    min_depth=args.fastq_min_depth,
                    min_mapq=args.fastq_min_mapq,
                    min_baseq=args.fastq_min_baseq,
                )

                if not consensus_dict:
                    print(
                        f"    No consensus could be built for {name}."
                    )
                    continue

                record.update(mod.run_fastq(consensus_dict))

            except Exception as e:
                print(f"Error running {name} on FASTQ: {e}")
                record[f"{name}_error"] = "Fail"

    return [record]


def write_report(args, all_results):
    if not all_results:
        sys.exit("No results generated.")

    df = pd.DataFrame(all_results).fillna("-")

    summary_cols = [
        "Sample", "Species", "Mash_distance", "Total_size", "N_contig", "N50", "QC",
        "ST", "arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL",
        "spa_type",
        "cap_type", "cap_completeness", "cap_genes",
        "sccmec_type", "sccmec_subtype", "sccmec_genes",
        "agr_type", "agr_confidence", "agr_frameshifts", "agr_operon_status",
        "res_score", "res_gene_count", "res_class_count", "Amino_res", "Bla_res",
        "Flq_res", "Gly_res", "Mec_res", "MLSB_res", "Oxa_res", "Rif_res", "Tet_res",
        "spurious_resistance_hits", "truncated_resistance_hits",
        "biofilm_score", "cna", "clfAB", "clf_genes", "fnbAB", "fnb_genes",
        "icaADBC", "ica_genes", "icaR_mutations", "biofilm_spurious_hits", "biofilm_truncated_hits",
        "vir_score", "vir_pvl", "vir_tsst", "vir_et", "vir_lukED", "vir_se",
        "spurious_virulence_hits", "truncated_virulence_hits",
    ]

    final_cols = [col for col in summary_cols if col in df.columns]

    filename = args.report or "staphscan_summary.tsv"
    if not filename.lower().endswith(".tsv"):
        filename += ".tsv"

    out_path = Path(args.outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    output_file = out_path / filename
    df[final_cols].to_csv(output_file, sep="\t", index=False)

    print(f"\nReport saved: {output_file}")
    print("Analysis complete.")


def main():
    available = get_available_modules()
    args = parse_arguments(available)

    if args.list_modules:
        print("Available StaphScan Modules:")
        for module_name in available:
            print(f"  - {module_name}")
        return

    if args.mlst_update:
        print("--- Initializing MLST Database Update ---")
        try:
            from staphscan.modules.mlst.update import run_update

            run_update(db_dir=args.db_mlst)
            return
        except Exception as e:
            sys.exit(f"Failed to run MLST update: {e}")

    if args.modules.lower() == "all":
        modules_to_run = available
    else:
        modules_to_run = [m.strip() for m in args.modules.split(",")]

        invalid = [
            module_name
            for module_name in modules_to_run
            if module_name not in available
        ]

        if invalid:
            sys.exit(
                f"Error: Unknown module(s): {', '.join(invalid)}\n"
                f"Available: {', '.join(available)}"
            )

    print("--- StaphScan Initialized ---")
    print(f"Modules: {', '.join(modules_to_run)}")

    loaded_modules = {}

    for module_name in modules_to_run:
        mod_instance = load_module(module_name, args)

        if (
            hasattr(mod_instance, "check_db")
            and not mod_instance.check_db()
        ):
            sys.exit(
                f"Error: Database check failed for module "
                f"'{module_name}'."
            )

        loaded_modules[module_name] = mod_instance

    if args.input:
        all_results = run_fasta_inputs(args, loaded_modules)
    else:
        all_results = run_fastq_input(args, loaded_modules)

    write_report(args, all_results)


if __name__ == "__main__":
    main()