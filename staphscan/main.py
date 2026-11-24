import argparse
import re
import sys
import tempfile
import pandas as pd
from pathlib import Path
import importlib
from importlib.metadata import version, PackageNotFoundError
from staphscan.fastq_utils import FastqProcessor

if sys.version_info < (3, 10):
    sys.exit("StaphScan requires Python 3.10+")

FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

def get_version():
    try:
        return version("staphscan")
    except PackageNotFoundError:
        return "unknown"    

def get_available_modules():
    modules_dir = Path(__file__).parent / "modules"
    if not modules_dir.exists():
        sys.exit(f"Error: 'modules' directory not found at {modules_dir}")
    return sorted([
        d.name for d in modules_dir.iterdir()
        if d.is_dir() and (d / f"{d.name}.py").exists()
    ])

def load_module(module_name, args):
    try:
        mod_pkg = importlib.import_module(f"staphscan.modules.{module_name}.{module_name}")
        
        if module_name == "mlst":
            return mod_pkg.Module(min_id=args.min_id_mlst, min_cov=args.min_cov_mlst, db_dir=args.db_mlst)
        elif module_name == "capsule":
           return mod_pkg.Module(min_id=args.min_id_capsule, min_cov=args.min_cov_capsule)
        elif module_name == "virulence":
           return mod_pkg.Module(min_id=args.min_id_vir, min_cov=args.min_cov_vir)
        elif module_name == "resistance":
            return mod_pkg.Module(min_id=args.min_id_res, min_cov=args.min_cov_res)
        elif module_name == "biofilm":
            return mod_pkg.Module(min_id=args.min_id_biofilm, min_cov=args.min_cov_biofilm)
        else:
            return mod_pkg.Module()
            
    except Exception as e:
        sys.exit(f"Failed to import module '{module_name}': {e}")

def is_fastq_file(fpath: Path) -> bool:
    name_lower = fpath.name.lower()
    return any(name_lower.endswith(ext) for ext in FASTQ_EXTENSIONS)

def clean_sample_name(fpath: Path) -> str:
    name = fpath.name
    for ext in (".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fna.gz", ".fastq", ".fq", ".fasta", ".fa", ".fna"):
        if name.lower().endswith(ext):
            name = name[:-len(ext)]
            break
    # Strip trailing _R1/_R2 or _1/_2
    name = re.sub(r"(_R?[12])(_001)?$", "", name, flags=re.IGNORECASE)
    return name

def parse_sample_inputs(args):
    """
    Organizes inputs into a list of dicts:
    [{'sample': str, 'files': [Path, ...], 'is_fastq': bool, 'read_type': str}, ...]
    """
    samples = []

    # 1. Explicit -1 / -2 short reads
    if args.read1:
        r1 = Path(args.read1)
        files = [r1]
        if args.read2:
            files.append(Path(args.read2))
        samples.append({
            "sample": clean_sample_name(r1),
            "files": files,
            "is_fastq": True,
            "read_type": args.read_type
        })

    # 2. Explicit -lr / --long_reads
    if args.long_reads:
        for lr_file in args.long_reads:
            lr_path = Path(lr_file)
            r_type = args.read_type if args.read_type in ("ont", "pb") else "ont"
            samples.append({
                "sample": clean_sample_name(lr_path),
                "files": [lr_path],
                "is_fastq": True,
                "read_type": r_type
            })

    # 3. Standard -i / --input (FASTA or auto-paired FASTQ)
    if args.input:
        fastq_by_sample = {}
        for item in args.input:
            p = Path(item)
            if is_fastq_file(p):
                sname = clean_sample_name(p)
                fastq_by_sample.setdefault(sname, []).append(p)
            else:
                samples.append({
                    "sample": p.stem,
                    "files": [p],
                    "is_fastq": False,
                    "read_type": "fasta"
                })

        for sname, fq_list in fastq_by_sample.items():
            fq_list_sorted = sorted(fq_list, key=lambda x: x.name)
            samples.append({
                "sample": sname,
                "files": fq_list_sorted[:2],
                "is_fastq": True,
                "read_type": args.read_type
            })

    return samples

def parse_arguments(available_modules):
    parser = argparse.ArgumentParser(
        description=f"Staphylococcus aureus Surveillance through Comprehensive Analysis and staNdardized reporting (v{get_version()})",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"staphscan {get_version()}"
    )
    parser.add_argument("--list-modules", action="store_true")
    parser.add_argument("--mlst_update", action="store_true", help="Authenticate and update the local PubMLST database")
    parser.add_argument("--db_mlst", type=str, default=None, help="Path to custom db folder for MLST")

    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("-i", "--input", nargs="+", help="Input FASTA or FASTQ file(s) (R1/R2 auto-paired by name)")
    io_group.add_argument("-1", "--read1", type=str, default=None, help="Forward FASTQ read (R1) or single-end short reads")
    io_group.add_argument("-2", "--read2", type=str, default=None, help="Reverse FASTQ read (R2) for paired-end reads")
    io_group.add_argument("-lr", "--long_reads", nargs="+", default=None, help="Single-end long-read FASTQ file(s) (e.g. Nanopore/PacBio)")
    io_group.add_argument("--read_type", choices=["sr", "ont", "pb"], default="sr", help="Read type for FASTQ mapping: 'sr' (Illumina), 'ont' (Nanopore), 'pb' (PacBio)")
    io_group.add_argument("--min_depth", type=int, default=5, help="Minimum read depth to call a base in FASTQ consensus")
    io_group.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for FASTQ mapping")
    io_group.add_argument("-o", "--outdir", help="Output directory")

    mod_group = parser.add_argument_group("Modules")
    mod_group.add_argument("-m", "--modules",
                           help=f"Comma-separated list of modules to run. Available: {', '.join(available_modules)}",
                           default="all")

    thresh_group = parser.add_argument_group("Thresholds")
    thresh_group.add_argument("--min_id_mlst", type=float, default=95.0, help="Min identity for MLST")
    thresh_group.add_argument("--min_cov_mlst", type=float, default=95.0, help="Min coverage for MLST")
    thresh_group.add_argument("--min_id_capsule", type=float, default=90.0, help="Min identity for Capsule")
    thresh_group.add_argument("--min_cov_capsule", type=float, default=80.0, help="Min coverage for Capsule")
    thresh_group.add_argument("--min_id_vir", type=float, default=90.0, help="Min identity for Virulence")
    thresh_group.add_argument("--min_cov_vir", type=float, default=80.0, help="Min coverage for Virulence")
    thresh_group.add_argument("--min_id_res", type=float, default=90.0, help="Min identity for Resistance")
    thresh_group.add_argument("--min_cov_res", type=float, default=80.0, help="Min coverage for Resistance")
    thresh_group.add_argument("--min_id_biofilm", type=float, default=90.0, help="Min identity for Biofilm")
    thresh_group.add_argument("--min_cov_biofilm", type=float, default=80.0, help="Min coverage for Biofilm")

    rep_group = parser.add_argument_group("Reporting")
    rep_group.add_argument("--report", type=str, default=None, help="Custom filename for the report output")
    
    args = parser.parse_args()

    has_input = any([args.input, args.read1, args.long_reads])
    if not (args.list_modules or args.mlst_update):
        if not has_input or not args.outdir:
            parser.error("Input file(s) (-i, -1/-2, or -lr) and -o/--outdir are required.")

    return args

def main():
    available = get_available_modules()
    args = parse_arguments(available)

    if args.list_modules:
        print("Available StaphScan Modules:")
        for m in available:
            print(f"  - {m}")
        sys.exit(0)

    if args.mlst_update:
        print("--- Initializing MLST Database Update ---")
        try:
            from staphscan.modules.mlst.update import run_update
            run_update(db_dir=args.db_mlst)
            sys.exit(0)    
        except Exception as e:
            sys.exit(f"Failed to run MLST update: {e}")    

    if args.modules.lower() == "all":
        modules_to_run = available
    else:
        requested = [m.strip() for m in args.modules.split(',')]
        invalid = [m for m in requested if m not in available]
        if invalid:
            sys.exit(f"Error: Unknown module(s): {', '.join(invalid)}\nAvailable: {', '.join(available)}")
        modules_to_run = requested

    samples = parse_sample_inputs(args)

    if any(s["is_fastq"] for s in samples):
        ok, err_msg = FastqProcessor.check_dependencies()
        if not ok:
            sys.exit(f"Error: {err_msg}")

    print(f"--- StaphScan Initialized ---")
    print(f"Modules: {', '.join(modules_to_run)}")
    print(f"Samples: {len(samples)}")

    loaded_modules = {}
    for m in modules_to_run:
        mod_instance = load_module(m, args)
        if hasattr(mod_instance, 'check_db') and not mod_instance.check_db():
            sys.exit(f"Error: Database check failed for module '{m}'.")
        loaded_modules[m] = mod_instance

    out_path = Path(args.outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    modules_dir = Path(__file__).parent / "modules"
    fq_processor = FastqProcessor(
        modules_dir=modules_dir,
        db_mlst=args.db_mlst,
        min_depth=args.min_depth,
        threads=args.threads
    )

    all_results = []

    for s_info in samples:
        sample_name = s_info["sample"]
        files = s_info["files"]
        is_fastq = s_info["is_fastq"]
        read_type = s_info["read_type"]

        missing = [str(p) for p in files if not p.exists()]
        if missing:
            print(f"Warning: File(s) not found for {sample_name}: {', '.join(missing)}")
            continue

        mode_label = f"FASTQ ({read_type})" if is_fastq else "FASTA"
        print(f"Processing: {sample_name} [{mode_label}]...")
        record = {'Sample': sample_name}

        perform_downstream = True

        if "assembly" in loaded_modules:
            try:
                asm_input = files if is_fastq else files[0]
                asm_res = loaded_modules["assembly"].run(asm_input, is_fastq=is_fastq, read_type=read_type)
                record.update(asm_res)
                species = asm_res.get("Species", "Unknown")
                if species != "S. aureus":
                    print(f"Species identified as '{species}'. Skipping downstream analyses.")
                    perform_downstream = False
            except Exception as e:
                print(f"Error running assembly: {e}")
                record["assembly_error"] = "Fail"
                perform_downstream = False

        downstream_modules = {k: v for k, v in loaded_modules.items() if k != "assembly"}
        if not perform_downstream or not downstream_modules:
            all_results.append(record)
            continue

        # Run downstream modules (on temporary consensus FASTA if input is FASTQ)
        with tempfile.TemporaryDirectory() as tmp_dir:
            if is_fastq:
                try:
                    target_fasta = Path(tmp_dir) / f"{sample_name}_consensus.fasta"
                    fq_processor.generate_consensus_fasta(
                        read_files=files,
                        out_fasta=target_fasta,
                        work_dir=Path(tmp_dir),
                        read_type=read_type,
                        active_modules=list(downstream_modules.keys())
                    )
                except Exception as e:
                    print(f"Error generating FASTQ consensus for {sample_name}: {e}")
                    record["fastq_mapping_error"] = "Fail"
                    all_results.append(record)
                    continue
            else:
                target_fasta = files[0]

            for name, mod in downstream_modules.items():
                try:
                    record.update(mod.run(target_fasta))
                except Exception as e:
                    print(f"Error running {name}: {e}")
                    record[f"{name}_error"] = "Fail"

        all_results.append(record)

    if not all_results:
        sys.exit("No results generated.")

    df = pd.DataFrame(all_results).fillna("-")

    summary_cols = [
        "Sample", "Species", "Mash_distance", "Total_size", "N_contig", "N50", "QC", #assembly module
        "ST", "arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL", # st module
        "spa_type", # spa module
        "cap_type", "cap_completeness", "cap_genes", #capsule module
        "sccmec_type", "sccmec_subtype", "sccmec_genes", #sccmec module
        "agr_type", "agr_confidence", "agr_frameshifts", "agr_operon_status", #agr module
        "res_score", "res_gene_count", "res_class_count", "Amino_res", "Bla_res", "Flq_res", "Gly_res", "Mec_res", "MLSB_res",  #resistance module 
        "Oxa_res", "Rif_res", "Tet_res", "spurious_resistance_hits", "truncated_resistance_hits", #resistance module
        "biofilm_score", "cna", "clfAB", "clf_genes", "fnbAB", "fnb_genes", "icaADBC", "ica_genes", "icaR_mutations", "biofilm_spurious_hits", "biofilm_truncated_hits", #biofilm module
        "vir_score", "vir_pvl", "vir_tsst", "vir_et", "vir_lukED", "vir_se", "spurious_virulence_hits", "truncated_virulence_hits" #virulence module
    ]

    final_cols = [c for c in summary_cols if c in df.columns]
    default_filename = "staphscan_summary.tsv"

    if args.report: 
        filename = args.report
        if not filename.lower().endswith('.tsv'):
            filename += '.tsv'
    else:
        filename = default_filename

    output_file = out_path / filename

    df[final_cols].to_csv(output_file, sep='\t', index=False)
    print(f"\nReport saved: {output_file}")
    print("Analysis complete.")

if __name__ == "__main__":
    main()