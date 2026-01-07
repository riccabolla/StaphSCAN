import argparse
import sys
import pandas as pd
from pathlib import Path
import importlib

if sys.version_info < (3, 8):
    sys.exit("StaphScan requires Python 3.8+")

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
        mod_pkg = importlib.import_module(f"modules.{module_name}.{module_name}")
        
        if module_name == "mlst":
            return mod_pkg.Module(min_id=args.min_id_mlst, min_cov=args.min_cov_mlst)
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

def parse_arguments(available_modules):
    parser = argparse.ArgumentParser(
        description="StaphScan: Staphylococcus aureus Genomic Typer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--list-modules", action="store_true")

    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("-i", "--input", nargs="+")
    io_group.add_argument("-o", "--outdir")

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
    rep_group.add_argument("--complete", action="store_true")
    rep_group.add_argument("--report", type=str, default=None, help="Custom filename for the report output")

    args = parser.parse_args()

    if not args.list_modules:
        if not args.input or not args.outdir:
            parser.error("The following arguments are required: -i/--input, -o/--outdir")

    return args

def main():
    available = get_available_modules()
    args = parse_arguments(available)

    if args.list_modules:
        print("Available StaphScan Modules:")
        for m in available:
            print(f"  - {m}")
        sys.exit(0)

    if args.modules.lower() == "all":
        modules_to_run = available
    else:
        requested = [m.strip() for m in args.modules.split(',')]
        invalid = [m for m in requested if m not in available]
        if invalid:
            sys.exit(f"Error: Unknown module(s): {', '.join(invalid)}\nAvailable: {', '.join(available)}")
        modules_to_run = requested

    print(f"--- StaphScan Initialized ---")
    print(f"Modules: {', '.join(modules_to_run)}")
    print(f"Inputs : {len(args.input)} file(s)")

    loaded_modules = {}
    for m in modules_to_run:
        mod_instance = load_module(m, args)
        if hasattr(mod_instance, 'check_db') and not mod_instance.check_db():
            sys.exit(f"Error: Database check failed for module '{m}'.")
        loaded_modules[m] = mod_instance

    out_path = Path(args.outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    all_results = []

    for fasta_file in args.input:
        fpath = Path(fasta_file)
        if not fpath.exists():
            print(f"Warning: File not found {fpath}")
            continue

        print(f"Processing: {fpath.stem}...")
        record = {'Sample': fpath.stem}

        for name, mod in loaded_modules.items():
            try:
                record.update(mod.run(fpath))
            except Exception as e:
                print(f"  [!] Error running {name}: {e}")
                record[f"{name}_error"] = "Fail"

        all_results.append(record)

    if not all_results:
        sys.exit("No results generated.")

    df = pd.DataFrame(all_results).fillna("-")

    summary_cols = [
        "Sample", "Species", "Total_size", "QC", "ST", "arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL", "spa_type", #assembly module
        "cap_type", "cap_completeness", #capsule module
        "sccmec_type", #sccmec module
        "agr_type", #agr module
        "res_score", "Mec_RES", "Beta_lactamases", "Fluoroquinolones", "Tetracyclines", "Vancomycin", "Other_RES", #resistance module 
        "truncated_resistance_hits", "spurious_resistance_hits", #resistance module
        "biofilm_score", "clfAB", "clf_genes", "fnbAB", "fnb_genes", "icaADBC", "ica_genes", "icaR_mutations", "biofilm_truncated_hits", #biofilm module
        "vir_pvl", "vir_tsst", "vir_et", "spurious_virulence_hits" #virulence module
    ]

    detailed_priority = [
        "Sample", "Species", "Mash_distance", "ST", "spa_type", "spa_repeats",
        "cap_type", "cap_genes", "sccmec_type", "sccmec_genes",
        "agr_type", "Mec_RES", "Mec_AA_Found", "Mec_AA_Ref", "Beta_lactamases", "Fluoroquinolones", "Tetracyclines", "Other_RES",
        "truncated_resistance_hits", "spurious_resistance_hits",
        "biofilm_score", "biofilm_genes", "biofilm_truncated_hits", "clfAB", "clf_genes", "clfA", "clfB",
        "fnbAB", "fnb_genes","fnbA", "fnbB", "icaADBC", "ica_genes", "icaA", "icaB", "icaC", "icaD", "icaR_mutations",
        "vir_pvl", "vir_tsst", "vir_genes", "vir_spurious"
    ]

    if args.complete:
        final_cols = [c for c in detailed_priority if c in df.columns]
        remaining = [c for c in df.columns if c not in final_cols and c not in summary_cols]
        final_cols.extend(remaining)
        default_filename = "staphscan_detailed.tsv"
    else:
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
