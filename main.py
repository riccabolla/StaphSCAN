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

def load_module(module_name):
    try:
        mod = importlib.import_module(f"modules.{module_name}.{module_name}")
        return mod.Module()
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

    rep_group = parser.add_argument_group("Reporting")
    rep_group.add_argument("--complete", action="store_true")

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
        mod_instance = load_module(m)
        if not mod_instance.check_db():
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
        "Sample", "Species", "Total_size", "QC", "ST", "spa_type",
        "cap_type", "cap_completeness", "sccmec_type", "agr_type",
        "res_score", "Mec_RES", "Beta_lactamases", "Fluoroquinolones", "Other_RES",
        "biofilm_score", "clfAB", "clf_genes", "fnbAB", "fnb_genes", "icaADBC", "ica_genes", "icaR_mutations",
        "vir_pvl", "vir_tsst"
    ]

    detailed_priority = [
        "Sample", "Species", "Mash_distance", "ST", "spa_type", "spa_repeats",
        "cap_type", "cap_genes", "sccmec_type", "sccmec_genes",
        "agr_type",
        "res_phenotype", "res_genes", "res_score", "res_mutations",
        "Mec_AA_Found", "Mec_AA_Ref", "truncated_resistance_hits", "spurious_resistance_hits",
        "biofilm_score", "biofilm_genes", "biofilm_truncated_hits",
        "clfAB", "clf_genes", "clfA", "clfB",
        "fnbAB", "fnb_genes","fnbA", "fnbB",
        "icaADBC", "ica_genes", "icaA", "icaB", "icaC", "icaD",
        "icaR_mutations",
        "vir_pvl", "vir_tsst", "vir_genes", "vir_spurious"
    ]

    if args.complete:
        final_detailed_cols = [c for c in detailed_priority if c in df.columns]
        remaining = [c for c in df.columns if c not in final_detailed_cols and c not in summary_cols]
        final_detailed_cols.extend(remaining)
        detailed_file = out_path / "staphscan_detailed.tsv"
        df[final_detailed_cols].to_csv(detailed_file, sep='\t', index=False)
        print(f"\n Complete report saved: {detailed_file}")
    else:
        final_summary_cols = [c for c in summary_cols if c in df.columns]
        summary_file = out_path / "staphscan_summary.tsv"
        df[final_summary_cols].to_csv(summary_file, sep='\t', index=False)
        print(f"\n Summary report saved: {summary_file}")

    print("Analysis complete.")

if __name__ == "__main__":
    main()
