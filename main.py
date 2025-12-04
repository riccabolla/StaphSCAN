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
    return sorted([d.name for d in modules_dir.iterdir() if d.is_dir() and (d / f"{d.name}.py").exists()])

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

    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input FASTA files")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-m", "--modules", default="all",
                        help=f"Comma-separated list of modules to run. Available: {', '.join(available_modules)}")
    parser.add_argument("--polish", action="store_true", help="Generate polished report")

    return parser.parse_args()

def main():
    available = get_available_modules()
    args = parse_arguments(available)

    modules_to_run = available if args.modules.lower() == "all" else [m.strip() for m in args.modules.split(',')]
    invalid = [m for m in modules_to_run if m not in available]
    if invalid:
        sys.exit(f"Error: Unknown module(s): {', '.join(invalid)}\nAvailable: {', '.join(available)}")

    print(f"--- StaphScan Initialized ---")
    print(f"Modules: {', '.join(modules_to_run)}")
    print(f"Inputs : {len(args.input)} file(s)")

    loaded_modules = {m: load_module(m) for m in modules_to_run}
    for name, mod in loaded_modules.items():
        if not mod.check_db():
            sys.exit(f"Error: Database check failed for module '{name}'.")

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

    polished_cols = ["Sample", "Species", "Total_size", "QC", "ST", "spa_type", "cap_type", "cap_completeness", "res_score",
                     "sccmec_type", "agr_type", "Mec_RES", "Beta_lactamases", "Fluoroquinolones", "Other_RES",
                     "spurious_resistance_hits", "vir_pvl", "vir_tsst"]

    detailed_cols = ["Sample", "Species", "Mash_distance", "ST", "spa_type", "spa_repeats", "cap_type", "cap_completeness",
                     "cap_genes", "sccmec_type", "agr_type", "res_phenotype", "res_genes", "res_score", 
                     "Drug_class", "Resistance_mechanism", "Total_size", "QC", "contig_count", "N50", 
                     "largest_contig", "Reference_accession", "Sequence_identity", "Coverage", "res_mutations", 
                     "Mec_AA_Found", "Mec_AA_Ref", "vir_pvl", "vir_tsst", "sccmec_genes", "truncated_resistance_hits"]

    if args.polish:
        final_cols = [c for c in polished_cols if c in df.columns]
        out_file = out_path / "staphscan_polished.tsv"
        print(f"Saving polished report to {out_file}")
    else:
        final_cols = [c for c in detailed_cols if c in df.columns]
        remaining = [c for c in df.columns if c not in final_cols and c not in polished_cols]
        final_cols.extend(remaining)
        out_file = out_path / "staphscan_detailed.tsv"
        print(f"Saving detailed report to {out_file}")

    df[final_cols].to_csv(out_file, sep='\t', index=False)
    print("Analysis complete.")

if __name__ == "__main__":
    main()
