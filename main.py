import argparse
import sys
import os
import pandas as pd
from pathlib import Path
import importlib

if sys.version_info < (3, 8):
    sys.exit("StaphScan requires Python 3.8+")

def get_available_modules():
    modules_dir = Path(__file__).parent / "modules"
    module_names = []
    if not modules_dir.exists():
        sys.exit("Critical Error: 'modules' directory not found.")
    for item in modules_dir.iterdir():
        if item.is_dir() and not item.name.startswith('_'):
            expected_file = item / f"{item.name}.py"
            if expected_file.exists():
                module_names.append(item.name)
    return sorted(module_names)

def load_module(module_name):
    try:
        import_path = f"modules.{module_name}.{module_name}"
        mod = importlib.import_module(import_path)
        return mod.Module() 
    except Exception as e:
        sys.exit(f"Failed to import module '{module_name}': {e}")

def parse_arguments(available_modules):
    parser = argparse.ArgumentParser(
        description="StaphScan: Staphylococcus aureus Genomic Typer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("-i", "--input", nargs="+", required=True, help="Input FASTA files")
    io_group.add_argument("-o", "--outdir", required=True, help="Output directory")
    mod_group = parser.add_argument_group("Modules")
    mod_group.add_argument("-m", "--modules", 
                           help=f"Modules to run (default: all)",
                           default="all")
    return parser.parse_args()

def main():
    available_modules = get_available_modules()
    args = parse_arguments(available_modules)

    if args.modules == "all":
        modules_to_run = available_modules
    else:
        requested = args.modules.split(',')
        modules_to_run = [m for m in requested if m in available_modules]

    print(f"--- StaphScan initialized ---")
    print(f"Modules: {', '.join(modules_to_run)}")

    loaded_modules = {}
    for m in modules_to_run:
        loaded_modules[m] = load_module(m)
        if not loaded_modules[m].check_db():
            sys.exit(f"Error: Database missing for module '{m}'.")

    out_path = Path(args.outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    all_results = []
    
    for fasta_file in args.input:
        fasta_path = Path(fasta_file)
        if not fasta_path.exists(): continue

        print(f"Processing: {fasta_path.stem}...")
        record = {'Sample': fasta_path.stem}

        for mod_name, module in loaded_modules.items():
            try:
                mod_result = module.run(fasta_path)
                record.update(mod_result)
            except Exception as e:
                record[f"{mod_name}_error"] = "Fail"

        all_results.append(record)

    if not all_results:
        sys.exit("No results generated.")

    df = pd.DataFrame(all_results)
    
    col_order = [
        "Sample",
        "sccmec_type", "agr_type", "agr_match_confidence",
        "res_phenotype", "res_genes", "res_mutations", "res_spurious",
        "vir_pvl", "vir_tsst", "vir_et", "vir_spurious", 
        "sccmec_genes"
    ]
    
    final_cols = [c for c in col_order if c in df.columns]
    
    remaining = [c for c in df.columns if c not in final_cols]
    final_cols.extend(remaining)
    
    df = df[final_cols]
    
   df = df.fillna("-")

    output_file = out_path / "staphscan_results.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Done. Results: {output_file}")

if __name__ == "__main__":
    main()
