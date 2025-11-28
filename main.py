import argparse
import sys
import pandas as pd
from pathlib import Path
import importlib

if sys.version_info < (3, 8): sys.exit("StaphScan requires Python 3.8+")

def get_available_modules():
    modules_dir = Path(__file__).parent / "modules"
    return sorted([d.name for d in modules_dir.iterdir() if d.is_dir() and (d/f"{d.name}.py").exists()])

def load_module(module_name):
    try:
        mod = importlib.import_module(f"modules.{module_name}.{module_name}")
        return mod.Module() 
    except Exception as e:
        sys.exit(f"Failed to import module '{module_name}': {e}")

def main():
    parser = argparse.ArgumentParser(description="StaphScan: S. aureus Genomic Typer")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input FASTA files")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--polish", action="store_true", help="Generate simplified report")
    args = parser.parse_args()

    modules = get_available_modules()
    print(f"--- StaphScan initialized ---")
    print(f"Loaded Modules: {', '.join(modules)}")

    loaded_mods = {m: load_module(m) for m in modules}
    
    for m, mod in loaded_mods.items():
        if not mod.check_db():
            sys.exit(f"Error: Database check failed for module '{m}'.")

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    all_results = []

    for fasta_file in args.input:
        fpath = Path(fasta_file)
        if not fpath.exists(): continue
        
        print(f"Processing: {fpath.stem}...")
        record = {'Sample': fpath.stem}
        
        for name, mod in loaded_mods.items():
            try:
                record.update(mod.run(fpath))
            except Exception as e:
                record[f"{name}_error"] = "Fail"
        
        all_results.append(record)

    if not all_results: sys.exit("No results.")
    df = pd.DataFrame(all_results).fillna("-")

    polished_cols = [
        "Sample",
        "Species",
        "Total_size",
        "QC",
        "ST",
        "res_score",
        "sccmec_type",
        "agr_type",
        "Mec_RES",
        "Beta_lactamases",
        "Fluoroquinolones",
        "Other_RES",
        "spurious_resistance_hits",
        "vir_pvl",
        "vir_tsst",
        "vir_spurious" 
    ]

    detailed_cols = [
        "Sample",
        "Species", 
        "Mash_distance",
        "ST",
        "arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL",
        "sccmec_type", 
        "agr_type", 
        "res_phenotype", 
        "res_genes",
        "res_score",
        "Drug_class",
        "Resistance_mechanism",
        "Total_size",
        "QC",
        "contig_count",
        "N50",
        "largest_contig",
        "ambiguous_bases",
        "Reference_accession",
        "Sequence_identity",
        "Coverage",
        "res_mutations",
        "vir_pvl", "vir_tsst", "vir_et", "vir_spurious",
        "sccmec_genes",
        "truncated_resistance_hits", 
        "spurious_resistance_hits"
    ]

    if args.polish:
        final_cols = [c for c in polished_cols if c in df.columns]
        out_file = Path(args.outdir) / "staphscan_polished.tsv"
        df[final_cols].to_csv(out_file, sep='\t', index=False)
        print(f"Done: {out_file}")

    else:
        final_cols = [c for c in detailed_cols if c in df.columns]
        remaining = [c for c in df.columns if c not in final_cols and c not in polished_cols]
        final_cols.extend(remaining)
        
        out_file = Path(args.outdir) / "staphscan_detailed.tsv"
        df[final_cols].to_csv(out_file, sep='\t', index=False)
        print(f"Done: {out_file}")

if __name__ == "__main__":
    main()