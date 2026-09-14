import argparse
import sys
import pandas as pd
from pathlib import Path
import importlib
from importlib.metadata import version, PackageNotFoundError

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

    parser.add_argument("--mlst_update", action="store_true", help="Authenticate and update the local PubMLST database") #mlst update option

    parser.add_argument("--db_mlst", type=str, default=None, help="Path to custom db folder for MLST")

    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument("-i", "--input", nargs="+")
    io_group.add_argument("--r1", type=str, help="Input reads 1")
    io_group.add_argument("--r2", type=str, help="Input reads 2 (not required if single-end or long reads)")
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
    #rep_group.add_argument("--complete", action="store_true") # removed in this version, maybe later
    rep_group.add_argument("--report", type=str, default=None, help="Custom filename for the report output")
    
    args = parser.parse_args()

    if not (args.list_modules or args.mlst_update):
        if not args.outdir:
            parser.error("The following arguments are required: -o/--outdir")
        if not args.input and not args.r1:
            parser.error("One of the following arguments is required: -i/--input, --r1")
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
            # imported only if called
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

    print(f"--- StaphScan Initialized ---")
    print(f"Modules: {', '.join(modules_to_run)}")

    loaded_modules = {}
    for m in modules_to_run:
        mod_instance = load_module(m, args)
        if hasattr(mod_instance, 'check_db') and not mod_instance.check_db():
            sys.exit(f"Error: Database check failed for module '{m}'.")
        loaded_modules[m] = mod_instance

    out_path = Path(args.outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    all_results = []

    if args.input:
        print(f"Inputs : {len(args.input)} file(s)")
        for fasta_file in args.input:
            fpath = Path(fasta_file)
            if not fpath.exists():
                print(f"Warning: File not found {fpath}")
                continue

            print(f"Processing: {fpath.stem}...")
            record = {'Sample': fpath.stem}

            perform_downstream = True

            if "assembly" in loaded_modules:
                try: 
                    asm_res = loaded_modules["assembly"].run(fpath)
                    record.update(asm_res)
                    species = asm_res.get("Species", "Unknown")
                    if species != "S. aureus":
                        print(f"Species identified as '{species}'. Skipping downstream analyses.") #species filtering
                        perform_downstream = False
                except Exception as e:
                    print(f"Error running assembly: {e}")
                    record["assembly_error"] = "Fail"
                    perform_downstream = False

            for name, mod in loaded_modules.items():
                if name == "assembly":
                    continue
                if not perform_downstream:
                    continue
                try:
                    record.update(mod.run(fpath))
                except Exception as e:
                    print(f"Error running {name}: {e}")
                    record[f"{name}_error"] = "Fail"

            all_results.append(record)

    elif args.r1:
        r1_path = Path(args.r1)
        r2_path = Path(args.r2) if args.r2 else None

        print(f"Processing Raw Reads: {r1_path.stem}...")
        sample_name = r1_path.stem.replace('.fastq', '').replace('.fq', '').replace('_R1', '').replace('_1', '')
        record = {'Sample': sample_name}
        record["Species"] = "S. aureus" # Placeholder until Mash is integrated

        import tempfile

        import staphscan.utils.fastq_mapper
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for name, mod in loaded_modules.items():
                if name == "assembly":
                    continue
                try:
                    if hasattr(mod, 'run_fastq'):
                        print(f" -> Mapping reads independently for {name.upper()} module...")
                        mod_db = Path(tmpdir) / f"{name}_refs_all.fasta"

                        if not staphscan.utils.fastq_mapper.build_module_db(mod.data_dir, mod_db):
                            continue

                        # Pass 1: map against ALL alleles to let minimap2's own scoring
                        # pick the best-supported allele per gene family.
                        bam_pass1 = Path(tmpdir) / f"{name}_pass1.bam"
                        staphscan.utils.fastq_mapper.align_reads(r1_path, r2_path, mod_db, bam_pass1)
                        best_targets = staphscan.utils.fastq_mapper.select_best_targets(bam_pass1)

                        # Pass 2: remap against only the winning allele per family — removes
                        #   the multi-allele read-splitting that was capping breadth.
                        focused_db = Path(tmpdir) / f"{name}_refs_focused.fasta"
                        if not staphscan.utils.fastq_mapper.build_focused_db(mod.data_dir, best_targets, focused_db):
                            continue

                        bam_out = Path(tmpdir) / f"{name}.bam"
                        staphscan.utils.fastq_mapper.align_reads(r1_path, r2_path, focused_db, bam_out)

                        module_min_cov = getattr(mod, "min_cov", 80.0)
                        consensus_dict = staphscan.utils.fastq_mapper.extract_consensus_from_bam(bam_out, min_depth=1, min_breadth=module_min_cov)

                        record.update(mod.run_fastq(consensus_dict))
                    else:
                        print(f" -> Warning: Module '{name}' does not yet support FASTQ reads.")
                except Exception as e:
                    print(f"Error running {name} on FASTQ: {e}")
                    record[f"{name}_error"] = "Fail"
                    
        all_results.append(record)

    if not all_results:
        sys.exit("No results generated.")

    df = pd.DataFrame(all_results).fillna("-")

    summary_cols = [
        "Sample", "Species", "Mash_distance","Total_size", "N_contig", "N50", "QC", #assembly module
        "ST", "arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL", # st module
        "spa_type", # spa module
        "cap_type", "cap_completeness", "cap_genes", #capsule module
        "sccmec_type", "sccmec_subtype", "sccmec_genes", #sccmec module
        "agr_type", "agr_confidence", "agr_frameshifts", "agr_operon_status", #agr module
        "res_score", "res_gene_count", "res_class_count", "Amino_res", "Bla_res", "Flq_res","Gly_res","Mec_res", "MLSB_res",  #resistance module 
        "Oxa_res", "Rif_res", "Tet_res", "spurious_resistance_hits", "truncated_resistance_hits", #resistance module
        "biofilm_score", "cna","clfAB", "clf_genes", "fnbAB", "fnb_genes", "icaADBC", "ica_genes", "icaR_mutations", "biofilm_spurious_hits", "biofilm_truncated_hits", #biofilm module
        "vir_score","vir_pvl", "vir_tsst", "vir_et", "vir_lukED","vir_se", "spurious_virulence_hits", "truncated_virulence_hits" #virulence module
    ]
    
    # currently removed the detailed report option
    # detailed_priority = [
    #     "Sample", "Species", "Mash_distance", "ST", "spa_type", "spa_repeats",
    #     "cap_type", "cap_genes", "sccmec_type", "sccmec_genes",
    #     "agr_type", "Mec_RES", "Mec_AA_Found", "Mec_AA_Ref", "Beta_lactamases", "Fluoroquinolones", "Tetracyclines", "Other_RES",
    #     "truncated_resistance_hits", "spurious_resistance_hits",
    #     "biofilm_score", "biofilm_genes", "biofilm_truncated_hits", "clfAB", "clf_genes", "clfA", "clfB",
    #     "fnbAB", "fnb_genes","fnbA", "fnbB", "icaADBC", "ica_genes", "icaA", "icaB", "icaC", "icaD", "icaR_mutations",
    #     "vir_pvl", "vir_tsst", "vir_genes", "vir_spurious"
    # ]

    # if args.complete:
    #     final_cols = [c for c in detailed_priority if c in df.columns]
    #     remaining = [c for c in df.columns if c not in final_cols and c not in summary_cols]
    #     final_cols.extend(remaining)
    #     default_filename = "staphscan_detailed.tsv"
    # else:
    #     final_cols = [c for c in summary_cols if c in df.columns]
    #     default_filename = "staphscan_summary.tsv"
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
