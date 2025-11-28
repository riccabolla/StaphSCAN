import pandas as pd
import subprocess
import io
from pathlib import Path
from Bio.Seq import Seq

class Module:
    def __init__(self):
        self.name = "resistance"
        self.module_dir = Path(__file__).parent
        self.db_fasta = self.module_dir / "data" / "targets.fasta"
        self.blast_db = self.module_dir / "data" / "res_db"
        self.aro_csv = self.module_dir / "data" / "aro_index.csv"
        
        self.card_data = {}
        if self.aro_csv.exists():
            try:
                df_meta = pd.read_csv(self.aro_csv)
                for _, row in df_meta.iterrows():
                    acc = str(row['ARO Accession']).replace('ARO:', '')
                    self.card_data[acc] = {
                        "drug_class": row.get('Drug Class', '-'),
                        "mechanism": row.get('Resistance Mechanism', '-')
                    }
            except Exception as e:
                print(f"Warning: Could not load ARO index: {e}")

        self.mutation_targets = ["gyrA", "parC", "grlA", "gyrB", "rpoB"]
        
        self.known_mutations = {
            "gyrA": {84: ('S', ['L']), 88: ('S', ['P'])}, 
            "parC": {80: ('S', ['F', 'Y']), 84: ('E', ['K', 'G', 'V'])},
            "grlA": {80: ('S', ['F', 'Y'])},
            "rpoB": {481: ('H', ['N', 'Y']), 527: ('I', ['M'])}
        }

    def check_db(self):
        if not self.db_fasta.exists(): return False
        if not (self.module_dir / "data" / "res_db.nhr").exists():
            cmd = ["makeblastdb", "-in", str(self.db_fasta), "-dbtype", "nucl", "-out", str(self.blast_db)]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        return True

    def translate_dna(self, dna_seq, sstart, strand):
        seq_obj = Seq(dna_seq)
        if strand == "minus":
            seq_obj = seq_obj.reverse_complement()
        
        frame = (sstart - 1) % 3
        trim_amt = (3 - frame) % 3
        if trim_amt == 3: trim_amt = 0
        
        clean_seq = str(seq_obj)[trim_amt:]
        clean_seq = clean_seq[:len(clean_seq)-(len(clean_seq)%3)]
        
        if not clean_seq: return "", True

        prot_seq = str(Seq(clean_seq).translate())
        if "*" in prot_seq[:-1]:
            return prot_seq, True 
        return prot_seq.strip("*"), False

    def run(self, assembly_path):
        headers = [
            "Input_sequence_ID", "Input_gene_start", "Input_gene_stop", "Input_gene_length",
            "Reference_accession", "Sequence_identity", "Coverage",
            "Drug_class", "Resistance_mechanism", 
            "truncated_resistance_hits", "spurious_resistance_hits",
            "res_genes", "res_phenotype", "res_score", "res_mutations",
            "Mec_RES", "Beta_lactamases", "Fluoroquinolones", "Other_RES"
        ]
        defaults = {k: "-" for k in headers}
        defaults["res_score"] = "0"

        cmd = [
            "blastn", "-query", str(assembly_path), "-db", str(self.blast_db),
            "-outfmt", "6 sseqid pident length slen bitscore qseq sstart qstart qend sstrand",
            "-max_target_seqs", "200"
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if not res.stdout: return defaults

            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", 
                             names=["sseqid", "pident", "length", "slen", "bitscore", "qseq", "sstart", "qstart", "qend", "sstrand"])
            
            df = df[df['pident'] >= 40.0]
            if df.empty: return defaults

            df['family'] = df['sseqid'].apply(lambda x: x.split('_')[0])
            df = df.sort_values('bitscore', ascending=False)
            best_hits = df.drop_duplicates(subset=['family'])

            strong_genes, mutations_found = [], []
            trunc_list, spur_list, phenotypes = [], [], []
            cat_mec, cat_bla, cat_fq, cat_other = [], [], [], []

            for _, row in best_hits.iterrows():
                gene_full = row['sseqid']
                parts = gene_full.split('_')
                fam = parts[0]
                aro_id = parts[1] if len(parts) > 1 and parts[1].isdigit() else "Unknown"
                meta = self.card_data.get(aro_id, {})
                d_class = meta.get("drug_class", "-")
                
                pid = row['pident']
                cov = (row['length'] / row['slen']) * 100
                strand = row['sstrand']
                clean_name = "_".join(parts[:-1]) if len(parts) > 2 else gene_full

                if fam in self.mutation_targets:
                    prot_seq, has_stop = self.translate_dna(row['qseq'], int(row['sstart']), strand)
                    if not prot_seq: continue
                    aa_start_gap = (int(row['sstart']) - 1) // 3
                    
                    fam_muts = []
                    if fam in self.known_mutations:
                        for pos, (ref, res_list) in self.known_mutations[fam].items():
                            local_idx = pos - aa_start_gap - 1
                            if 0 <= local_idx < len(prot_seq):
                                obs_aa = prot_seq[local_idx]
                                if obs_aa != ref and obs_aa in res_list:
                                    mut_str = f"{ref}{pos}{obs_aa}"
                                    fam_muts.append(mut_str)
                                    if fam in ["gyrA", "parC", "grlA"]:
                                        phenotypes.append("Fluoroquinolones-R")
                                        cat_fq.append(f"{fam} [{mut_str}]")
                                    elif fam == "rpoB":
                                        phenotypes.append("Rifampicin-R")
                                        cat_other.append(f"{fam} [{mut_str}]")
                    if fam_muts:
                        mutations_found.append(f"{fam} [{','.join(fam_muts)}]")
                    continue

                prot_seq, has_internal_stop = self.translate_dna(row['qseq'], int(row['sstart']), strand)

                is_spurious = (pid < 90.0)
                is_truncated = (cov < 90.0)
                is_pseudogene = has_internal_stop
                
                display_str = f"{clean_name}"
                if pid < 100.0: display_str += "*" 
                
                info_tag = f"(Id:{int(pid)}% Cov:{int(cov)}%)"

                if is_pseudogene and fam in ["mecA", "mecC", "blaZ"] and cov > 90.0:
                    is_pseudogene = False
                    display_str += "(fs)"

                if is_spurious:
                    spur_list.append(f"{display_str} {info_tag}")
                    continue 
                
                if is_pseudogene:
                    trunc_list.append(f"{display_str} {info_tag} [Pseudogene]")
                    continue 

                if is_truncated:
                    trunc_list.append(f"{display_str} {info_tag} [Truncated]")
                    continue 

                strong_genes.append(display_str)

                if fam in ["mecA", "mecC"]:
                    cat_mec.append(display_str)
                    phenotypes.append("MRSA")
                elif fam == "blaZ":
                    cat_bla.append(display_str)
                    phenotypes.append("Penicillin-R")
                elif fam == "vanA":
                    cat_other.append(display_str)
                    phenotypes.append("VRSA")
                else:
                    cat_other.append(display_str)
                    if d_class and d_class != "-": phenotypes.append(d_class)

            res_score = 0
            is_mrsa = "MRSA" in phenotypes
            is_vrsa = "VRSA" in phenotypes
            is_bla = "Penicillin-R" in phenotypes
            
            if is_vrsa: res_score = 3
            elif is_mrsa: res_score = 2
            elif is_bla: res_score = 1
            
            final_out = defaults.copy()
            final_out["res_genes"] = "; ".join(strong_genes) if strong_genes else "-"
            final_out["res_mutations"] = "; ".join(mutations_found) if mutations_found else "-"
            final_out["truncated_resistance_hits"] = "; ".join(trunc_list) if trunc_list else "-"
            
            final_out["spurious_resistance_hits"] = "; ".join(spur_list) if spur_list else "-"
            
            final_out["Mec_RES"] = "; ".join(cat_mec) if cat_mec else "-"
            final_out["Beta_lactamases"] = "; ".join(cat_bla) if cat_bla else "-"
            final_out["Fluoroquinolones"] = "; ".join(cat_fq) if cat_fq else "-"
            final_out["Other_RES"] = "; ".join(cat_other) if cat_other else "-"
            
            uniq_pheno = sorted(list(set(phenotypes)))
            final_out["res_phenotype"] = " + ".join(uniq_pheno) if uniq_pheno else "Susceptible"
            final_out["res_score"] = str(res_score)

            return final_out

        except Exception as e:
            defaults["res_score"] = "Error"
            return defaults