import pandas as pd
import subprocess
import io
import sys
import gzip
import shutil
import tempfile
from pathlib import Path
#from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align

class Module:
    def __init__(self):
        self.name = "agr"
        self.module_dir = Path(__file__).parent
        self.data_dir = self.module_dir / "data"
        self.db_fasta = self.data_dir / "targets.fasta"
        self.blast_db = self.data_dir / "agr_db"
        
        self.ref_data = {} 
        self.load_references()

        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.open_gap_score = -5
        self.aligner.extend_gap_score = -0.5
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1

    def check_db(self):
        if not self.db_fasta.exists(): return False
        if not (self.data_dir / "agr_db.nhr").exists(): 
            return False
        return True

    def load_references(self):
        gbk_map = {
            "gp1-operon_ref.gbk": "I",
            "gp2-operon_ref.gbk": "II",
            "gp3-operon_ref.gbk": "III",
            "gp4-operon_ref.gbk": "IV"
        }
        for fname, type_id in gbk_map.items():
            fpath = self.data_dir / fname
            if not fpath.exists(): continue
            try:
                rec = SeqIO.read(fpath, "genbank")
                full_dna = rec.seq
                proteins = {}
                for feat in rec.features:
                    if feat.type == "CDS":
                        gene = feat.qualifiers.get("gene", ["unknown"])[0]
                        translation = feat.qualifiers.get("translation", [""])[0]
                        if gene in ["agrA", "agrB", "agrC", "agrD"]:
                            proteins[gene] = translation
                self.ref_data[type_id] = {"dna": full_dna, "proteins": proteins}
            except Exception as e:
                print(f"Warning: Failed to load {fname}: {e}", file=sys.stderr)

    def extract_sequence(self, input_fasta, ref_seq):
        """Extract operon region using BLAST with padding."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_ref:
            tmp_ref.write(f">ref\n{str(ref_seq)}\n")
            tmp_ref_path = tmp_ref.name

        cmd = [
            "blastn", "-query", tmp_ref_path, "-subject", str(input_fasta),
            "-outfmt", "6 sseqid sstart send length pident bitscore",
            "-perc_identity", "75", "-qcov_hsp_perc", "40"
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            Path(tmp_ref_path).unlink() 
            if not res.stdout: return None

            df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["sseqid", "sstart", "send", "len", "pident", "score"])
            best = df.sort_values("score", ascending=False).iloc[0]
            
            seqs = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
            contig = seqs[best['sseqid']].seq
            
            pad = 1000
            start, end = best['sstart'], best['send']
            c_len = len(contig)
            
            if start < end:
                s = max(0, start - 1 - pad)
                e = min(c_len, end + pad)
                return contig[s:e]
            else:
                s = max(0, end - 1 - pad)
                e = min(c_len, start + pad)
                return contig[s:e].reverse_complement()
        except:
            if Path(tmp_ref_path).exists(): Path(tmp_ref_path).unlink()
            return None

    def analyze_gene(self, ref_prot, sample_frames): #removed gene_name and sample_len_bp

        best_score = -1.0
        best_seq = ""
        #best_frame_idx = -1
        ref_len = len(ref_prot)
        
        for prot in sample_frames: #prot in enumerate(sample_frames):
            if not prot: continue
            score = self.aligner.score(ref_prot, prot)
            if score > best_score:
                best_score = score
                best_seq = prot
                #best_frame_idx = i

        if best_score <= 0: return "Missing"

        alignments = self.aligner.align(ref_prot, best_seq)
        if not alignments: return "Missing"
        
        aln = alignments[0]
        
        covered_bases = 0
        for start, end in aln.aligned[0]:
            covered_bases += (end - start)
        coverage_pct = (covered_bases / ref_len) * 100
        
        if coverage_pct < 75.0:
            q_start = aln.aligned[1][0][0] 
            q_end = aln.aligned[1][-1][1]  
            
            seq_len_aa = len(best_seq)
            
            # assembly gap check
            dist_start = q_start
            dist_end = seq_len_aa - q_end
            
            if dist_start < 5 or dist_end < 5:

                return "Incomplete (Contig Edge)"
                
            return "Truncated/Frameshift"
            
        q_start = aln.aligned[1][0][0]
        q_end = aln.aligned[1][-1][1]
        aligned_segment = best_seq[q_start:q_end]
        
        if "*" in aligned_segment[:-1]: 
            return "Premature Stop"
            
        return "Intact"

    def check_frameshifts(self, sample_dna, type_id):
        ref_prots = self.ref_data[type_id]["proteins"]
        mutated_genes = []
        targets = ["agrC", "agrA"] 
        
        def safe_frames(seq_obj):
            f = []
            for i in range(3):
                s = seq_obj[i:]
                valid = (len(s) // 3) * 3
                if valid > 0: f.append(str(s[:valid].translate(table=11)))
                else: f.append("")
            return f

        frames = safe_frames(sample_dna) + safe_frames(sample_dna.reverse_complement())
        #sample_len_bp = len(sample_dna)

        for gene in targets:
            if gene not in ref_prots: continue
            
            status = self.analyze_gene(ref_prots[gene], frames) #removed gene_name and sample_len_bp
            
            if status != "Intact":
                mutated_genes.append(f"{gene} ({status})")
                
        if not mutated_genes:
            return "No"
        else:
            return f"Yes ({', '.join(mutated_genes)})"

    def run(self, assembly_path):
        out = {
            "agr_type": "Negative", 
            "agr_confidence": "0",  
            "agr_frameshifts": "-",       
            "agr_operon_status": "-"      
        }

        if not assembly_path.exists(): return out
        
        with tempfile.TemporaryDirectory() as td:
            input_path = assembly_path
            if assembly_path.suffix == ".gz":
                input_path = Path(td) / "unzipped.fasta"
                try:
                    with gzip.open(assembly_path, 'rb') as f_in:
                        with open(input_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                except: return out
            
            cmd = [
                "blastn", "-task", "blastn-short", 
                "-query", str(input_path), "-db", str(self.blast_db),
                "-outfmt", "6 sseqid pident length slen bitscore", 
                "-perc_identity", "90",   
                "-max_target_seqs", "5000"
            ]
            
            try:
                res = subprocess.run(cmd, capture_output=True, text=True)
                if not res.stdout: return out
                
                df = pd.read_csv(io.StringIO(res.stdout), sep="\t", names=["sseqid", "pident", "length", "slen", "bitscore"])
                
                df['coverage'] = (df['length'] / df['slen']) * 100
                df = df[df['pident'] >= 90.0] 
                
                if df.empty: return out

                df = df.sort_values(by=['pident', 'coverage', 'bitscore'], ascending=[False, False, False])
                best_hit = df.iloc[0]
                
                raw_id = str(best_hit['sseqid'])
                best_group_tag = raw_id.split('_')[0]
                
                type_map = {'gp1': 'I', 'gp2': 'II', 'gp3': 'III', 'gp4': 'IV'}
                final_type = type_map.get(best_group_tag, best_group_tag)
                
                out["agr_type"] = f"agr {final_type}"
                
                df['group'] = df['sseqid'].apply(lambda x: str(x).split('_')[0])
                probes_found = df[df['group'] == best_group_tag]['sseqid'].nunique()
                out["agr_confidence"] = str(probes_found) #confidence score, number of probes matching

                type_key = type_map.get(best_group_tag, None)
                
                if type_key and type_key in self.ref_data:
                    ref_seq = self.ref_data[type_key]["dna"]
                    extracted_operon = self.extract_sequence(input_path, ref_seq)
                    
                    if extracted_operon:
                        frameshift_result = self.check_frameshifts(extracted_operon, type_key)
                        out["agr_frameshifts"] = frameshift_result #frameshift detection
                        
                        if frameshift_result == "No":
                            out["agr_operon_status"] = "Intact"
                        elif "Incomplete" in frameshift_result:
                            out["agr_operon_status"] = "Assembly Gap" # reduce false positives
                        else:
                            out["agr_operon_status"] = "Pseudogene"
                    else:
                        out["agr_operon_status"] = "Missing/Fragmented"
                else:
                    out["agr_operon_status"] = "Ref Missing"

                return out

            except Exception as e:
                print(f"Error in agr module: {e}", file=sys.stderr)
                return out