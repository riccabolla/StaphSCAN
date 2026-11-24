import gzip
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from Bio import SeqIO
from Bio.Seq import Seq

_SAM_COMPLEMENT = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")
_CIGAR_OP_RE = re.compile(rb"(\d+)([MIDNSHP=X])")

# cna_M81736 (3,552 bp) single-copy A-domain and 3x561-bp tandem B-repeat coordinates (0-based, half-open)
_CNA_A_START, _CNA_A_END = 100, 1550
_CNA_B_START, _CNA_B_END = 1596, 3279
_CNA_REF_B_COPIES = 3

_MSCR_BIO_FAMS = ("clfa", "clfb", "fnba", "fnbb")
_MSCR_BIO_PREFIXES = (b"bio__clfa__", b"bio__clfb__", b"bio__fnba__", b"bio__fnbb__")


class FastqProcessor:
    def __init__(self, modules_dir: Path, db_mlst: Optional[str] = None, min_depth: int = 5, threads: int = 4):
        self.modules_dir = Path(modules_dir)
        self.db_mlst = Path(db_mlst) / "mlst" if db_mlst else self.modules_dir / "mlst" / "data"
        self.min_depth = int(min_depth)
        self.threads = int(threads)

    @staticmethod
    def check_dependencies() -> Tuple[bool, str]:
        for tool in ("minimap2", "samtools"):
            if shutil.which(tool) is None:
                return False, f"Required external tool '{tool}' not found in PATH."
        return True, ""

    @staticmethod
    def _classify_sccmec_target(raw_id: str) -> str:
        lower_id = raw_id.lower()
        for tag in (
            "ccra1", "ccrb1", "ccra2", "ccrb2", "ccra3", "ccrb3",
            "ccra4", "ccrb4", "ccrb6", "ccrc1", "ccrc2", "ccrc",
            "meci", "mecr1", "is1272", "is431", "is12960d",
            "meca", "mecc", "blaz"
        ):
            if tag in lower_id:
                return "ccrc1" if tag == "ccrc" else tag
        return lower_id.split("_")[0]

    def build_combined_reference(self, dest_fasta: Path, active_modules: Optional[List[str]] = None, read_type: str = "sr") -> Dict[str, str]:
        """
        Builds a non-redundant target reference FASTA (1 representative per gene family,
        except divergent fnbA/fnbB allelic types which only share ~35% k-mers) so primary
        alignments receive full read depth and MAPQ > 0 without multi-allele splitting.
        """
        if active_modules is None:
            active_modules = ["mlst", "capsule", "sccmec", "agr", "biofilm", "virulence", "resistance"]

        seen_families = set()
        seen_headers = set()
        seen_seqs = set()
        ref_dict: Dict[str, str] = {}

        def write_unique(out_handle, raw_header: str, seq_str: str):
            if seq_str in seen_seqs:
                return
            seen_seqs.add(seq_str)
            clean_id = raw_header.split()[0]
            unique_id = clean_id
            counter = 1
            while unique_id in seen_headers:
                counter += 1
                unique_id = f"{clean_id}_dup{counter}"
            seen_headers.add(unique_id)
            out_handle.write(f">{unique_id}\n{seq_str}\n")
            ref_dict[unique_id] = seq_str

        with open(dest_fasta, "w") as out_f:
            # 1. MLST representative alleles (1 allele per locus)
            if "mlst" in active_modules:
                mlst_refs = self.db_mlst / "refs.fasta"
                if mlst_refs.exists():
                    for rec in SeqIO.parse(mlst_refs, "fasta"):
                        write_unique(out_f, f"mlst__{rec.id}", str(rec.seq).upper())

            # 2. Capsule targets (1 copy of shared core genes A-G, L-P + both cap5H-K and cap8H-K)
            if "capsule" in active_modules:
                cap_files = list((self.modules_dir / "capsule" / "data").glob("capsule_targets*.fasta"))
                if cap_files:
                    seen_cap_core = set()
                    for rec in SeqIO.parse(cap_files[0], "fasta"):
                        m = re.search(r"cap([58])([A-P])", rec.id)
                        if m:
                            _, locus = m.group(1), m.group(2)
                            if locus not in ("H", "I", "J", "K"):
                                if locus in seen_cap_core:
                                    continue
                                seen_cap_core.add(locus)
                        write_unique(out_f, f"cap__{rec.id}", str(rec.seq).upper())

            # 3. AGR operons
            if "agr" in active_modules:
                agr_data = self.modules_dir / "agr" / "data"
                gbk_files = sorted(agr_data.glob("*.gbk"))
                if gbk_files:
                    for gbk in gbk_files:
                        for rec in SeqIO.parse(gbk, "genbank"):
                            write_unique(out_f, f"agrop__{gbk.stem}_{rec.id}", str(rec.seq).upper())
                else:
                    for fpath in agr_data.glob("targets*.fasta"):
                        for rec in SeqIO.parse(fpath, "fasta"):
                            write_unique(out_f, f"agr__{rec.id}", str(rec.seq).upper())

            # 4. Biofilm, Virulence, and Resistance (1 representative per gene family;
            # keep distinct fnbA_1/fnbA_2 and fnbB_1/fnbB_2 types because they diverge by >60%)
            rep_sources = []
            if "biofilm" in active_modules:
                bio_files = list((self.modules_dir / "biofilm" / "data").glob("biofilm*.fasta"))
                if bio_files:
                    rep_sources.append(("bio", bio_files[0]))

            if "virulence" in active_modules:
                vir_files = list((self.modules_dir / "virulence" / "data").glob("targets*.fasta"))
                if vir_files:
                    rep_sources.append(("vir", vir_files[0]))

            if "resistance" in active_modules:
                res_data = self.modules_dir / "resistance" / "data"
                for res_f in sorted(res_data.glob("*_res*.fasta")):
                    rep_sources.append(("res", res_f))

            for prefix, fpath in rep_sources:
                for rec in SeqIO.parse(fpath, "fasta"):
                    gene_family = rec.id.split("_")[0].lower()
                    seq_str = str(rec.seq).upper()
                    if prefix == "bio" and gene_family in _MSCR_BIO_FAMS and read_type == "sr":
                        # Include all distinct alleles of hypervariable MSCRAMM genes so
                        # minimap2 baits complete paired reads for local de Bruijn assembly.
                        write_unique(out_f, f"{prefix}__{gene_family}__{rec.id}", seq_str)
                        continue
                    if prefix == "bio" and gene_family in ("fnba", "fnbb"):
                        # Keep fnbA_1 + fnbA_2 (which covers fnbA_3/4) and fnbB_1 + fnbB_2
                        sub_key = "_".join(rec.id.split("_")[:2]).lower()
                        if sub_key in ("fnba_3", "fnba_4"):
                            continue
                        if sub_key in seen_families:
                            continue
                        seen_families.add(sub_key)
                    else:
                        if gene_family in seen_families:
                            continue
                        seen_families.add(gene_family)
                    write_unique(out_f, f"{prefix}__{gene_family}__{rec.id}", seq_str)

            # 5. SCCmec targets (1 representative per target gene)
            if "sccmec" in active_modules:
                scc_data = self.modules_dir / "sccmec" / "data"
                t_files = list(scc_data.glob("targets*.fasta"))
                if t_files:
                    for rec in SeqIO.parse(t_files[0], "fasta"):
                        tag = self._classify_sccmec_target(rec.id)
                        if tag in seen_families:
                            continue
                        seen_families.add(tag)
                        write_unique(out_f, f"scctarg__{tag}__{rec.id}", str(rec.seq).upper())

        return ref_dict

    @staticmethod
    def _detect_pileup_frameshift(seq_dmin_ins: str, seq_dmin_noins: str, ref_seq: str) -> Optional[str]:
        """
        Checks for a genuine high-depth (>= min_depth) frameshift truncation across the full
        gene without truncating at the first N. Replaces N dips with ref_seq bases while
        retaining all high-depth insertions and deletions.
        """
        if not seq_dmin_ins or not ref_seq or len(seq_dmin_noins) != len(ref_seq):
            return None

        # Only test when high-depth reads contain an actual indel relative to the reference
        has_deletion = "*" in seq_dmin_ins
        has_insertion = len(seq_dmin_ins) != len(seq_dmin_noins)
        if not (has_deletion or has_insertion):
            return None

        # Walk seq_dmin_noins (1:1 with ref_seq) and seq_dmin_ins in lockstep:
        # fill 'N' positions from ref_seq, delete '*' positions, and keep inserted bases
        reconstructed: List[str] = []
        ref_idx = 0
        ref_len = len(ref_seq)

        for ch in seq_dmin_ins:
            if ref_idx < ref_len and ch == seq_dmin_noins[ref_idx]:
                if ch == "N":
                    reconstructed.append(ref_seq[ref_idx])
                elif ch == "*":
                    pass  # true high-depth deletion: omit base to preserve frameshift
                else:
                    reconstructed.append(ch)
                ref_idx += 1
            else:
                # Inserted base present in seq_dmin_ins but not in 1:1 seq_dmin_noins
                if ch not in ("N", "*"):
                    reconstructed.append(ch)

        candidate_dna = "".join(reconstructed)
        if len(candidate_dna) < 300 or abs(len(candidate_dna) - ref_len) % 3 == 0:
            return None

        trim = len(candidate_dna) - (len(candidate_dna) % 3)
        prot = str(Seq(candidate_dna[:trim]).translate(table=11))
        stop_pos = prot.find("*")
        ref_aa_len = ref_len // 3
        if stop_pos != -1 and stop_pos < (ref_aa_len * 0.90):
            return candidate_dna
        return None

    @staticmethod
    def _strip_deletions_in_frame(seq_with_stars: str, ref_seq: str) -> str:
        """
        Removes in-frame whole-codon deletions ('*' runs where length % 3 == 0) and restores
        isolated out-of-frame '*' characters from ref_seq when no true frameshift stop was detected.
        """
        if "*" not in seq_with_stars or len(seq_with_stars) != len(ref_seq):
            return seq_with_stars.replace("*", "")

        pieces: List[str] = []
        last_idx = 0
        for m in re.finditer(r"\*+", seq_with_stars):
            s, e = m.start(), m.end()
            del_len = e - s
            pieces.append(seq_with_stars[last_idx:s])
            if del_len % 3 != 0:
                pieces.append(ref_seq[s:e])
            last_idx = e
        pieces.append(seq_with_stars[last_idx:])
        return "".join(pieces)

    @staticmethod
    def _resolve_paired_repeats(
        unitigs: List[str],
        read_pairs: List[Tuple[str, str]],
        k: int,
        min_unitig_len: int = 350
    ) -> List[str]:
        """
        Resolves unambiguous 2-in, 2-out repeat unitigs (e.g. shared internal domains
        between tandem fnbA/fnbB paralogs) using paired-end read links.
        """
        if not unitigs or not read_pairs:
            return sorted([u for u in unitigs if len(u) >= min_unitig_len], key=len, reverse=True)

        rc_trans = str.maketrans("ACGT", "TGCA")

        def rc(s: str) -> str:
            return s.translate(rc_trans)[::-1]

        def rev_node(n: Tuple[int, int]) -> Tuple[int, int]:
            return (n[0], 1 - n[1])

        k_ov = k - 1
        k_idx = 31
        seqs: Dict[Tuple[int, int], str] = {}
        for uid, u in enumerate(unitigs):
            seqs[(uid, 0)] = u
            seqs[(uid, 1)] = rc(u)

        by_pref: Dict[str, List[Tuple[int, int]]] = {}
        for n, s in seqs.items():
            if len(s) >= k_ov:
                by_pref.setdefault(s[:k_ov], []).append(n)

        out_adj: Dict[Tuple[int, int], List[Tuple[int, int]]] = {n: [] for n in seqs}
        in_adj: Dict[Tuple[int, int], List[Tuple[int, int]]] = {n: [] for n in seqs}
        for n, s in seqs.items():
            if len(s) < k_ov:
                continue
            for nxt in by_pref.get(s[-k_ov:], ()):
                if nxt == n or nxt == rev_node(n):
                    continue
                out_adj[n].append(nxt)
                in_adj[nxt].append(n)

        resolved_repeats: Set[int] = set()
        bridge_map: Dict[Tuple[Tuple[int, int], Tuple[int, int]], Tuple[int, int]] = {}

        def match_seq(seq: str, tbl: Dict[str, Tuple[int, int]]) -> Optional[Tuple[int, int]]:
            for off in range(0, len(seq) - k_idx + 1, 15):
                hit = tbl.get(seq[off:off + k_idx])
                if hit is not None:
                    return hit
            return None

        for uid, u in enumerate(unitigs):
            if len(u) > 500 or len(u) < k:
                continue
            rnode = (uid, 0)
            ins = in_adj[rnode]
            outs = out_adj[rnode]
            if len(ins) != 2 or len(outs) != 2 or len(set(ins)) != 2 or len(set(outs)) != 2:
                continue
            if rnode in ins or rnode in outs or rev_node(rnode) in ins or rev_node(rnode) in outs:
                continue

            in_km: Dict[str, Set[Tuple[int, int]]] = {}
            for n in ins:
                s = seqs[n]
                end_w = max(0, len(s) - k_ov)
                win = s[max(0, end_w - 650):end_w]
                for i in range(len(win) - k_idx + 1):
                    in_km.setdefault(win[i:i + k_idx], set()).add(n)

            out_km: Dict[str, Set[Tuple[int, int]]] = {}
            for n in outs:
                s = seqs[n]
                start_w = min(len(s), k_ov)
                win = s[start_w:min(len(s), start_w + 650)]
                for i in range(len(win) - k_idx + 1):
                    out_km.setdefault(win[i:i + k_idx], set()).add(n)

            in_uniq = {km: next(iter(v)) for km, v in in_km.items() if len(v) == 1}
            out_uniq = {km: next(iter(v)) for km, v in out_km.items() if len(v) == 1}
            if not in_uniq or not out_uniq:
                continue

            votes = {(i, o): 0 for i in ins for o in outs}
            for a, b in read_pairs:
                for fwd, rev in ((a, rc(b)), (b, rc(a))):
                    hi = match_seq(fwd, in_uniq)
                    if hi is not None:
                        ho = match_seq(rev, out_uniq)
                        if ho is not None:
                            votes[(hi, ho)] += 1

            i1, i2 = ins
            o1, o2 = outs
            s_dir = votes[(i1, o1)] + votes[(i2, o2)]
            s_crs = votes[(i1, o2)] + votes[(i2, o1)]
            min_dir = min(votes[(i1, o1)], votes[(i2, o2)])
            min_crs = min(votes[(i1, o2)], votes[(i2, o1)])

            chosen = None
            if s_dir >= 4 and min_dir >= 2 and s_dir >= 4 * (s_crs + 1):
                chosen = ((i1, o1), (i2, o2))
            elif s_crs >= 4 and min_crs >= 2 and s_crs >= 4 * (s_dir + 1):
                chosen = ((i1, o2), (i2, o1))

            if chosen is not None:
                resolved_repeats.add(uid)
                for ia, oa in chosen:
                    bridge_map[(ia, rnode)] = oa
                    bridge_map[(rev_node(oa), rev_node(rnode))] = rev_node(ia)

        if not resolved_repeats:
            return sorted([u for u in unitigs if len(u) >= min_unitig_len], key=len, reverse=True)

        next_link: Dict[Tuple[int, int], Tuple[Tuple[int, int], Tuple[int, int]]] = {}
        has_prev: Set[Tuple[int, int]] = set()
        for (ia, rn), oa in bridge_map.items():
            if len(out_adj[ia]) == 1 and len(in_adj[oa]) == 1:
                next_link[ia] = (rn, oa)
                has_prev.add(oa)

        merged_seqs: List[str] = []
        used_uids: Set[int] = set(resolved_repeats)

        for n in sorted(seqs):
            if n[0] in resolved_repeats or n in has_prev or n not in next_link:
                continue
            chain_parts = [seqs[n]]
            cur = n
            seen_chain = {cur}
            used_uids.add(cur[0])
            while cur in next_link:
                rn, nxt = next_link[cur]
                if nxt in seen_chain:
                    break
                seen_chain.add(nxt)
                used_uids.add(nxt[0])
                chain_parts.append(seqs[rn][k_ov:])
                chain_parts.append(seqs[nxt][k_ov:])
                cur = nxt
            merged_seqs.append("".join(chain_parts))

        for uid, u in enumerate(unitigs):
            if uid not in used_uids:
                merged_seqs.append(u)

        seen_canon: Set[str] = set()
        final_unitigs: List[str] = []
        for s in sorted(merged_seqs, key=len, reverse=True):
            if len(s) < min_unitig_len:
                continue
            sr = rc(s)
            if s in seen_canon or sr in seen_canon:
                continue
            seen_canon.add(s)
            seen_canon.add(sr)
            final_unitigs.append(s)
        return final_unitigs

    @staticmethod
    def _assemble_local_debruijn(
        reads: List[str],
        min_depth: int = 5,
        min_unitig_len: int = 350,
        read_pairs: Optional[List[Tuple[str, str]]] = None
    ) -> List[str]:
        """
        Iterative multi-k de Bruijn unitig assembler for baited hypervariable
        MSCRAMM reads (clfA, clfB, fnbA, fnbB) to avoid reference-template bias.
        """
        clean_reads = [r.strip().upper() for r in reads if len(r.strip()) >= 55]
        if not clean_reads:
            return []

        rc_trans = str.maketrans("ACGT", "TGCA")

        def rc(s: str) -> str:
            return s.translate(rc_trans)[::-1]

        def build_graph(valid_kmers: Dict[str, int]) -> Tuple[Dict[str, List[str]], Dict[str, int]]:
            out_e: Dict[str, List[str]] = {}
            in_d: Dict[str, int] = {}
            for km in valid_kmers:
                u, b, v = km[:-1], km[-1], km[1:]
                out_e.setdefault(u, []).append(b)
                in_d[v] = in_d.get(v, 0) + 1
            return out_e, in_d

        max_rlen = max(len(r) for r in clean_reads)
        if max_rlen >= 135:
            ks = (55, 77, 101)
        elif max_rlen >= 100:
            ks = (41, 61, 81)
        else:
            ks = (31, 45, 55)

        prev_unitigs: List[str] = []
        final_deduped: List[str] = []
        final_k: int = ks[-1]

        for stage_idx, k in enumerate(ks):
            min_kcov = max(3, min_depth) if stage_idx == 0 else max(2, min_depth - 2)
            kmer_counts: Dict[str, int] = {}

            for r in clean_reads:
                if len(r) < k:
                    continue
                for seq in (r, rc(r)):
                    for i in range(len(seq) - k + 1):
                        km = seq[i:i + k]
                        if not km.strip("ACGT"):
                            kmer_counts[km] = kmer_counts.get(km, 0) + 1

            for u in prev_unitigs:
                if len(u) < k:
                    continue
                for seq in (u, rc(u)):
                    for i in range(len(seq) - k + 1):
                        km = seq[i:i + k]
                        if not km.strip("ACGT"):
                            kmer_counts[km] = kmer_counts.get(km, 0) + min_kcov

            valid = {km: c for km, c in kmer_counts.items() if c >= min_kcov}
            if not valid:
                break

            out_edges, in_deg = build_graph(valid)

            for _ in range(2):
                changed = False

                # Clean 1: Low-coverage error edge pruning
                to_remove: Set[str] = set()
                for u, bs in out_edges.items():
                    if len(bs) > 1:
                        max_c = max(valid[u + b] for b in bs)
                        cutoff = max_c * 0.20
                        for b in bs:
                            e = u + b
                            if valid[e] < cutoff:
                                to_remove.add(e)
                                to_remove.add(rc(e))
                if to_remove:
                    for e in to_remove:
                        valid.pop(e, None)
                    out_edges, in_deg = build_graph(valid)
                    changed = True

                # Clean 2: Dead-end tip clipping (2 passes)
                for _tip_pass in range(2):
                    to_remove = set()
                    for u, bs in out_edges.items():
                        if len(bs) <= 1:
                            continue
                        branches = []
                        for b in bs:
                            e = u + b
                            cov0 = valid[e]
                            curr = e[1:]
                            path = [e]
                            steps = 1
                            while steps < 2 * k and in_deg.get(curr, 0) == 1 and len(out_edges.get(curr, ())) == 1:
                                nb = out_edges[curr][0]
                                ne = curr + nb
                                if ne in path:
                                    break
                                path.append(ne)
                                curr = ne[1:]
                                steps += 1
                            is_dead = len(out_edges.get(curr, ())) == 0 and steps < 2 * k
                            branches.append((cov0, is_dead, path))
                        for cov0, is_dead, path in branches:
                            if is_dead and any((not odead) or (ocov0 > cov0) for ocov0, odead, _ in branches):
                                for pe in path:
                                    to_remove.add(pe)
                                    to_remove.add(rc(pe))
                    if not to_remove:
                        break
                    for e in to_remove:
                        valid.pop(e, None)
                    out_edges, in_deg = build_graph(valid)
                    changed = True

                # Clean 3: Simple bubble popping
                to_remove = set()
                max_bubble_steps = 2 * k + 15
                for u, bs in out_edges.items():
                    if len(bs) <= 1:
                        continue
                    by_dest: Dict[str, List[Tuple[float, int, int, List[str]]]] = {}
                    for b in bs:
                        e = u + b
                        if e not in valid:
                            continue
                        curr = e[1:]
                        path = [e]
                        steps = 1
                        while steps < max_bubble_steps and in_deg.get(curr, 0) == 1 and len(out_edges.get(curr, ())) == 1:
                            nb = out_edges[curr][0]
                            ne = curr + nb
                            if ne in path:
                                break
                            path.append(ne)
                            curr = ne[1:]
                            steps += 1
                        if curr != u and (in_deg.get(curr, 0) > 1 or len(out_edges.get(curr, ())) != 1):
                            covs = [valid[pe] for pe in path]
                            avg_cov = sum(covs) / len(covs)
                            min_cov_b = min(covs)
                            by_dest.setdefault(curr, []).append((avg_cov, min_cov_b, len(path), path))
                    for dest_branches in by_dest.values():
                        if len(dest_branches) >= 2:
                            dest_branches.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
                            for _, _, _, drop_path in dest_branches[1:]:
                                for pe in drop_path:
                                    to_remove.add(pe)
                                    to_remove.add(rc(pe))
                if to_remove:
                    for e in to_remove:
                        valid.pop(e, None)
                    out_edges, in_deg = build_graph(valid)
                    changed = True

                if not changed:
                    break

            # Clean 4 (final stage only): Resolve unambiguous 1-in, 1-out simple tandem repeat loops (<= 350 bp)
            if stage_idx == len(ks) - 1:
                for _loop_pass in range(2):
                    to_remove = set()
                    for u, bs in out_edges.items():
                        if len(bs) != 2:
                            continue
                        for b in bs:
                            e = u + b
                            if e not in valid or e in to_remove:
                                continue
                            curr = e[1:]
                            loop_path = [e]
                            steps = 1
                            while steps < 350 and in_deg.get(curr, 0) == 1 and len(out_edges.get(curr, ())) == 1:
                                nb = out_edges[curr][0]
                                ne = curr + nb
                                if ne in loop_path:
                                    break
                                loop_path.append(ne)
                                curr = ne[1:]
                                steps += 1
                            w = curr
                            is_tandem_loop = (w == u and in_deg.get(u, 0) == 2)
                            if not is_tandem_loop and w != u and in_deg.get(w, 0) == 2 and len(out_edges.get(w, ())) == 1:
                                r_curr = (w + out_edges[w][0])[1:]
                                r_steps = 1
                                seen_r = {w}
                                while r_steps < 350 and r_curr != u and in_deg.get(r_curr, 0) == 1 and len(out_edges.get(r_curr, ())) == 1:
                                    if r_curr in seen_r:
                                        break
                                    seen_r.add(r_curr)
                                    r_curr = (r_curr + out_edges[r_curr][0])[1:]
                                    r_steps += 1
                                if r_curr == u and in_deg.get(u, 0) == 1:
                                    is_tandem_loop = True
                            if is_tandem_loop:
                                for pe in loop_path:
                                    to_remove.add(pe)
                                    to_remove.add(rc(pe))
                                break
                    if not to_remove:
                        break
                    for e in to_remove:
                        valid.pop(e, None)
                    out_edges, in_deg = build_graph(valid)

            # Unitig compaction
            visited_edges: Set[str] = set()
            stage_unitigs: List[str] = []

            def walk_unitig(start_u: str, start_b: str) -> None:
                e = start_u + start_b
                if e in visited_edges:
                    return
                visited_edges.add(e)
                bases = [e]
                curr = e[1:]
                while in_deg.get(curr, 0) == 1 and len(out_edges.get(curr, ())) == 1:
                    nb = out_edges[curr][0]
                    ne = curr + nb
                    if ne in visited_edges:
                        break
                    visited_edges.add(ne)
                    bases.append(nb)
                    curr = ne[1:]
                stage_unitigs.append("".join(bases))

            for u in sorted(out_edges):
                bs = out_edges[u]
                if in_deg.get(u, 0) != 1 or len(bs) != 1:
                    for b in sorted(bs):
                        walk_unitig(u, b)

            for u in sorted(out_edges):
                for b in sorted(out_edges[u]):
                    if (u + b) not in visited_edges:
                        walk_unitig(u, b)

            seen_canon: Set[str] = set()
            deduped: List[str] = []
            for seq in sorted(stage_unitigs, key=len, reverse=True):
                seq_rc = rc(seq)
                if seq in seen_canon or seq_rc in seen_canon:
                    continue
                seen_canon.add(seq)
                seen_canon.add(seq_rc)
                deduped.append(seq)

            final_deduped = deduped
            final_k = k
            prev_unitigs = [u for u in deduped if len(u) >= k + 20]

        if read_pairs and final_deduped:
            return FastqProcessor._resolve_paired_repeats(
                final_deduped, read_pairs, k=final_k, min_unitig_len=min_unitig_len
            )
        return sorted([u for u in prev_unitigs if len(u) >= min_unitig_len], key=len, reverse=True)

    def _extend_mscr_reads(
        self,
        seed_unitigs: List[str],
        read_files: List[Path],
        work_dir: Path,
        baited_qnames: Set[bytes],
        tag_prefix: str = "main"
    ) -> Tuple[List[str], List[Tuple[str, str]]]:
        """
        Recruits additional read pairs overlapping the tips of Round-1 MSCRAMM
        unitigs so divergent internal domains are bridged before final assembly.
        """
        seed_fasta = work_dir / f"{tag_prefix}_mscr_seed.fasta"
        with open(seed_fasta, "w") as fh:
            for idx, u_seq in enumerate(seed_unitigs, 1):
                fh.write(f">seed_{idx}\n{u_seq}\n")

        mm2_cmd = [
            "minimap2", "-ax", "sr",
            "-t", str(self.threads),
            "--secondary=no",
            str(seed_fasta)
        ] + [str(f) for f in read_files]

        p_mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        extra_reads: List[str] = []
        extra_pairs: List[Tuple[str, str]] = []
        curr_qname: Optional[bytes] = None
        curr_pair_seqs: List[bytes] = []
        curr_hit = False

        try:
            for raw_line in p_mm2.stdout:
                if raw_line.startswith(b"@"):
                    continue
                parts = raw_line.rstrip(b"\r\n").split(b"\t")
                if len(parts) < 11:
                    continue
                flag = int(parts[1])
                qname = parts[0]

                if qname != curr_qname:
                    if curr_hit and curr_qname is not None and curr_qname not in baited_qnames:
                        baited_qnames.add(curr_qname)
                        pair_strs = [s.decode("ascii", errors="ignore").upper() for s in curr_pair_seqs]
                        extra_reads.extend(pair_strs)
                        if len(pair_strs) == 2:
                            extra_pairs.append((pair_strs[0], pair_strs[1]))
                    curr_qname = qname
                    curr_pair_seqs = []
                    curr_hit = False

                if not (flag & 2304) and parts[9] != b"*":
                    seq_b = parts[9].translate(_SAM_COMPLEMENT)[::-1] if (flag & 16) else parts[9]
                    curr_pair_seqs.append(seq_b)

                if not (flag & 4):
                    curr_hit = True

            if curr_hit and curr_qname is not None and curr_qname not in baited_qnames:
                baited_qnames.add(curr_qname)
                pair_strs = [s.decode("ascii", errors="ignore").upper() for s in curr_pair_seqs]
                extra_reads.extend(pair_strs)
                if len(pair_strs) == 2:
                    extra_pairs.append((pair_strs[0], pair_strs[1]))
        finally:
            if p_mm2.stdout:
                try:
                    p_mm2.stdout.close()
                except Exception:
                    pass
                p_mm2.stdout = None
            p_mm2.communicate()

        return extra_reads, extra_pairs

    def _map_and_consensus(
        self,
        ref_fasta: Path,
        ref_dict: Dict[str, str],
        read_files: List[Path],
        work_dir: Path,
        read_type: str,
        tag_prefix: str = "main"
    ) -> List[Tuple[str, str]]:
        if not ref_dict:
            return []

        preset_map = {"sr": "sr", "ont": "map-ont", "pb": "map-hifi"}
        preset = preset_map.get(read_type, "sr")
        bam_path = work_dir / f"{tag_prefix}.sorted.bam"
        cons_d1_noins_path = work_dir / f"{tag_prefix}_d1_noins.fasta"
        cons_dmin_ins_path = work_dir / f"{tag_prefix}_dmin_ins.fasta"
        cons_dmin_path = work_dir / f"{tag_prefix}_dmin.fasta"

        if read_type == "sr":
            mm2_cmd = [
                "minimap2", "-ax", preset,
                "-t", str(self.threads),
                "-k", "13", "-w", "5", "-f", "10000",
                "-Y", "-N", "100", "-p", "0.70",
                str(ref_fasta)
            ] + [str(f) for f in read_files]
        else:
            mm2_cmd = [
                "minimap2", "-ax", preset,
                "-t", str(self.threads),
                "-Y", "-N", "25", "-p", "0.70",
                str(ref_fasta)
            ] + [str(f) for f in read_files]

        sort_cmd = ["samtools", "sort", "-@", str(self.threads), "-o", str(bam_path), "-"]

        p_mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p_sort = subprocess.Popen(sort_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        primary_cache: Dict[Tuple[bytes, int], Tuple[int, bytes, bytes, Set[bytes]]] = {}
        curr_qname: Optional[bytes] = None
        curr_pair_seqs: List[bytes] = []
        curr_is_mscr = False
        mscr_baited_reads: List[str] = []
        mscr_baited_pairs: List[Tuple[str, str]] = []
        baited_qnames: Set[bytes] = set()
        cna_a_bases = 0
        cna_b_bases = 0

        try:
            for raw_line in p_mm2.stdout:
                if raw_line.startswith(b"@"):
                    p_sort.stdin.write(raw_line)
                    continue

                line_clean = raw_line.rstrip(b"\r\n")
                parts = line_clean.split(b"\t")
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                qname = parts[0]

                if qname != curr_qname:
                    if curr_is_mscr and read_type == "sr":
                        if curr_qname is not None:
                            baited_qnames.add(curr_qname)
                        pair_strs = [s.decode("ascii", errors="ignore").upper() for s in curr_pair_seqs]
                        mscr_baited_reads.extend(pair_strs)
                        if len(pair_strs) == 2:
                            mscr_baited_pairs.append((pair_strs[0], pair_strs[1]))
                    curr_qname = qname
                    curr_pair_seqs = []
                    curr_is_mscr = False

                if read_type == "sr" and not (flag & 2304) and parts[9] != b"*":
                    seq_b = parts[9].translate(_SAM_COMPLEMENT)[::-1] if (flag & 16) else parts[9]
                    curr_pair_seqs.append(seq_b)

                if flag & 4:
                    continue

                mate_bits = flag & 192
                is_rev = flag & 16
                rname = parts[2]
                key = (qname, mate_bits)

                if read_type == "sr" and rname.startswith(_MSCR_BIO_PREFIXES):
                    curr_is_mscr = True
                    continue

                if not (flag & 2304):  # primary alignment
                    seq_bytes = parts[9]
                    if seq_bytes != b"*":
                        if len(primary_cache) > 8192:
                            primary_cache.clear()
                        primary_cache[key] = (is_rev, seq_bytes, parts[10], {rname})

                    if rname.startswith(b"bio__cna__") and parts[5] != b"*":
                        ref_pos = int(parts[3]) - 1
                        for n_len_b, op in _CIGAR_OP_RE.findall(parts[5]):
                            n_len = int(n_len_b)
                            if op in (b"M", b"=", b"X"):
                                ref_end = ref_pos + n_len
                                if ref_end > _CNA_A_START and ref_pos < _CNA_A_END:
                                    cna_a_bases += min(ref_end, _CNA_A_END) - max(ref_pos, _CNA_A_START)
                                if ref_end > _CNA_B_START and ref_pos < _CNA_B_END:
                                    cna_b_bases += min(ref_end, _CNA_B_END) - max(ref_pos, _CNA_B_START)
                                ref_pos = ref_end
                            elif op in (b"D", b"N"):
                                ref_pos += n_len

                    # For non-sr MSCRAMM paralogs, require MAPQ > 0 on primary alignments
                    mapq = int(parts[4])
                    if rname.startswith(_MSCR_BIO_PREFIXES) and mapq == 0:
                        continue

                    p_sort.stdin.write(raw_line)
                else:
                    # Never copy secondary alignments onto biofilm MSCRAMM genes (cna, clfA, clfB, fnbA, fnbB)
                    if rname.startswith(b"bio__"):
                        continue

                    cached = primary_cache.get(key)
                    if cached is None:
                        continue
                    prim_rev, prim_seq, prim_qual, seen_rnames = cached
                    if rname in seen_rnames:
                        continue
                    seen_rnames.add(rname)

                    if parts[9] == b"*":
                        if is_rev == prim_rev:
                            parts[9] = prim_seq
                            parts[10] = prim_qual
                        else:
                            parts[9] = prim_seq.translate(_SAM_COMPLEMENT)[::-1]
                            parts[10] = prim_qual[::-1]
                    if flag & 256:
                        flag &= ~256
                        parts[1] = str(flag).encode("ascii")
                    p_sort.stdin.write(b"\t".join(parts) + b"\n")

            if curr_is_mscr and read_type == "sr":
                if curr_qname is not None:
                    baited_qnames.add(curr_qname)
                pair_strs = [s.decode("ascii", errors="ignore").upper() for s in curr_pair_seqs]
                mscr_baited_reads.extend(pair_strs)
                if len(pair_strs) == 2:
                    mscr_baited_pairs.append((pair_strs[0], pair_strs[1]))
        except BrokenPipeError:
            pass
        finally:
            if p_sort.stdin:
                try:
                    p_sort.stdin.close()
                except Exception:
                    pass
                p_sort.stdin = None
            if p_mm2.stdout:
                try:
                    p_mm2.stdout.close()
                except Exception:
                    pass
                p_mm2.stdout = None

        _, sort_err = p_sort.communicate()
        _, mm2_err = p_mm2.communicate()

        if p_mm2.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {mm2_err.decode('utf-8', errors='replace').strip()}")
        if p_sort.returncode != 0 or not bam_path.exists():
            raise RuntimeError(f"samtools sort failed: {sort_err.decode('utf-8', errors='replace').strip()}")

        # 1. 1:1 coordinate consensus at depth >= 1 (used only to check 5'/3' read-overhang tips)
        cmd_d1_noins = [
            "samtools", "consensus", "-m", "simple", "-c", "0.51",
            "-f", "fasta", "-a", "--show-ins", "no", "--show-del", "yes",
            "-d", "1", "-o", str(cons_d1_noins_path), str(bam_path)
        ]
        subprocess.run(cmd_d1_noins, capture_output=True, text=True, check=True)

        # 2. High-depth insertion-aware consensus at depth >= min_depth (--show-ins yes)
        cmd_dmin_ins = [
            "samtools", "consensus", "-m", "simple", "-c", "0.51",
            "-f", "fasta", "-a", "--show-ins", "yes", "--show-del", "yes",
            "-d", str(self.min_depth), "-o", str(cons_dmin_ins_path), str(bam_path)
        ]
        subprocess.run(cmd_dmin_ins, capture_output=True, text=True, check=True)

        # 3. High-depth 1:1 coordinate consensus at depth >= min_depth (--show-ins no)
        cmd_dmin = [
            "samtools", "consensus", "-m", "simple", "-c", "0.51",
            "-f", "fasta", "-a", "--show-ins", "no", "--show-del", "yes",
            "-d", str(self.min_depth), "-o", str(cons_dmin_path), str(bam_path)
        ]
        subprocess.run(cmd_dmin, capture_output=True, text=True, check=True)

        d1_noins_seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(cons_d1_noins_path, "fasta")}
        dmin_ins_seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(cons_dmin_ins_path, "fasta")}
        dmin_seqs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(cons_dmin_path, "fasta")}

        best_by_group: Dict[str, Tuple[Tuple[float, float, int], str, str]] = {}

        for seq_id, seq_dmin in dmin_seqs.items():
            total_len = len(seq_dmin)
            if total_len == 0:
                continue
            seq_d1 = d1_noins_seqs.get(seq_id, seq_dmin)
            seq_dmin_ins = dmin_ins_seqs.get(seq_id, seq_dmin)
            ref_seq = ref_dict.get(seq_id, "")

            non_n_dmin = total_len - seq_dmin.count("N")
            non_n_d1 = len(seq_d1) - seq_d1.count("N")
            cov_dmin = non_n_dmin / total_len
            cov_d1 = non_n_d1 / len(seq_d1) if seq_d1 else 0.0

            if cov_dmin < 0.25 or cov_d1 < 0.35:
                continue

            id_parts = seq_id.split("__")
            gene_fam = id_parts[1] if len(id_parts) >= 3 else ""

            # Check for a genuine high-depth frameshift truncation across the full pileup (e.g. cna-25% in 03)
            fs_dna = self._detect_pileup_frameshift(seq_dmin_ins, seq_dmin, ref_seq) if cov_dmin >= 0.80 else None
            if fs_dna is not None:
                healed_seq = fs_dna
            else:
                seq_list = list(seq_dmin)

                if ref_seq and len(ref_seq) == total_len:
                    # Step A: Estimate B-repeat copy number in cna_M81736 from primary aligned depth
                    # relative to the single-copy A-domain. Alleles with < 2.5 estimated B copies
                    # (1- or 2-repeat variants) have single-HSP coverage < 80% and are masked in B2/B3
                    # so blastn reflects true < 80% coverage rather than reference-biased healing.
                    if gene_fam == "cna" and total_len == 3552:
                        mean_a_depth = cna_a_bases / (_CNA_A_END - _CNA_A_START)
                        mean_b_depth = cna_b_bases / (_CNA_B_END - _CNA_B_START)
                        if mean_a_depth >= self.min_depth:
                            est_b_copies = (_CNA_REF_B_COPIES * mean_b_depth) / mean_a_depth
                            if est_b_copies < 2.5:
                                seq_list[2157:3279] = ["N"] * (3279 - 2157)
                                cov_dmin = (total_len - seq_list.count("N")) / total_len

                    # Step B: Restore 5'/3' terminal bases that are covered at d>=1 (seq_d1) when the gene
                    # has strong overall high-depth coverage (cov_dmin >= 0.85), preventing taper-off at contig ends
                    if cov_dmin >= 0.85 and len(seq_d1) == total_len:
                        tip_limit = 45
                        for idx in range(min(tip_limit, total_len)):
                            if seq_list[idx] == "N" and seq_d1[idx] != "N":
                                seq_list[idx] = seq_d1[idx]
                        for idx in range(max(0, total_len - tip_limit), total_len):
                            if seq_list[idx] == "N" and seq_d1[idx] != "N":
                                seq_list[idx] = seq_d1[idx]

                    # Step C: Strictly capped reference healing (<= 10% of total gene length)
                    # ONLY when unhealed cov_dmin >= 0.83 (so <80% spurious hits are NEVER healed into Complete calls)
                    if cov_dmin >= 0.83 and gene_fam != "cna" and id_parts[0] != "agrop":
                        max_heal_budget = int(total_len * 0.10)
                        healed_so_far = 0
                        cur_seq = "".join(seq_list)
                        for m in re.finditer(r"N+", cur_seq):
                            s, e = m.start(), m.end()
                            n_len = e - s
                            if healed_so_far + n_len > max_heal_budget:
                                continue
                            # Heal short 5'/3' overhang tips (<= 25 bp) or internal MAPQ=0 repeat dropouts
                            if (s <= 15 or e >= total_len - 15) and n_len <= 25:
                                seq_list[s:e] = list(ref_seq[s:e])
                                healed_so_far += n_len
                            elif s > 0 and e < int(total_len * 0.95):
                                seq_list[s:e] = list(ref_seq[s:e])
                                healed_so_far += n_len

                healed_seq = self._strip_deletions_in_frame("".join(seq_list), ref_seq)

            if len(id_parts) >= 3 and id_parts[0] in ("bio", "vir", "res", "scctarg"):
                group_key = f"{id_parts[0]}__{id_parts[1]}"
            elif id_parts[0] == "agrop":
                group_key = "agrop"
            else:
                group_key = seq_id

            rank_score = (round(cov_dmin, 4), round(cov_d1, 4), non_n_dmin)
            if group_key not in best_by_group or rank_score > best_by_group[group_key][0]:
                best_by_group[group_key] = (rank_score, seq_id, healed_seq)

        records = [(seq_id, healed_seq) for _, seq_id, healed_seq in best_by_group.values()]
        if read_type == "sr" and mscr_baited_reads:
            seed_unitigs = self._assemble_local_debruijn(
                mscr_baited_reads,
                min_depth=self.min_depth,
                min_unitig_len=200,
                read_pairs=mscr_baited_pairs
            )
            if seed_unitigs:
                extra_reads, extra_pairs = self._extend_mscr_reads(
                    seed_unitigs=seed_unitigs,
                    read_files=read_files,
                    work_dir=work_dir,
                    baited_qnames=baited_qnames,
                    tag_prefix=tag_prefix
                )
                if extra_reads:
                    mscr_baited_reads.extend(extra_reads)
                    mscr_baited_pairs.extend(extra_pairs)
                    unitigs = self._assemble_local_debruijn(
                        mscr_baited_reads,
                        min_depth=self.min_depth,
                        min_unitig_len=350,
                        read_pairs=mscr_baited_pairs
                    )
                else:
                    unitigs = [u for u in seed_unitigs if len(u) >= 350]
            else:
                unitigs = []
            for idx, u_seq in enumerate(unitigs, 1):
                records.append((f"bio_local__unitig_{idx}", u_seq))
        return records

    def _map_sccmec_regions(
        self,
        read_files: List[Path],
        work_dir: Path,
        read_type: str
    ) -> Optional[Tuple[str, str]]:
        scc_data = self.modules_dir / "sccmec" / "data"
        r_files = list(scc_data.glob("regions*.fasta"))
        if not r_files:
            return None

        reg_ref = r_files[0]
        preset_map = {"sr": "sr", "ont": "map-ont", "pb": "map-hifi"}
        preset = preset_map.get(read_type, "sr")
        bam_path = work_dir / "sccreg.sorted.bam"
        cons_path = work_dir / "sccreg_dmin.fasta"

        mm2_cmd = [
            "minimap2", "-ax", preset,
            "-t", str(self.threads),
            "-Y", "-N", "25", "-p", "0.80",
            str(reg_ref)
        ] + [str(f) for f in read_files]

        sort_cmd = ["samtools", "sort", "-@", str(self.threads), "-o", str(bam_path), "-"]

        p_mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p_sort = subprocess.Popen(sort_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        primary_cache: Dict[Tuple[bytes, int], Tuple[int, bytes, bytes]] = {}
        try:
            for raw_line in p_mm2.stdout:
                if raw_line.startswith(b"@"):
                    p_sort.stdin.write(raw_line)
                    continue
                line_clean = raw_line.rstrip(b"\r\n")
                parts = line_clean.split(b"\t")
                if len(parts) < 11:
                    continue
                flag = int(parts[1])
                if flag & 4:
                    continue
                key = (parts[0], flag & 192)
                is_rev = flag & 16
                if not (flag & 2304):
                    if parts[9] != b"*":
                        if len(primary_cache) > 8192:
                            primary_cache.clear()
                        primary_cache[key] = (is_rev, parts[9], parts[10])
                    p_sort.stdin.write(raw_line)
                else:
                    cached = primary_cache.get(key)
                    if cached is None:
                        continue
                    prim_rev, prim_seq, prim_qual = cached
                    if parts[9] == b"*":
                        if is_rev == prim_rev:
                            parts[9] = prim_seq
                            parts[10] = prim_qual
                        else:
                            parts[9] = prim_seq.translate(_SAM_COMPLEMENT)[::-1]
                            parts[10] = prim_qual[::-1]
                    if flag & 256:
                        flag &= ~256
                        parts[1] = str(flag).encode("ascii")
                    p_sort.stdin.write(b"\t".join(parts) + b"\n")
        except BrokenPipeError:
            pass
        finally:
            if p_sort.stdin:
                try:
                    p_sort.stdin.close()
                except Exception:
                    pass
                p_sort.stdin = None
            if p_mm2.stdout:
                try:
                    p_mm2.stdout.close()
                except Exception:
                    pass
                p_mm2.stdout = None

        p_sort.communicate()
        p_mm2.communicate()
        if not bam_path.exists():
            return None

        cmd_cons = [
            "samtools", "consensus", "-m", "simple", "-c", "0.51",
            "-f", "fasta", "-a", "--show-ins", "no", "--show-del", "yes",
            "-d", str(self.min_depth), "-o", str(cons_path), str(bam_path)
        ]
        res = subprocess.run(cmd_cons, capture_output=True, text=True)
        if res.returncode != 0 or not cons_path.exists():
            return None

        best_score = None
        best_record = None
        for rec in SeqIO.parse(cons_path, "fasta"):
            seq_str = str(rec.seq).upper()
            total_len = len(seq_str)
            if total_len == 0:
                continue
            non_n = total_len - seq_str.count("N")
            cov = non_n / total_len
            if cov < 0.68:
                continue
            score = (round(cov, 4), non_n)
            if best_score is None or score > best_score:
                best_score = score
                best_record = (f"sccreg__{rec.id}", seq_str.replace("*", ""))

        return best_record

    def generate_consensus_fasta(
        self,
        read_files: List[Path],
        out_fasta: Path,
        work_dir: Path,
        read_type: str = "sr",
        active_modules: Optional[List[str]] = None
    ) -> Path:
        if active_modules is None:
            active_modules = ["mlst", "spa", "capsule", "sccmec", "agr", "biofilm", "virulence", "resistance"]

        ref_fasta = work_dir / "combined_refs.fasta"
        ref_dict = self.build_combined_reference(ref_fasta, active_modules=active_modules, read_type=read_type)

        consensus_records = self._map_and_consensus(
            ref_fasta=ref_fasta,
            ref_dict=ref_dict,
            read_files=read_files,
            work_dir=work_dir,
            read_type=read_type,
            tag_prefix="targets"
        )

        with open(out_fasta, "w") as final_f:
            for seq_id, seq_str in consensus_records:
                final_f.write(f">{seq_id}\n{seq_str}\n")

            if "sccmec" in active_modules:
                best_reg = self._map_sccmec_regions(read_files, work_dir=work_dir, read_type=read_type)
                if best_reg:
                    final_f.write(f">{best_reg[0]}\n{best_reg[1]}\n")

            if "spa" in active_modules:
                spa_contig = self.reconstruct_spa_locus(read_files, read_type=read_type)
                if spa_contig:
                    final_f.write(f">spa_reconstructed_amplicon\n{spa_contig}\n")

        return out_fasta

    def _iter_fastq_seqs(self, read_files: List[Path]):
        for fpath in read_files:
            open_fn = gzip.open if fpath.name.lower().endswith(".gz") else open
            with open_fn(fpath, "rt", errors="ignore") as fh:
                for line_idx, line in enumerate(fh):
                    if line_idx % 4 == 1:
                        yield line.strip().upper()

    def reconstruct_spa_locus(self, read_files: List[Path], read_type: str = "sr") -> Optional[str]:
        primer_sets = [
            ('TAAAGACGATCCTTCGGTGAG', 'CAGCAGTAGTGCCGTTTGCTT'),        # Set 1
            ('AGACGATCCTTCGGTGAGC', 'GCTTTTGCAATGTCATTTACTG'),         # Set 2
            ('ATAGCGTGATTTTGCGGTT', 'CTAAATATAAATAATGTTGTCACTTGGA'),   # Set 3
            ('CAACGCAATGGTTTCATCCA', 'GCTTTTGCAATGTCATTTACTG')         # Set 4
        ]

        all_seeds = set()

        def add_sliding_kmers(s: str, k_seed: int = 12):
            s_rc = str(Seq(s).reverse_complement())
            for target in (s, s_rc):
                for idx in range(len(target) - k_seed + 1):
                    all_seeds.add(target[idx:idx + k_seed])

        spa_data_dir = self.modules_dir / "spa" / "data"
        repeat_files = list(spa_data_dir.glob("repeats*.fasta"))
        if repeat_files:
            for rec in SeqIO.parse(repeat_files[0], "fasta"):
                rseq = str(rec.seq).upper()
                if len(rseq) >= 12:
                    add_sliding_kmers(rseq, 12)

        for fwd, rev in primer_sets:
            add_sliding_kmers(fwd, 12)
            add_sliding_kmers(rev, 12)

        baited_reads = []

        for seq in self._iter_fastq_seqs(read_files):
            for fwd, rev in primer_sets:
                rev_rc = str(Seq(rev).reverse_complement())
                if fwd in seq and rev_rc in seq:
                    s_idx = seq.find(fwd)
                    e_idx = seq.find(rev_rc, s_idx)
                    if 50 <= (e_idx - (s_idx + len(fwd))) <= 5000:
                        return seq[s_idx:e_idx + len(rev_rc)]
                seq_rc = str(Seq(seq).reverse_complement())
                if fwd in seq_rc and rev_rc in seq_rc:
                    s_idx = seq_rc.find(fwd)
                    e_idx = seq_rc.find(rev_rc, s_idx)
                    if 50 <= (e_idx - (s_idx + len(fwd))) <= 5000:
                        return seq_rc[s_idx:e_idx + len(rev_rc)]

            for i in range(0, len(seq) - 12 + 1, 6):
                if seq[i:i + 12] in all_seeds:
                    baited_reads.append(seq)
                    baited_reads.append(str(Seq(seq).reverse_complement()))
                    break

        if not baited_reads:
            return None

        for k in (135, 121, 107, 95, 77, 63, 51, 39, 27):
            for min_cov in (2, 1):
                for fwd, rev in primer_sets:
                    rev_rc = str(Seq(rev).reverse_complement())
                    amplicon = self._debruijn_walk(baited_reads, fwd, rev_rc, k=k, min_kmer_cov=min_cov)
                    if amplicon:
                        return amplicon
        return None

    @staticmethod
    def _debruijn_walk(
        reads: List[str],
        fwd_primer: str,
        rev_primer_rc: str,
        k: int = 63,
        min_kmer_cov: int = 2
    ) -> Optional[str]:
        kmer_counts: Dict[str, int] = {}
        for r in reads:
            if len(r) < k:
                continue
            for i in range(len(r) - k + 1):
                kmer = r[i:i + k]
                if "N" not in kmer:
                    kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

        valid_kmers = {km: c for km, c in kmer_counts.items() if c >= min_kmer_cov}
        if not valid_kmers:
            return None

        adj: Dict[str, List[Tuple[str, int]]] = {}
        for km, count in valid_kmers.items():
            prefix = km[:-1]
            adj.setdefault(prefix, []).append((km, count))

        for prefix in adj:
            adj[prefix].sort(key=lambda x: x[1])

        start_kmers = [km for km in valid_kmers if km.startswith(fwd_primer)]
        if not start_kmers:
            start_kmers = [km for km in valid_kmers if fwd_primer in km]
        if not start_kmers:
            return None
        start_kmers.sort(key=lambda km: valid_kmers[km], reverse=True)

        stack = [(km, km, {km: 1}) for km in reversed(start_kmers[:5])]
        max_len = 2500
        steps = 0

        while stack and steps < 25000:
            steps += 1
            curr_kmer, seq_so_far, visited = stack.pop()

            if rev_primer_rc in seq_so_far:
                s_idx = seq_so_far.find(fwd_primer)
                e_idx = seq_so_far.find(rev_primer_rc, s_idx)
                if s_idx != -1 and e_idx != -1 and (e_idx - (s_idx + len(fwd_primer))) >= 42:
                    return seq_so_far[s_idx:e_idx + len(rev_primer_rc)]

            if len(seq_so_far) > max_len:
                continue

            suffix = curr_kmer[1:]
            for next_kmer, _ in adj.get(suffix, []):
                if visited.get(next_kmer, 0) < 4:
                    new_visited = dict(visited)
                    new_visited[next_kmer] = new_visited.get(next_kmer, 0) + 1
                    stack.append((next_kmer, seq_so_far + next_kmer[-1], new_visited))

        return None