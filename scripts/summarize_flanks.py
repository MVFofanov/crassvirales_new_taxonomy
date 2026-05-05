#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import pandas as pd


# ---------- Helpers to read IDs ----------

def read_prophage_ids(ann_path: Path, limit: int) -> List[str]:
    """Read first `limit` unique genome_id values, preserving file order."""
    df = pd.read_csv(ann_path, sep="\t", dtype=str)
    if "genome_id" not in df.columns:
        raise ValueError(f"'genome_id' column not found in {ann_path}")
    seen: Set[str] = set()
    ids: List[str] = []
    for val in df["genome_id"].dropna().astype(str):
        if val not in seen:
            seen.add(val)
            ids.append(val)
            if len(ids) >= limit:
                break
    return ids


# ---------- Finding BLAST files (handles | vs & in folder names) ----------

def _blast_dir_candidates(wd: Path, prophage_id: str) -> List[Path]:
    raw = wd / prophage_id / "blast"
    amp = wd / prophage_id.replace("|", "&") / "blast"
    enc = wd / prophage_id.replace("|", "%7C") / "blast"
    return [raw, amp, enc]


def find_flank_files(wd: Path, prophage_id: str, debug: bool = False) -> List[Path]:
    for cand in _blast_dir_candidates(wd, prophage_id):
        if debug:
            print(f"[debug] Full path: {cand}")
        if cand.is_dir():
            return sorted(p for p in cand.glob("*-Alignment.txt") if p.is_file())
    print(f"[warn] No blast dir for {prophage_id}")
    return []


# ---------- Parsing & union logic ----------

def _to_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except Exception:
        return None


def _to_int(s: str) -> Optional[int]:
    try:
        return int(float(s))
    except Exception:
        return None


def _union_length(intervals: List[Tuple[int, int]], merge_gap: int = 0) -> int:
    """
    Length of union of closed intervals [s,e] with optional gap-merging.
    merge_gap=0 merges overlaps only; merge_gap=1 also merges touching intervals.
    """
    if not intervals:
        return 0
    ivals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    cs, ce = ivals[0]
    for s, e in ivals[1:]:
        if s <= ce + merge_gap:  # overlap or touch
            if e > ce:
                ce = e
        else:
            total += (ce - cs + 1)
            cs, ce = s, e
    total += (ce - cs + 1)
    return total


_side_left = re.compile(r"\bleft\b", re.IGNORECASE)
_side_right = re.compile(r"\bright\b", re.IGNORECASE)

def _flank_side(qid: str) -> Optional[str]:
    """Return 'left', 'right', or None based on query ID text."""
    q = qid.lower()
    # strong signals with underscores (your naming scheme)
    if "_left_" in q or q.endswith("_left") or q.startswith("left_"):
        return "left"
    if "_right_" in q or q.endswith("_right") or q.startswith("right_"):
        return "right"
    # fallbacks for other separators just in case
    if re.search(r'(?<![a-z0-9])left(?![a-z0-9])', q):
        return "left"
    if re.search(r'(?<![a-z0-9])right(?![a-z0-9])', q):
        return "right"
    return None

def _strip_ver(acc: str) -> str:
    return acc.split(".", 1)[0] if "." in acc else acc


def accumulate_union_query_stats(
    files: Sequence[Path],
    query_col: int = 1,      # 1-based: query acc.ver
    subject_col: int = 2,    # 1-based: subject acc.ver
    ident_col: int = 3,      # 1-based: % identity
    aln_len_col: int = 4,    # 1-based: alignment length
    qstart_col: int = 7,     # 1-based: q.start
    qend_col: int = 8,       # 1-based: q.end
    min_identity: Optional[float] = None,
    min_aln_len: Optional[int] = None,
    merge_touching: bool = False,
) -> Dict[str, Tuple[int, int, int, int]]:
    """
    Build subject -> (total_union_len_on_query, n_flanks, n_left, n_right).

    Steps:
      - collect HSP intervals on the query for each (subject, query) after filters,
      - union per (subject, query),
      - sum unions over queries per subject,
      - count how many distinct queries had ≥1 passing HSP overall, and split by left/right.
    """
    q_idx = query_col - 1
    s_idx = subject_col - 1
    i_idx = ident_col - 1
    l_idx = aln_len_col - 1
    qs_idx = qstart_col - 1
    qe_idx = qend_col - 1

    # subject -> query -> list of (qstart, qend)
    subj_query_intervals: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    # sets for counting distinct flanks hit per subject
    subj_seen: Dict[str, Set[str]] = {}
    subj_seen_left: Dict[str, Set[str]] = {}
    subj_seen_right: Dict[str, Set[str]] = {}

    for path in files:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if max(q_idx, s_idx, l_idx, qs_idx, qe_idx) >= len(parts):
                    continue

                subj = parts[s_idx].strip()
                qid = parts[q_idx].strip()
                if not subj or not qid:
                    continue

                # Filters
                aln_len = _to_int(parts[l_idx])
                if aln_len is None:
                    continue
                if min_aln_len is not None and aln_len < min_aln_len:
                    continue

                if min_identity is not None and i_idx < len(parts):
                    ident = _to_float(parts[i_idx])
                    if ident is None or ident < min_identity:
                        continue

                qstart = _to_int(parts[qs_idx])
                qend = _to_int(parts[qe_idx])
                if qstart is None or qend is None:
                    continue

                s, e = (qstart, qend) if qstart <= qend else (qend, qstart)

                subj_query_intervals.setdefault(subj, {}).setdefault(qid, []).append((s, e))

                # record counts per flank (distinct queries)
                subj_seen.setdefault(subj, set()).add(qid)
                side = _flank_side(qid)
                if side == "left":
                    subj_seen_left.setdefault(subj, set()).add(qid)
                elif side == "right":
                    subj_seen_right.setdefault(subj, set()).add(qid)

    # Compute union per query, then sum per subject; add counts
    gap = 1 if merge_touching else 0
    stats: Dict[str, Tuple[int, int, int, int]] = {}
    for subj, qmap in subj_query_intervals.items():
        tot_len = 0
        for qid, ivals in qmap.items():
            tot_len += _union_length(ivals, merge_gap=gap)

        n_all = len(subj_seen.get(subj, set()))
        n_left = len(subj_seen_left.get(subj, set()))
        n_right = len(subj_seen_right.get(subj, set()))
        stats[subj] = (tot_len, n_all, n_left, n_right)

    return stats


# ---------- End-to-end summarization ----------

def summarize_by_union_len(
    wd: Path,
    prophage_ids: Sequence[str],
    top_n: int = 50,
    query_col: int = 1,
    subject_col: int = 2,
    ident_col: int = 3,
    aln_len_col: int = 4,
    qstart_col: int = 7,
    qend_col: int = 8,
    min_identity: Optional[float] = None,
    min_aln_len: Optional[int] = None,
    merge_touching: bool = False,
    debug: bool = False,
    # NEW:
    min_total_len: int = 10000,
) -> pd.DataFrame:
    """
    For each prophage, pool all flank BLASTs, rank subjects by TOTAL union length
    on query coords, keep Top-N, and return a TSV-ready DataFrame with columns:

      [prophage_id, bacterial_contig_id, total_alignment_length,
       n_flanks, n_flanks_left, n_flanks_right]

    Only subjects with total_alignment_length >= min_total_len are kept.
    """
    rows: List[Tuple[str, str, int, int, int, int]] = []

    for pid in prophage_ids:
        files = find_flank_files(wd, pid, debug=debug)
        subj_stats = accumulate_union_query_stats(
            files,
            query_col=query_col,
            subject_col=subject_col,
            ident_col=ident_col,
            aln_len_col=aln_len_col,
            qstart_col=qstart_col,
            qend_col=qend_col,
            min_identity=min_identity,
            min_aln_len=min_aln_len,
            merge_touching=merge_touching,
        )

        # NEW: filter by minimum total unioned length (bp)
        if min_total_len and min_total_len > 0:
            subj_stats = {s: v for s, v in subj_stats.items() if v[0] >= min_total_len}

        # rank by total union length desc
        ranked = sorted(subj_stats.items(), key=lambda kv: kv[1][0], reverse=True)[:top_n]
        print(f"{pid}\tunique_matches={len(ranked)}")
        for subj, (tot, n_all, n_left, n_right) in ranked:
            rows.append((pid, subj, tot, n_all, n_left, n_right))

    df = pd.DataFrame(
        rows,
        columns=[
            "prophage_id",
            "bacterial_contig_id",
            "total_alignment_length",
            "n_flanks",
            "n_flanks_left",
            "n_flanks_right",
        ],
    )
    if not df.empty:
        df = df.sort_values(
            by=["prophage_id", "total_alignment_length", "n_flanks"],
            ascending=[True, False, False],
            kind="mergesort",
        )
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Pool all flank BLASTs per prophage; for each subject, sum the UNION "
            "of HSP lengths on query coordinates (per query), then keep Top-N per "
            "prophage by that total. Adds the number of flanks hit (total/left/right)."
        )
    )
    ap.add_argument("--wd", required=True, type=Path, help="Working dir with {prophage_id}/blast/*-Alignment.txt")
    ap.add_argument("--ann_file", required=True, type=Path, help="Annotation table (TSV by default; CSV if .csv)")
    ap.add_argument("--top_matches", type=int, default=50, help="Top-N subjects per prophage [50]")
    ap.add_argument("--limit", type=int, default=16, help="Number of first genome_id to use [16]")

    # Column indices (1-based), matching your BLAST header
    ap.add_argument("--query_col", type=int, default=1, help="1-based column index for query acc.ver [1]")
    ap.add_argument("--subject_col", type=int, default=2, help="1-based column index for subject acc.ver [2]")
    ap.add_argument("--ident_col", type=int, default=3, help="1-based column index for % identity [3]")
    ap.add_argument("--aln_len_col", type=int, default=4, help="1-based column index for alignment length [4]")
    ap.add_argument("--qstart_col", type=int, default=7, help="1-based column index for q.start [7]")
    ap.add_argument("--qend_col", type=int, default=8, help="1-based column index for q.end [8]")

    # Filters
    ap.add_argument("--min_identity", type=float, default=None, help="Min % identity per HSP (optional)")
    ap.add_argument("--min_aln_len", type=int, default=None, help="Min alignment length per HSP (optional)")
    ap.add_argument(
        "--merge_touching",
        action="store_true",
        help="Merge HSPs that just touch on query (treat end and next start as continuous)",
    )

    ap.add_argument("--out_tsv", type=Path, default=Path("flank_blast_topN_union_query.tsv"), help="Output TSV")
    ap.add_argument("--debug", action="store_true", help="Verbose path checking")

    ap.add_argument(
        "--min_total_len",
        type=int,
        default=10000,
        help="Minimum total unioned alignment length (bp) per subject to keep [10000]. Use 0 to disable."
    )

    # --- add near other argparse options ---
    ap.add_argument(
        "--out_id_list",
        type=Path,
        required=True,
        help="Path to write a unique list of contig IDs (one per line).",
    )
    ap.add_argument(
        "--strip_version_in_list",
        action="store_true",
        help="Strip version suffix (e.g., '.1') from accessions in the ID list.",
    )

    args = ap.parse_args()

    prophage_ids = read_prophage_ids(args.ann_file, limit=args.limit)
    if not prophage_ids:
        raise SystemExit("No prophage_ids found in annotation table.")

    df = summarize_by_union_len(
        wd=args.wd,
        prophage_ids=prophage_ids,
        top_n=args.top_matches,
        query_col=args.query_col,
        subject_col=args.subject_col,
        ident_col=args.ident_col,
        aln_len_col=args.aln_len_col,
        qstart_col=args.qstart_col,
        qend_col=args.qend_col,
        min_identity=args.min_identity,
        min_aln_len=args.min_aln_len,
        merge_touching=args.merge_touching,
        debug=args.debug,
        min_total_len=args.min_total_len,
    )

    total_unique = df["bacterial_contig_id"].nunique() if not df.empty else 0
    print(f"TOTAL_UNIQUE_bacterial_contig_id\t{total_unique}")

    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"Wrote TSV: {args.out_tsv.resolve()}")

    # host contigs from the selected prophage_ids
    host_contigs = {pid.split("|", 1)[0] for pid in prophage_ids}

    # subject contigs from the output table (may be empty)
    subject_contigs = set(df["bacterial_contig_id"].astype(str)) if not df.empty else set()

    all_ids = host_contigs | subject_contigs

    # optional: strip version
    if args.strip_version_in_list:
        all_ids = {_strip_ver(x) for x in all_ids}

    # write one ID per line
    args.out_id_list.parent.mkdir(parents=True, exist_ok=True)
    with args.out_id_list.open("w") as fh:
        for acc in sorted(all_ids):
            fh.write(f"{acc}\n")

    print(f"Wrote {len(all_ids)} unique contig IDs to: {args.out_id_list.resolve()}")


if __name__ == "__main__":
    main()
