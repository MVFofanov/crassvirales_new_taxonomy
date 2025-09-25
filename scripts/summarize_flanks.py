#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd


def read_prophage_ids(ann_path: Path, limit: int) -> List[str]:
    """
    Read the first `limit` unique genome_id values from the annotation table,
    preserving file order. Supports TSV or CSV by extension.
    """
    sep = "," if ann_path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(ann_path, sep="\t", dtype=str)
    if "genome_id" not in df.columns:
        raise ValueError(f"'genome_id' column not found in {ann_path}")
    ids: List[str] = []
    seen: Set[str] = set()
    for val in df["genome_id"].dropna().astype(str):
        if val not in seen:
            seen.add(val)
            ids.append(val)
            if len(ids) >= limit:
                break
    return ids


# --- add this helper near the top ---
def _resolve_prophage_dir(wd: Path, pid: str) -> Optional[Path]:
    """
    Return an existing directory for this prophage id, trying common
    filesystem-safe variants (e.g., replacing '|' with '&' or '_').
    """
    # candidates = [
    #     pid,
    #     pid.replace("|", "&")
    #     # pid.replace("|", "_"),
    #     # pid.replace("|", "%7C"),
    # ]

    cand = pid.replace("|", "&")
    # print(f"[debug] wd: {wd}")
    # print(f"[debug] Trying candidate dir: {cand}")
    d = wd / cand / "blast"
    print(f"[debug] Full path: {d}")
    if d.is_dir():
        return d
    return None


# --- replace your current find_flank_files with this ---
def find_flank_files(wd: Path, prophage_id: str, file_glob: str = "*-Alignment.txt") -> List[Path]:
    """
    Return all flank BLAST files for a given prophage, e.g.:
    {wd}/{prophage_id or sanitized}/blast/*-Alignment.txt
    """
    blast_dir = _resolve_prophage_dir(wd, prophage_id)
    if blast_dir is None:
        # optional: print a one-line hint for debugging
        print(f"[warn] No blast dir for {prophage_id}")
        return []
    files = sorted(p for p in blast_dir.glob(file_glob) if p.is_file())
    if not files:
        print(f"[warn] No files matched {file_glob} under {blast_dir}")
    return files


def _parse_float_safe(s: str) -> Optional[float]:
    try:
        return float(s)
    except Exception:
        return None


def top_unique_subjects_from_blast(
    path: Path,
    subject_col: int = 2,
    bitscore_col: Optional[int] = 12,
    evalue_col: Optional[int] = 11,
    top_n: int = 50,
) -> List[str]:
    """
    Parse a BLAST tabular file (comment lines starting with '#'
    are ignored, fields are tab-separated). Picks the top-N UNIQUE
    subject IDs by best bitscore (descending), ties by lowest evalue.
    Column numbers are 1-based (like `cut -f2`), per BLAST outfmt 6/7.
    Falls back to first-seen order if score columns are not available.
    """
    s_idx = subject_col - 1
    b_idx = bitscore_col - 1 if bitscore_col else None
    e_idx = evalue_col - 1 if evalue_col else None

    # subject_id -> (best_bitscore or -inf, best_evalue or +inf)
    best: Dict[str, Tuple[float, float]] = {}
    insertion_order: Dict[str, int] = {}

    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if s_idx >= len(parts):
                continue
            subj = parts[s_idx].strip()
            if not subj:
                continue

            bs = None
            ev = None
            if b_idx is not None and b_idx < len(parts):
                bs = _parse_float_safe(parts[b_idx])
            if e_idx is not None and e_idx < len(parts):
                ev = _parse_float_safe(parts[e_idx])

            if subj not in insertion_order:
                insertion_order[subj] = i

            # initialize with worst possible if missing
            cur_bs, cur_ev = best.get(subj, (float("-inf"), float("inf")))
            cand_bs = bs if bs is not None else cur_bs
            cand_ev = ev if ev is not None else cur_ev

            # better if higher bitscore; if equal, lower evalue; if both missing, keep first
            improved = False
            if bs is not None and cand_bs > cur_bs:
                improved = True
            elif bs is not None and cand_bs == cur_bs and ev is not None and cand_ev < cur_ev:
                improved = True

            if subj not in best:
                best[subj] = (cand_bs, cand_ev)
            elif improved:
                best[subj] = (cand_bs, cand_ev)

    # Sort by best score if available; otherwise by first-seen insertion order
    def sort_key(subj: str) -> Tuple[int, float, float, int]:
        bs, ev = best.get(subj, (float("-inf"), float("inf")))
        # primary: whether we have a real score (1) vs not (0) -> scored first
        has_score = int(bs != float("-inf") or ev != float("inf"))
        # then: bitscore desc (-bs), evalue asc (ev), then insertion order asc
        return (has_score, -bs, ev, insertion_order.get(subj, 1_000_000))

    ranked = sorted(best.keys(), key=sort_key, reverse=False)
    return ranked[:top_n]


def summarize_prophage_matches(
    wd: Path,
    prophage_ids: Sequence[str],
    top_n: int = 50,
    subject_col: int = 2,
    bitscore_col: Optional[int] = 12,
    evalue_col: Optional[int] = 11,
) -> pd.DataFrame:
    """
    For each prophage_id, aggregate the union of top-N unique subjects
    across all flank BLAST files. Return a tidy DataFrame with columns:
    prophage_id, bacterial_contig_id (unique per prophage).
    """
    rows: List[Tuple[str, str]] = []
    for pid in prophage_ids:
        flank_files = find_flank_files(wd, pid)
        per_prophage: Set[str] = set()
        for f in flank_files:
            subjects = top_unique_subjects_from_blast(
                f,
                subject_col=subject_col,
                bitscore_col=bitscore_col,
                evalue_col=evalue_col,
                top_n=top_n,
            )
            per_prophage.update(subjects)
        # add unique rows per prophage
        for subj in sorted(per_prophage):
            rows.append((pid, subj))
        # print per-prophage unique count
        print(f"{pid}\tunique_matches={len(per_prophage)}")
    df = pd.DataFrame(rows, columns=["prophage_id", "bacterial_contig_id"])
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Summarize BLAST flank matches for prophages: for each flanking "
            "region file, take top-N unique subject matches, aggregate per "
            "prophage, and output a tidy table."
        )
    )
    ap.add_argument("--wd", required=True, type=Path, help="Working dir containing {prophage_id}/blast/*-Alignment.txt")
    ap.add_argument("--ann_file", required=True, type=Path, help="Prophage annotation table (TSV by default; CSV if .csv)")
    ap.add_argument("--top_matches", type=int, default=50, help="Top-N unique subjects per flank [default: 50]")
    ap.add_argument("--limit", type=int, default=16, help="Number of first genome_id to use [default: 16]")
    ap.add_argument("--subject_col", type=int, default=2, help="1-based column index for subject ID [default: 2]")
    ap.add_argument("--bitscore_col", type=int, default=12, help="1-based column index for bitscore (outfmt 6) [default: 12]")
    ap.add_argument("--evalue_col", type=int, default=11, help="1-based column index for evalue (outfmt 6) [default: 11]")
    ap.add_argument("--out_tsv", type=Path, default=Path("flank_blast_summary.csv"), help="Output CSV path")

    args = ap.parse_args()

    prophage_ids = read_prophage_ids(args.ann_file, limit=args.limit)
    if not prophage_ids:
        raise SystemExit("No prophage_ids found in annotation table.")

    df = summarize_prophage_matches(
        wd=args.wd,
        prophage_ids=prophage_ids,
        top_n=args.top_matches,
        subject_col=args.subject_col,
        bitscore_col=args.bitscore_col,
        evalue_col=args.evalue_col,
    )

    # Print grand total unique bacterial_contig_id across whole table
    total_unique = df["bacterial_contig_id"].nunique()
    print(f"TOTAL_UNIQUE_bacterial_contig_id\t{total_unique}")

    # Save output as TSV
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"Wrote TSV: {args.out_tsv.resolve()}")


if __name__ == "__main__":
    main()
