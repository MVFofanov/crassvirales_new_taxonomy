#!/usr/bin/env python3
"""
Make per-genome figures comparing genomic maps produced with different genetic codes (11, TAG, TGA).

Inputs
------
- --crassus_output: base directory; GFFs are read from 4_ORF/0_all_codings_renamed under it
- --genomes: comma-separated list OR path to a text file with one genome ID per line
- --out_dir: where to save figures (default: <base>/4_ORF/0_all_codings_renamed/genetic_code_maps)
- --dpi: image DPI (default 300)
- --svg: also save an SVG next to each PNG

Expected GFF filenames (already renamed to old_id):
  {genome_id}_tbl-11.gff
  {genome_id}_tbl-TAG.gff
  {genome_id}_tbl-TGA.gff

Each output file is named:
  {genome_id}_genetic_code_annotations.png
(We also produce a safe/sanitized filename on disk to avoid OS issues with characters like '|'.)
"""

import argparse
import pandas as pd
from pathlib import Path
import os
import re
from typing import List, Optional, Tuple, Dict

import matplotlib
matplotlib.use("Agg")  # headless
os.environ["QT_QPA_PLATFORM"] = "offscreen"
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord

# ---- Codes to show, in fixed order ----
CODES = ["11", "TAG", "TGA"]
CODE_TITLES = {
    "11": "11 (standard)",
    "TAG": "TAG reassigned → Gln",
    "TGA": "TGA reassigned → Trp",
}

# ---- A tiny color scheme (all CDS in neutral gray by default) ----
DEFAULT_COLOR = "#bfbfbf"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot genomic maps for multiple genetic codes per genome (from GFF).")
    p.add_argument(
        "--crassus_output",
        type=Path,
        required=True,
        help="Base directory of CrassUS results; GFFs will be read from 4_ORF/0_all_codings_renamed",
    )
    p.add_argument(
        "--genomes",
        type=str,
        required=True,
        help="Comma-separated list of genome IDs OR path to a file with one ID per line.",
    )
    p.add_argument(
        "--out_dir",
        type=Path,
        default=None,
        help="Where to save output figures (default: <base>/4_ORF/0_all_codings_renamed/genetic_code_maps)",
    )
    p.add_argument("--dpi", type=int, default=300, help="DPI for PNG outputs")
    p.add_argument("--svg", action="store_true", help="Also save SVG next to each PNG")
    return p.parse_args()


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def parse_genomes_arg(genomes_arg: str) -> List[str]:
    path = Path(genomes_arg)
    if path.exists() and path.is_file():
        return [line.strip() for line in path.read_text().splitlines() if line.strip()]
    # else treat as comma-separated list
    return [x.strip() for x in genomes_arg.split(",") if x.strip()]


def safe_filename(name: str) -> str:
    """Replace characters that are problematic on some filesystems (e.g. '|' on NTFS)."""
    return re.sub(r"[^A-Za-z0-9._-]", "__", name)


def gff_path_for(genome_id: str, code: str, gff_dir: Path) -> Path:
    return gff_dir / f"{genome_id}_tbl-{code}.gff"


def read_gff_features(gff_path: Path) -> List[GraphicFeature]:
    """Parse a minimal GFF (Prodigal-style). Keep only CDS, create GraphicFeature for each."""
    feats: List[GraphicFeature] = []
    if not gff_path.exists():
        return feats
    for line in gff_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = parts
        if ftype != "CDS":
            continue
        try:
            s = int(start)
            e = int(end)
        except ValueError:
            continue
        st = 1 if strand == "+" else -1
        feats.append(GraphicFeature(start=s, end=e, strand=st, color=DEFAULT_COLOR, label=None))
    return feats

def _extract_code_from_filename(name: str) -> Optional[str]:
    m = re.search(r"_tbl-(11|TAG|TGA)\.gff$", name)
    return m.group(1) if m else None


def detect_best_code(genome_id: str, base: Path) -> Optional[str]:
    """
    Look only in 4_ORF/1_best_coding_renamed for {genome_id}_tbl-<CODE>.gff
    and return <CODE> if found; else None.
    """
    d = base / "4_ORF" / "1_best_coding_renamed"
    if not d.exists():
        return None
    candidates = list(d.glob(f"{genome_id}_tbl-*.gff"))
    if not candidates:
        return None
    return _extract_code_from_filename(candidates[0].name)

def load_functional_annotations(base: Path) -> pd.DataFrame:
    """
    Load 4_ORF/3_functional_annot_table_renamed.tsv.
    Returns an empty DataFrame if the file is missing.
    """
    tsv = base / "4_ORF" / "3_functional_annot_table_renamed.tsv"
    if not tsv.exists():
        return pd.DataFrame(columns=["protein_id","genome","start","end","strand","partial","yutin","cluster"])
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    # normalize types
    if "start" in df.columns:
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    if "end" in df.columns:
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    if "strand" in df.columns:
        df["strand"] = df["strand"].astype(str).str.strip()
    if "yutin" in df.columns:
        # Keep proper missing values (<NA>) instead of the literal "nan" string
        df["yutin"] = df["yutin"].astype("string").str.strip()
        # Treat case-insensitive "nan" and empty strings as missing
        df.loc[df["yutin"].str.lower().eq("nan"), "yutin"] = pd.NA
        df["yutin"] = df["yutin"].replace("", pd.NA)

    return df


def yutin_index_for_genome(func_df: pd.DataFrame, genome_id: str) -> Dict[tuple, str]:
    """
    Build a lookup: (start, end, strand_sign) -> yutin_label
    Only includes rows for this genome with non-empty yutin.
    strand_sign is  1 for '+', -1 for anything else.
    """
    if func_df is None or func_df.empty:
        return {}
    sub = func_df[(func_df["genome"] == genome_id)]
    if "yutin" not in sub.columns:
        return {}
    # After normalization above, just drop missing
    sub = sub[sub["yutin"].notna()]
    # ensure ints
    sub = sub.dropna(subset=["start", "end"])
    idx: Dict[tuple, str] = {}
    for _, r in sub.iterrows():
        try:
            s = int(r["start"]); e = int(r["end"])
        except Exception:
            continue
        strand_sign = 1 if str(r.get("strand", "+")) == "+" else -1
        idx[(s, e, strand_sign)] = str(r["yutin"])
    return idx


def annotate_best_track_features(features: List[GraphicFeature], yidx: Dict[tuple, str]) -> None:
    """
    In-place: for features whose (start,end,strand) match yidx, add label and color red.
    Leave all other features untouched.
    """
    if not features or not yidx:
        return
    for f in features:
        key = (int(f.start), int(f.end), int(f.strand))
        label = yidx.get(key)
        if label is not None and str(label).strip():
            f.label = str(label)
            f.color = "red"


def plot_genome(genome_id: str, gff_dir: Path, out_dir: Path, base: Path,
                func_df: Optional[pd.DataFrame] = None,
                dpi: int = 300, save_svg: bool = False) -> Tuple[Path, Optional[Path]]:
    # Detect best code from 4_ORF/1_best_coding_renamed
    best_code = detect_best_code(genome_id, base)

    # Load features for each code and compute a common sequence length
    code_to_feats: Dict[str, List[GraphicFeature]] = {}
    max_end = 0
    for code in CODES:
        path = gff_path_for(genome_id, code, gff_dir)
        feats = read_gff_features(path)
        code_to_feats[code] = feats
        if feats:
            max_end = max(max_end, max(f.end for f in feats))

    # If we have a best code and a functional table, annotate only that track
    if best_code in code_to_feats and func_df is not None and not func_df.empty:
        yidx = yutin_index_for_genome(func_df, genome_id)
        annotate_best_track_features(code_to_feats[best_code], yidx)

    seq_len = max(1000, max_end + 100)  # minimal baseline if empty

    # Build records per code
    records = {code: GraphicRecord(sequence_length=seq_len, features=feats)
               for code, feats in code_to_feats.items()}

    # Figure layout: one row per code
    height_per_row = 2.0
    fig_h = max(2.5, height_per_row * len(CODES))
    fig, axes = plt.subplots(len(CODES), 1, figsize=(14, fig_h))
    if len(CODES) == 1:
        axes = [axes]

    # Suptitle with genome id (include best code if known)
    sup = genome_id if best_code is None else f"{genome_id}   (CrassUS best code: {CODE_TITLES.get(best_code, best_code)})"
    fig.suptitle(sup, fontsize=12, y=0.995)

    for ax, code in zip(axes, CODES):
        records[code].plot(ax=ax)
        title = f"Genetic code: {CODE_TITLES.get(code, code)}"
        if best_code is not None and code == best_code:
            title += " — CrassUS BEST"
        ax.set_title(title, loc="left", fontsize=10)

    plt.tight_layout(rect=[0, 0, 1, 0.97])

    # Save
    ensure_dir(out_dir)
    out_png = out_dir / f"{safe_filename(genome_id)}_genetic_code_annotations.png"
    fig.savefig(out_png, dpi=dpi)

    out_svg = None
    if save_svg:
        out_svg = out_dir / f"{safe_filename(genome_id)}_genetic_code_annotations.svg"
        fig.savefig(out_svg)

    plt.close(fig)
    return out_png, out_svg



def main() -> None:
    args = parse_args()
    base = args.crassus_output
    gff_dir = base / "4_ORF" / "0_all_codings_renamed"
    out_dir = args.out_dir or (gff_dir / "genetic_code_maps")

    genomes = parse_genomes_arg(args.genomes)

    print(f"[INFO] Reading GFFs from: {gff_dir}")
    print(f"[INFO] Saving figures to: {out_dir}")

    func_df = load_functional_annotations(base)

    for g in genomes:
        png, svg = plot_genome(g, gff_dir, out_dir, base=base, func_df=func_df, dpi=args.dpi, save_svg=args.svg)
        if svg is None:
            print(f"✅ Saved: {png}")
        else:
            print(f"✅ Saved: {png} and {svg}")


if __name__ == "__main__":
    main()
