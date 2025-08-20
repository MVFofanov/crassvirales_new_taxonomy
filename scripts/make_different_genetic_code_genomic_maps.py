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


def plot_genome(genome_id: str, gff_dir: Path, out_dir: Path, dpi: int = 300, save_svg: bool = False) -> Tuple[Path, Optional[Path]]:
    # Load features for each code and compute a common sequence length
    code_to_feats: Dict[str, List[GraphicFeature]] = {}
    max_end = 0
    for code in CODES:
        path = gff_path_for(genome_id, code, gff_dir)
        feats = read_gff_features(path)
        code_to_feats[code] = feats
        if feats:
            max_end = max(max_end, max(f.end for f in feats))

    seq_len = max(1000, max_end + 100)  # minimal baseline if empty

    # Build records per code
    records = {code: GraphicRecord(sequence_length=seq_len, features=feats) for code, feats in code_to_feats.items()}

    # Figure layout: one row per code
    height_per_row = 2.0
    fig_h = max(2.5, height_per_row * len(CODES))
    fig, axes = plt.subplots(len(CODES), 1, figsize=(14, fig_h))
    if len(CODES) == 1:
        axes = [axes]

    # Suptitle with genome id
    fig.suptitle(genome_id, fontsize=12, y=0.995)

    for ax, code in zip(axes, CODES):
        records[code].plot(ax=ax)
        ax.set_title(f"Genetic code: {CODE_TITLES.get(code, code)}", loc="left", fontsize=10)

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
    for g in genomes:
        png, svg = plot_genome(g, gff_dir, out_dir, dpi=args.dpi, save_svg=args.svg)
        if svg is None:
            print(f"✅ Saved: {png}")
        else:
            print(f"✅ Saved: {png} and {svg}")


if __name__ == "__main__":
    main()
