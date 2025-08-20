#!/usr/bin/env python3
"""
Plot per-base read coverage over bacterial contigs with prophage overlays
AND report mean coverage for the whole contig and for the prophage.

Directory layout example (wd):
  NZ_FOJF01000001.1&provirus_3692071_3771634/
    basecov.SRR4235163.txt
  NZ_JAAFZH010000004.1&provirus_416622_568155/
    basecov.SRR11574547.txt
  NZ_PVTE01000009.1&provirus_114479_228771/
    basecov.SRR6479507.txt
  Crassphage_prophage_analysis_annotation.tsv

Outputs: figures/<dir>_basecov.png for each prophage.
Prints:  contig and prophage mean coverage to stdout.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

# ---------- configurable defaults ----------
DEFAULT_DOWNSAMPLE_TARGET = 5000  # ~number of points to plot; 0 disables
FIG_WIDTH = 12
FIG_HEIGHT = 6
LINE_WIDTH = 0.9
ALPHA_SHADE = 0.18
COVERAGE_PCTL_FOR_YMAX = 99.5  # clip ylim to this percentile
# ------------------------------------------


@dataclass
class ProphageInfo:
    dir_name: str
    contig_accession: str          # exactly as in annotation
    start_1based: int
    end_1based: int
    contig_length: int
    prophage_length: int           # <-- NEW
    srr_id: Optional[str]
    family: Optional[str]
    origin: Optional[str]
    tax_lineage: Optional[str]

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot coverage with prophage overlay.")
    p.add_argument(
        "--wd",
        type=Path,
        default=Path("/mnt/c/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/integrated_prophages"),
        help="Working directory containing prophage subfolders.",
    )
    p.add_argument(
        "--annotation",
        type=Path,
        default=Path("Crassphage_prophage_analysis_annotation.tsv"),
        help="Annotation TSV with contig length and prophage coords.",
    )
    p.add_argument(
        "--figdir",
        type=Path,
        default=Path("figures"),
        help="Output directory for figures.",
    )
    p.add_argument(
        "--downsample_target",
        type=int,
        default=DEFAULT_DOWNSAMPLE_TARGET,
        help="Target points to plot (~bins). Use 0 to disable.",
    )
    return p.parse_args()


def load_annotation_table(path: Path) -> pd.DataFrame:
    """Load and type-cast the annotation table (tab-delimited)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = ["genome_id", "bacterial_id", "bacterial_contig_length", "prophage_start", "prophage_end"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Annotation table missing required column: {col}")
    df["bacterial_contig_length"] = df["bacterial_contig_length"].astype(int)
    df["prophage_start"] = df["prophage_start"].astype(int)
    df["prophage_end"] = df["prophage_end"].astype(int)
    return df


def iter_prophage_dirs(wd: Path) -> Iterable[Path]:
    """Yield subdirectories named like '<accession>&provirus_<start>_<end>'."""
    for sub in sorted(wd.iterdir()):
        if sub.is_dir() and "&provirus_" in sub.name:
            yield sub


def parse_dir_name(dir_name: str) -> Tuple[str, int, int]:
    """
    'NZ_FOJF01000001.1&provirus_3692071_3771634'
      -> ('NZ_FOJF01000001.1', 3692071, 3771634)
    """
    left, right = dir_name.split("&provirus_", 1)
    s_str, e_str, *_ = right.split("_")
    return left, int(s_str), int(e_str)


def find_annotation_row_for_dir(
    annot: pd.DataFrame, contig_acc: str, start_1b: int, end_1b: int
) -> Optional[pd.Series]:
    """
    Match by bacterial_id and exact coords.
    Fallback: genome_id equals '<contig_acc>|provirus_<start>_<end>'.
    """
    candidates = annot[
        (annot["bacterial_id"] == contig_acc)
        & (annot["prophage_start"] == start_1b)
        & (annot["prophage_end"] == end_1b)
    ]
    if len(candidates) == 0:
        genome_id_guess = f"{contig_acc}|provirus_{start_1b}_{end_1b}"
        candidates = annot[annot["genome_id"] == genome_id_guess]
    if len(candidates) == 0:
        return None
    return candidates.iloc[0]


def find_basecov_file(prophage_dir: Path) -> Optional[Path]:
    """Return the first basecov*.txt file in the dir, or None."""
    matches = sorted(prophage_dir.glob("basecov*.txt"))
    return matches[0] if matches else None


def parse_srr_from_filename(p: Path) -> Optional[str]:
    """Extract SRR/ERR/DRR from e.g. 'basecov.SRR4235163.txt'."""
    stem = p.stem  # 'basecov.SRR4235163'
    for token in stem.split("."):
        if token.startswith(("SRR", "ERR", "DRR")):
            return token
    return None


def load_basecov_filtered(basecov_path: Path, contig_acc: str) -> pd.DataFrame:
    """
    Load BBMap basecov (tab-delimited) and filter to the contig rows.
    Columns: 'RefName' (may contain spaces), 'Pos' (0-based), 'Coverage'.
    Compare accession to the first token of RefName before the first space.
    Returns DataFrame with columns: pos1 (1-based), cov (float).
    """
    df = pd.read_csv(basecov_path, sep="\t", engine="python", dtype=str)
    df.columns = [c.strip().lstrip("#") for c in df.columns]  # normalize header (#RefName -> RefName)
    expected = {"RefName", "Pos", "Coverage"}
    if not expected.issubset(df.columns):
        raise ValueError(f"Unexpected basecov columns in {basecov_path.name}: {df.columns.tolist()}")

    acc_token = df["RefName"].astype(str).str.split(" ", n=1, expand=True)[0]
    sub = df.loc[acc_token == contig_acc, ["Pos", "Coverage"]].copy()
    if sub.empty:
        raise ValueError(f"No rows in {basecov_path.name} match accession '{contig_acc}'")

    sub["pos1"] = sub["Pos"].astype(int) + 1  # 0-based -> 1-based
    sub["cov"] = sub["Coverage"].astype(float)
    sub = sub[["pos1", "cov"]].sort_values("pos1")
    return sub


def downsample_coverage(df: pd.DataFrame, contig_len: int, target_points: int) -> pd.DataFrame:
    """Adaptive binning to ~target_points for fast/clean plotting."""
    if target_points <= 0 or contig_len <= target_points:
        return df
    bin_size = max(1, int(np.ceil(contig_len / target_points)))
    bins = ((df["pos1"] - 1) // bin_size).astype(int)
    return df.groupby(bins, as_index=False).agg(pos1=("pos1", "median"), cov=("cov", "mean"))


def thousands(x: float, pos: int) -> str:
    return f"{int(x):,}".replace(",", " ")


def format_bp(n: int) -> str:
    """Human-friendly bp string, with a kb/Mb hint."""
    if n >= 1_000_000:
        return f"{n:,} bp ({n/1_000_000:.2f} Mb)".replace(",", " ")
    if n >= 1_000:
        return f"{n:,} bp ({n/1_000:.1f} kb)".replace(",", " ")
    return f"{n:,} bp".replace(",", " ")


def compute_means(cov_df: pd.DataFrame, start_1b: int, end_1b: int) -> Tuple[float, float]:
    """
    Compute mean coverage across the entire contig (cov_df)
    and within the prophage [start_1b, end_1b].
    """
    contig_mean = float(cov_df["cov"].mean())
    sub = cov_df[(cov_df["pos1"] >= start_1b) & (cov_df["pos1"] <= end_1b)]
    prophage_mean = float(sub["cov"].mean()) if not sub.empty else float("nan")
    return contig_mean, prophage_mean


def plot_prophage_coverage(
    info: ProphageInfo,
    cov_df: pd.DataFrame,
    out_png: Path,
    downsample_target: int = DEFAULT_DOWNSAMPLE_TARGET,
) -> Tuple[float, float]:
    """
    Make 2-row figure: coverage (top) + contig/prophage bars (bottom).
    Returns (contig_mean, prophage_mean).
    """
    # Means (computed on full-resolution data)
    contig_mean = float(cov_df["cov"].mean())
    sub = cov_df[(cov_df["pos1"] >= info.start_1based) & (cov_df["pos1"] <= info.end_1based)]
    prophage_mean = float(sub["cov"].mean()) if not sub.empty else float("nan")

    # Downsampled view for plotting
    cov_plot = downsample_coverage(cov_df, info.contig_length, downsample_target)

    fig, (ax_cov, ax_bar) = plt.subplots(
        2, 1, figsize=(FIG_WIDTH, FIG_HEIGHT), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )

    # Coverage line + shaded prophage
    ax_cov.plot(cov_plot["pos1"].values, cov_plot["cov"].values, lw=LINE_WIDTH)
    ax_cov.axvspan(info.start_1based, info.end_1based, color="red", alpha=ALPHA_SHADE, lw=0)

    # Add dashed mean lines
    ax_cov.axhline(contig_mean, linestyle="--", linewidth=0.9, color="0.35", label=f"Contig mean: {contig_mean:.1f}×")
    if np.isfinite(prophage_mean):
        ax_cov.axhline(prophage_mean, linestyle="--", linewidth=0.9, color="red",
                       label=f"Prophage mean: {prophage_mean:.1f}×")

    # y-limit clip for readability
    ymax = np.percentile(cov_plot["cov"].values, COVERAGE_PCTL_FOR_YMAX)
    ax_cov.set_ylim(0, max(ymax * 1.1, cov_plot["cov"].max() * 1.02, contig_mean * 1.2))
    ax_cov.set_ylabel("Coverage (×)")
    ax_cov.grid(True, alpha=0.25)
    ax_cov.legend(loc="upper right", fontsize=8, frameon=False)

    # Bottom bars
    ax_bar.add_patch(Rectangle((1, 0.25), max(0, info.contig_length - 1), 0.5,
                               facecolor="#DDDDDD", edgecolor="black", lw=1.0))
    ax_bar.add_patch(Rectangle((info.start_1based, 0.15),
                               max(1, info.end_1based - info.start_1based + 1), 0.7,
                               facecolor="red", edgecolor="black", lw=1.0, alpha=0.7))

    # NEW: label the prophage length centered over the red bar
    x_mid = info.start_1based + (info.prophage_length - 1) / 2.0
    ax_bar.text(
        x_mid, 0.90,
        f"Prophage length: {format_bp(info.prophage_length)}",
        ha="center", va="center", fontsize=9, color="red",
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1.5),
        clip_on=False
    )

    ax_bar.set_ylim(0, 1)
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("Position on contig (bp)")
    ax_bar.set_xlim(1, info.contig_length)

    # Titles
    bits = [f"{info.contig_accession} [{info.start_1based:,}-{info.end_1based:,}]"]
    if info.srr_id: bits.append(info.srr_id)
    if info.family: bits.append(info.family)
    if info.origin: bits.append(info.origin)
    plt.suptitle("   ·   ".join(bits), fontsize=12)

    if info.tax_lineage:
        ax_cov.set_title(info.tax_lineage, fontsize=9, loc="left")

    # Pretty x formatting
    ax_bar.xaxis.set_major_formatter(FuncFormatter(thousands))
    ax_cov.xaxis.set_major_formatter(FuncFormatter(thousands))

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    return contig_mean, prophage_mean


def build_prophage_info(prophage_dir: Path, annot: pd.DataFrame, basecov_path: Path) -> ProphageInfo:
    contig_acc, s1, e1 = parse_dir_name(prophage_dir.name)
    row = find_annotation_row_for_dir(annot, contig_acc, s1, e1)
    if row is None:
        raise RuntimeError(
            f"No annotation match for {prophage_dir.name} "
            f"(bacterial_id={contig_acc}, start={s1}, end={e1})."
        )

    contig_len = int(row["bacterial_contig_length"])
    # Prefer the table's prophage_length if present and valid; else compute.
    if "prophage_length" in row.index and pd.notna(row["prophage_length"]):
        try:
            proph_len = int(row["prophage_length"])
        except Exception:
            proph_len = e1 - s1 + 1
    else:
        proph_len = e1 - s1 + 1

    family = str(row["crassvirales_family"]) if "crassvirales_family" in row.index else None
    origin = str(row["genome_origin"]) if "genome_origin" in row.index else None
    tax   = str(row["bacterial_taxonomy_lineage"]) if "bacterial_taxonomy_lineage" in row.index else None
    # Normalize empties
    family = None if family in ("nan", "", "None") else family
    origin = None if origin in ("nan", "", "None") else origin
    tax    = None if tax    in ("nan", "", "None") else tax

    srr = parse_srr_from_filename(basecov_path)

    return ProphageInfo(
        dir_name=prophage_dir.name,
        contig_accession=contig_acc,
        start_1based=s1,
        end_1based=e1,
        contig_length=contig_len,
        prophage_length=proph_len,       # <-- NEW
        srr_id=srr,
        family=family,
        origin=origin,
        tax_lineage=tax,
    )


def main() -> None:
    args = parse_args()
    wd: Path = args.wd
    figdir: Path = args.figdir
    annot_path: Path = args.annotation if args.annotation.is_absolute() else (wd / args.annotation)

    if not wd.exists():
        raise SystemExit(f"[ERR] Working directory not found: {wd}")
    if not annot_path.exists():
        raise SystemExit(f"[ERR] Annotation table not found: {annot_path}")

    annot = load_annotation_table(annot_path)

    made_any = False
    for prophage_dir in iter_prophage_dirs(wd):
        basecov = find_basecov_file(prophage_dir)
        if basecov is None:
            print(f"[SKIP] No basecov*.txt in {prophage_dir.name}")
            continue

        try:
            info = build_prophage_info(prophage_dir, annot, basecov)
        except Exception as e:
            print(f"[WARN] {prophage_dir.name}: {e}")
            continue

        try:
            cov_df = load_basecov_filtered(basecov, info.contig_accession)
        except Exception as e:
            print(f"[WARN] {prophage_dir.name}: {e}")
            continue

        out_png = figdir / f"{prophage_dir.name}_basecov.png"
        contig_mean, prophage_mean = plot_prophage_coverage(info, cov_df, out_png, args.downsample_target)
        print(f"[OK] Saved: {out_png}")
        print(f"     Means — Contig: {contig_mean:.2f}×   Prophage: {prophage_mean:.2f}×")
        made_any = True

    if not made_any:
        print("[INFO] No figures produced. Check folder names and basecov files.")


if __name__ == "__main__":
    main()
