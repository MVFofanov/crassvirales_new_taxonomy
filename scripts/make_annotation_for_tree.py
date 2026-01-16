#!/usr/bin/env python3
"""
Create protein-level taxonomic annotation table from one or more phylogenetic trees
and a genome-level taxonomy table.

Only proteins whose leaf labels can be parsed AND whose genome_id
exists in the taxonomy table are included. Others are silently skipped.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from ete3 import Tree


# Regex: split {genome_id}_{protein_id}, where protein_id ends with "_CDS_number"
LEAF_REGEX = re.compile(r"^(?P<genome>.+?)_(?P<prot>[A-Za-z0-9]+_CDS_\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make protein-level annotation table for tree leaves (one or more trees)."
    )
    parser.add_argument(
        "--tree",
        type=Path,
        required=True,
        nargs="+",
        help="One or more Newick files with protein leaf labels.",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        required=True,
        help="TSV with columns: genome, family, subfamily, genus, species.",
    )

    out_group = parser.add_mutually_exclusive_group(required=True)
    out_group.add_argument(
        "--out",
        type=Path,
        nargs="+",
        help="One or more output TSV paths. If given, must match number of --tree files.",
    )
    out_group.add_argument(
        "--outdir",
        type=Path,
        help="Output directory. Output names are derived from tree filenames.",
    )

    parser.add_argument(
        "--suffix",
        type=str,
        default="_tree_protein_taxonomy.tsv",
        help="Suffix for outputs when using --outdir (default: _tree_protein_taxonomy.tsv).",
    )
    return parser.parse_args()


def load_taxonomy_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"genome", "family", "subfamily", "genus", "species"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Taxonomy table {path} missing columns: {missing}")
    return df


def split_leaf_label(label: str) -> Optional[Tuple[str, str]]:
    """Return (genome_id, protein_id) or None if cannot parse."""
    m = LEAF_REGEX.match(label)
    if not m:
        return None
    return m.group("genome"), m.group("prot")


def extract_leaf_labels(tree_path: Path) -> List[str]:
    t = Tree(tree_path.read_text(), format=1)
    return [leaf.name for leaf in t.iter_leaves()]


def build_protein_annotation(
    leaf_labels: List[str], taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """Build annotation for all parseable leaves whose genome_id is in taxonomy."""
    # Faster + simpler than iterrows->Series in a dict:
    tax_map: Dict[str, Dict[str, str]] = (
        taxonomy_df.set_index("genome")[["family", "subfamily", "genus", "species"]]
        .to_dict(orient="index")
    )

    rows: List[Dict[str, str]] = []

    for label in leaf_labels:
        parsed = split_leaf_label(label)
        if parsed is None:
            continue
        genome_id, protein_id = parsed

        tax_row = tax_map.get(genome_id)
        if tax_row is None:
            continue

        rows.append(
            {
                "leaf_label": label,
                "genome_id": genome_id,
                "protein_id": protein_id,
                "family": tax_row["family"],
                "subfamily": tax_row["subfamily"],
                "genus": tax_row["genus"],
                "species": tax_row["species"],
            }
        )

    return pd.DataFrame(rows)


def resolve_outputs(trees: List[Path], out: Optional[List[Path]], outdir: Optional[Path], suffix: str) -> List[Path]:
    if out is not None:
        if len(out) != len(trees):
            raise SystemExit(
                f"[ERROR] You provided {len(trees)} --tree files but {len(out)} --out files. They must match."
            )
        return out

    assert outdir is not None
    outdir.mkdir(parents=True, exist_ok=True)

    outs: List[Path] = []
    for tp in trees:
        # e.g. TerL_gappy0.5.treefile -> TerL_gappy0.5_tree_protein_taxonomy.tsv
        outs.append(outdir / f"{tp.stem}{suffix}")
    return outs


def main() -> None:
    args = parse_args()

    taxonomy_df = load_taxonomy_table(args.taxonomy)
    out_paths = resolve_outputs(args.tree, args.out, args.outdir, args.suffix)

    for tree_path, out_path in zip(args.tree, out_paths):
        leaf_labels = extract_leaf_labels(tree_path)
        df = build_protein_annotation(leaf_labels, taxonomy_df)
        df.to_csv(out_path, sep="\t", index=False)

        n_proteins = df["protein_id"].nunique() if not df.empty else 0
        n_genomes = df["genome_id"].nunique() if not df.empty else 0

        print(f"[INFO] Tree: {tree_path}")
        print(f"[INFO] Wrote {len(df)} annotated rows to: {out_path}")
        print(f"[INFO] Unique proteins annotated: {n_proteins}")
        print(f"[INFO] Unique genomes annotated:  {n_genomes}")
        print("")


if __name__ == "__main__":
    main()
