#!/usr/bin/env python3
"""
Create protein-level taxonomic annotation table from a phylogenetic tree
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make protein-level annotation table for tree leaves."
    )
    parser.add_argument(
        "--tree",
        type=Path,
        required=True,
        help="Newick file with protein leaf labels.",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        required=True,
        help="TSV: genome, family, subfamily, genus, species.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output TSV.",
    )
    return parser.parse_args()


def load_taxonomy_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"genome", "family", "subfamily", "genus", "species"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Taxonomy table {path} missing columns: {missing}")
    return df


# Regex: split {genome_id}_{protein_id}, where protein_id ends with "_CDS_number"
LEAF_REGEX = re.compile(r"^(?P<genome>.+?)_(?P<prot>[A-Za-z0-9]+_CDS_\d+)$")


def split_leaf_label(label: str) -> Optional[Tuple[str, str]]:
    """
    Return (genome_id, protein_id) or None if cannot parse.
    """
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
    """
    Build annotation for all parseable leaves whose genome_id is in taxonomy.
    """
    tax_map: Dict[str, pd.Series] = {
        row["genome"]: row for _, row in taxonomy_df.iterrows()
    }

    rows: List[Dict[str, str]] = []

    for label in leaf_labels:
        parsed = split_leaf_label(label)
        if parsed is None:
            continue  # skip silently
        genome_id, protein_id = parsed

        tax_row = tax_map.get(genome_id)
        if tax_row is None:
            continue  # skip silently

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


def main() -> None:
    args = parse_args()

    taxonomy_df = load_taxonomy_table(args.taxonomy)
    leaf_labels = extract_leaf_labels(args.tree)

    df = build_protein_annotation(leaf_labels, taxonomy_df)
    df.to_csv(args.out, sep="\t", index=False)

    # ---- Summary ----
    n_proteins = df["protein_id"].nunique()
    n_genomes = df["genome_id"].nunique()

    print(f"[INFO] Wrote {len(df)} annotated rows to: {args.out}")
    print(f"[INFO] Unique proteins annotated: {n_proteins}")
    print(f"[INFO] Unique genomes annotated:  {n_genomes}")


if __name__ == "__main__":
    main()
