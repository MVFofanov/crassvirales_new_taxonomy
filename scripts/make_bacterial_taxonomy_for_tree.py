#!/usr/bin/env python3
"""
Build a bacterial taxonomy + geNomad + bacterial length annotation table
for all leaves in a TerL tree.

Output columns:
  leaf_label
  seq_name
  accession
  taxonomy
  domain
  kingdom
  phylum
  class
  order
  family
  genus
  species
  genomad_length
  topology
  bacterial_length
  prophage_ratio   (genomad_length / bacterial_length, rounded to 2 d.p., or "missed")

Only bacterial genomes are kept.
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
        description=(
            "Make bacterial + geNomad + bacterial length annotation table "
            "for TerL tree leaves."
        )
    )
    parser.add_argument("--tree", type=Path, required=True, help="TerL Newick treefile")
    parser.add_argument(
        "--taxonomy",
        type=Path,
        required=True,
        help="TSV with bacterial taxonomy: accession, taxonomy",
    )
    parser.add_argument(
        "--genomad_summary",
        type=Path,
        required=True,
        help="geNomad summary TSV with columns: seq_name, length, topology",
    )
    parser.add_argument(
        "--bacterial_lengths",
        type=Path,
        required=True,
        help="TSV with bacterial contig lengths: contig, length",
    )
    parser.add_argument("--out", type=Path, required=True, help="Output TSV")
    return parser.parse_args()


def load_taxonomy(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"accession", "taxonomy"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Taxonomy file {path} missing columns: {missing}")
    return df


def load_genomad_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"seq_name", "length", "topology"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"geNomad summary {path} missing columns: {missing}")
    return df


def load_bacterial_lengths(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"contig", "length"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Bacterial length file {path} missing columns: {missing}")
    return df


# Regex to capture nucleotide accessions like NC_020054.1 / DAOMAL010000008.1
ACCESSION_REGEX = re.compile(r"^([A-Za-z0-9_]+?\d\.\d+)")


def extract_accession(label: str) -> Optional[str]:
    m = ACCESSION_REGEX.match(label)
    return m.group(1) if m else None


def extract_seq_name_from_leaf(label: str) -> str:
    """
    Example:
      NC_020054.1|provirus_3988657_4054034_BVMJILSW_CDS_9409
        →  NC_020054.1|provirus_3988657_4054034

      BAAHBE010000021.1_BVMJILSW_CDS_1234
        →  BAAHBE010000021.1
    """
    if "_BVMJILSW_CDS_" in label:
        return label.split("_BVMJILSW_CDS_")[0]

    m = re.match(r"^(.*)_CDS_\d+$", label)
    if m:
        return m.group(1)

    return label


def parse_taxonomy_string(
    tax_str: str,
) -> Tuple[str, str, str, str, str, str, str, str]:
    """
    Returns:
      domain, kingdom, phylum, class, order, family, genus, species

    If only domain is known (e.g. "Bacteria"), all ranks inherit domain.
    """
    if not isinstance(tax_str, str) or not tax_str.strip():
        return ("unknown",) * 8

    parts = [p.strip() for p in tax_str.split(";")]

    domain = parts[0] if parts else "unknown"
    kingdom = parts[1].strip() if len(parts) > 1 and parts[1].strip() else ""

    # Expected: phylum, class, order, family, genus, species
    raw: List[str] = []
    for idx in range(2, 8):
        raw.append(parts[idx].strip() if idx < len(parts) else "")

    last_known = ""
    filled: List[str] = []
    for r in raw:
        if r:
            last_known = r
            filled.append(r)
        else:
            filled.append(last_known)

    # If nothing below kingdom is known:
    if not kingdom and not last_known:
        if domain != "unknown":
            # everything inherits domain
            return (domain,) * 8
        else:
            return ("unknown",) * 8

    if not kingdom:
        kingdom = domain

    phylum, klass, order, family, genus, species = filled
    return domain, kingdom, phylum, klass, order, family, genus, species


def safe_prophage_ratio(genomad_length: str, bacterial_length: str) -> str:
    """
    Compute prophage_ratio = round(genomad_length / bacterial_length, 2).

    If any value is missing / non-numeric / zero → "missed".
    """
    try:
        gl = float(genomad_length)
        bl = float(bacterial_length)
        if bl <= 0:
            return "missed"
        ratio = gl / bl
        return f"{ratio:.2f}"
    except (TypeError, ValueError):
        return "missed"


def build_leaf_taxonomy(
    tree_path: Path,
    taxonomy_df: pd.DataFrame,
    genomad_df: pd.DataFrame,
    bacterial_len_df: pd.DataFrame,
) -> pd.DataFrame:

    tax_map: Dict[str, str] = dict(zip(taxonomy_df["accession"], taxonomy_df["taxonomy"]))

    # geNomad mapping: seq_name -> (genomad_length, topology)
    geno_map: Dict[str, Tuple[str, str]] = {
        row["seq_name"]: (row["length"], row["topology"])
        for _, row in genomad_df.iterrows()
    }

    # bacterial contig length mapping: accession/contig -> length
    bact_len_map: Dict[str, str] = dict(zip(bacterial_len_df["contig"], bacterial_len_df["length"]))

    t = Tree(tree_path.read_text(), format=1)

    rows: List[Dict[str, str]] = []
    total = 0
    kept = 0

    for leaf in t.iter_leaves():
        total += 1
        label = leaf.name
        acc = extract_accession(label)

        if acc is None:
            continue

        taxonomy_str = tax_map.get(acc)
        if taxonomy_str is None:
            continue

        domain, kingdom, phylum, klass, order, family, genus, species = \
            parse_taxonomy_string(taxonomy_str)

        # keep only bacterial genomes
        if domain != "Bacteria":
            continue

        seq_name = extract_seq_name_from_leaf(label)
        genomad_length, topo_val = geno_map.get(seq_name, ("unknown", "unknown"))

        bacterial_length = bact_len_map.get(acc, "unknown")
        prophage_ratio = (
            safe_prophage_ratio(genomad_length, bacterial_length)
            if genomad_length != "unknown" and bacterial_length != "unknown"
            else "missed"
        )

        kept += 1
        rows.append(
            {
                "leaf_label": label,
                "seq_name": seq_name,
                "accession": acc,
                "taxonomy": taxonomy_str,
                "domain": domain,
                "kingdom": kingdom,
                "phylum": phylum,
                "class": klass,
                "order": order,
                "family": family,
                "genus": genus,
                "species": species,
                "genomad_length": genomad_length,
                "topology": topo_val,
                "bacterial_length": bacterial_length,
                "prophage_ratio": prophage_ratio,
            }
        )

    print(f"[INFO] Total leaves: {total}")
    print(f"[INFO] Bacterial leaves kept: {kept}")
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()

    taxonomy_df = load_taxonomy(args.taxonomy)
    genomad_df = load_genomad_summary(args.genomad_summary)
    bacterial_len_df = load_bacterial_lengths(args.bacterial_lengths)

    df = build_leaf_taxonomy(args.tree, taxonomy_df, genomad_df, bacterial_len_df)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"[INFO] Wrote {len(df)} rows → {args.out}")
    print(f"[INFO] Unique bacterial accessions: {df['accession'].nunique()}")

    print("\n[INFO] Phylum counts:")
    print(df["phylum"].value_counts().sort_values(ascending=False))

    print("\n[INFO] Class counts:")
    print(df["class"].value_counts().sort_values(ascending=False))


if __name__ == "__main__":
    main()
