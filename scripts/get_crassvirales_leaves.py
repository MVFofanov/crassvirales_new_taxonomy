#!/usr/bin/env python3

import argparse
from pathlib import Path

from ete3 import Tree


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate phylogenetic tree using iTOL-style file and extract unannotated Crassvirales leaves."
    )
    parser.add_argument("--tree_file", type=Path, required=True, help="Newick format tree file")
    parser.add_argument("--itol_annotation", type=Path, required=True, help="iTOL DATASET_STYLE annotation file")
    parser.add_argument(
        "--crassvirales_output", type=Path, required=True, help="Output .txt file with unannotated Crassvirales leaves"
    )
    parser.add_argument(
        "--genome_output",
        type=Path,
        required=True,
        help="Output .txt file with genome IDs of unannotated Crassvirales leaves",
    )
    return parser.parse_args()


def parse_itol_annotation(file_path: Path) -> tuple[dict[str, str], set[str]]:
    leaf_to_family: dict[str, str] = {}
    crassvirales_families: set[str] = set()
    color_to_family: dict[str, str] = {}

    legend_labels = []
    legend_colors = []

    in_data_section = False
    with file_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("LEGEND_LABELS"):
                legend_labels = line.split("\t")[1:]
            elif line.startswith("LEGEND_COLORS"):
                legend_colors = line.split("\t")[1:]
            elif line.startswith("DATA"):
                in_data_section = True
                continue

            if not in_data_section:
                continue

            # Parse DATA section
            parts = line.split("\t")
            if len(parts) >= 4 and parts[1] == "label" and parts[2] == "node":
                leaf_name = parts[0]
                color = parts[3]
                if color in color_to_family:
                    family = color_to_family[color]
                    leaf_to_family[leaf_name] = family
                else:
                    # fallback in case LEGEND mapping isn't parsed yet
                    leaf_to_family[leaf_name] = None
            elif len(parts) == 4 and parts[1] == "branch" and parts[2] == "clade":
                crassvirales_families.add(parts[0])

    # Build color → family map
    for color, label in zip(legend_colors, legend_labels, strict=False):
        if label not in {"NA", "outgroup"}:
            color_to_family[color] = label

    # Remap leaf_to_family now that we have color mapping
    for leaf, fam in list(leaf_to_family.items()):
        if fam is None and leaf in leaf_to_family:
            # Try to resolve based on final mapping
            color = None
            # Find the color again (inefficient fallback)
            with file_path.open() as f:
                for line in f:
                    if line.startswith(leaf + "\t"):
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            color = parts[3]
                            break
            if color in color_to_family:
                leaf_to_family[leaf] = color_to_family[color]

    return leaf_to_family, crassvirales_families


def annotate_tree(tree: Tree, leaf_to_family: dict[str, str]) -> None:
    for leaf in tree.iter_leaves():
        if leaf.name in leaf_to_family and leaf_to_family[leaf.name] is not None:
            leaf.add_feature("family", leaf_to_family[leaf.name])


def find_crassvirales_mrca(tree: Tree, crassvirales_families: set[str]) -> Tree:
    crass_leaves = [
        leaf for leaf in tree.iter_leaves() if any(leaf.name.startswith(fam) for fam in crassvirales_families)
    ]
    if not crass_leaves:
        raise ValueError("No Crassvirales leaves found in the tree")
    return tree.get_common_ancestor(crass_leaves)


def find_unannotated_leaves(mrca_node: Tree, leaf_to_family: dict[str, str]) -> list[str]:
    return sorted([leaf.name for leaf in mrca_node.iter_leaves() if leaf.name not in leaf_to_family])


def main() -> None:
    args = parse_args()
    tree = Tree(str(args.tree_file), format=1)
    leaf_to_family, crassvirales_families = parse_itol_annotation(args.itol_annotation)
    crass_mrca = find_crassvirales_mrca(tree, crassvirales_families)
    unannotated = find_unannotated_leaves(crass_mrca, leaf_to_family)

    # Save unannotated leaf names
    with args.crassvirales_output.open("w") as out:
        for name in unannotated:
            out.write(name + "\n")

    # Extract unique genome IDs from leaf names
    genome_ids = sorted(set("|".join(name.split("|")[:-2]) for name in unannotated))

    with args.genome_output.open("w") as out:
        for genome_id in genome_ids:
            out.write(genome_id + "\n")

    print(f"Saved {len(unannotated)} unannotated Crassvirales leaves to {args.crassvirales_output}")
    print(f"Saved {len(genome_ids)} unique genome IDs to {args.genome_output}")


if __name__ == "__main__":
    main()
