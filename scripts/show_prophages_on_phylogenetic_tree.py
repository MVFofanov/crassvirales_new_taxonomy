#!/usr/bin/env python3

from collections import OrderedDict
import os
from pathlib import Path
from typing import Dict, Tuple, Set, Optional
from ete3 import Tree, TreeStyle, NodeStyle, RectFace, faces, TextFace, Face
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use('Agg')
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering


def read_tree(tree_file: Path) -> Tree:
    return Tree(str(tree_file), format=1)

def plot_prophages(df: pd.DataFrame, output_file: str):
    """
    Plot prophage coordinates per genome using matplotlib.

    df must have: 'genome', 'contig_id', 'start', 'end', 'contig_length'
    """
    df = df.copy()
    df["genome"] = df["contig_id"].str.split("|").str[0]
    df["length"] = df["end"].astype(int) - df["start"].astype(int)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["contig_length"] = df["contig_length"].astype(int)

    fig, ax = plt.subplots(figsize=(10, 0.4 * len(df["genome"].unique())))

    y_ticks = []
    y_labels = []
    for i, (genome, group) in enumerate(df.groupby("genome")):
        y_ticks.append(i)
        y_labels.append(genome)

        for _, row in group.iterrows():
            ax.hlines(y=i, xmin=0, xmax=row["contig_length"], color="lightgray", linewidth=6)
            ax.hlines(y=i, xmin=row["start"], xmax=row["end"], color="red", linewidth=6)

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Genomic position")
    ax.set_title("Prophage positions in contigs")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_prophage_positions(
    tree: Tree,
    prophage_dict: Dict[str, Dict[str, str]],
    leaf_to_family: Dict[str, str],
    family_to_color: Dict[str, str],
    output_file: str
) -> None:
    """
    Plot prophage positions ordered by tree leaves with Crassvirales family color box on the right.
    """
    from collections import OrderedDict

    # Get leaf order from tree
    leaf_order = [leaf.name for leaf in tree.iter_leaves()]

    # Reorder prophage_dict based on leaf order
    ordered_prophage_dict = OrderedDict()
    for leaf_name in leaf_order:
        genome_id = match_genome_id(leaf_name)
        for contig_id, prophage in prophage_dict.items():
            if prophage.get("prophage_id", "").startswith(genome_id):
                ordered_prophage_dict[contig_id] = prophage

    fig, ax = plt.subplots(figsize=(10, len(ordered_prophage_dict) * 0.4))
    y_labels = []
    y_positions = list(range(len(ordered_prophage_dict)))

    for i, (contig_id, prophage_info) in enumerate(ordered_prophage_dict.items()):
        start = int(prophage_info["start"])
        end = int(prophage_info["end"])
        contig_length = int(prophage_info["contig_length"])
        prophage_id = prophage_info["prophage_id"]
        genome_id = prophage_id.split("|")[0]

        # Label: contig and genome
        label = f"{contig_id}\n{genome_id}"
        y_labels.append(label)

        # Contig bar
        ax.broken_barh([(0, contig_length)], (i - 0.4, 0.8), facecolors='lightgrey')
        ax.broken_barh([(start, end - start)], (i - 0.4, 0.8), facecolors='red')

        # Family color box (right side)
        family = leaf_to_family.get(prophage_id)
        if family:
            color = family_to_color.get(family, "#999999")
            # Draw a small colored bar at the far right
            ax.broken_barh([(contig_length + 1000, 50000)], (i - 0.3, 0.6), facecolors=color)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=6)
    ax.set_xlabel("Genomic position")
    ax.set_title("Prophage positions in contigs (ordered by tree)")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()



def parse_itol_annotation(file_path: Path) -> Tuple[Dict[str, str], Dict[str, str], Set[str]]:
    leaf_to_family: Dict[str, str] = {}
    color_to_family: Dict[str, str] = {}
    family_to_color: Dict[str, str] = {}
    crassvirales_families: Set[str] = set()

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
            parts = line.split("\t")
            if len(parts) >= 4 and parts[1] == "label" and parts[2] == "node":
                leaf_to_family[parts[0]] = None  # placeholder
            elif len(parts) == 4 and parts[1] == "branch" and parts[2] == "clade":
                crassvirales_families.add(parts[0])

    for color, label in zip(legend_colors, legend_labels):
        if label not in {"NA", "outgroup"}:
            color_to_family[color] = label
            family_to_color[label] = color

    # Second pass: assign families to leaves based on color
    with file_path.open() as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and parts[1] == "label" and parts[2] == "node":
                leaf, color = parts[0], parts[3]
                if leaf in leaf_to_family and color in color_to_family:
                    leaf_to_family[leaf] = color_to_family[color]

    return leaf_to_family, family_to_color, crassvirales_families


def parse_prophage_table(path: Path) -> Tuple[Dict[str, Dict[str, str]], pd.DataFrame]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    prophage_dict: Dict[str, Dict[str, str]] = {}
    for _, row in df.iterrows():
        contig_id = row["contig_id"]
        prophage_dict[contig_id] = row.to_dict()  # 👈 ключ = contig_id
    return prophage_dict, df


def match_genome_id(leaf_name: str) -> str:
    return leaf_name.split("|")[0]


def annotate_tree_features(
    tree: Tree,
    leaf_to_family: Dict[str, str],
    family_to_color: Dict[str, str],
    prophage_dict: Dict[str, Dict[str, str]]
) -> None:
    matched = 0

    for leaf in tree.iter_leaves():
        matched_prophage = None
        for prophage in prophage_dict.values():
            prophage_id = prophage.get("prophage_id", "")
            if leaf.name.startswith(prophage_id):
                matched_prophage = prophage
                break

        family = leaf_to_family.get(leaf.name)
        leaf.add_features(family=family, is_prophage=matched_prophage is not None)

        if matched_prophage:
            leaf.add_features(
                start=int(matched_prophage["start"]),
                end=int(matched_prophage["end"]),
                contig_length=int(matched_prophage["contig_length"])
            )
            matched += 1
            print(f"[MATCHED] {leaf.name} ← {matched_prophage['prophage_id']}")

    print(f"[INFO] Matched {matched} prophage leaves out of {len(tree.get_leaves())}")


# class ProphageBarFace(Face):
#     def __init__(self, contig_len: int, start: int, end: int, width: int = 100, height: int = 10):
#         Face.__init__(self)
#         self.contig_len = contig_len
#         self.start = start
#         self.end = end
#         self.width = width
#         self.height = height
#         self.margin_left = 2
#         self.margin_right = 2
#         self.margin_top = 2
#         self.margin_bottom = 2

#     def update_pixels(self, node):
#         from PyQt5.QtGui import QPainter, QColor
#         self.image = self._get_image(self.width, self.height)
#         qp = QPainter(self.image)
#         qp.setRenderHint(QPainter.Antialiasing)

#         # Draw black background (full contig)
#         qp.setBrush(QColor("black"))
#         qp.drawRect(0, 0, self.width, self.height)

#         # Draw red region (prophage)
#         rel_start = int((self.start / self.contig_len) * self.width)
#         rel_end = int((self.end / self.contig_len) * self.width)
#         width = max(1, rel_end - rel_start)

#         qp.setBrush(QColor("red"))
#         qp.drawRect(rel_start, 0, width, self.height)

#         print(f"[DEBUG] Rendering prophage bar | start: {self.start}, end: {self.end}, contig_len: {self.contig_len}")

#         qp.end()

class ProphageBarFace(Face):
    def __init__(self, contig_len: int, start: int, end: int, width: int = 100, height: int = 10):
        super().__init__()
        self.contig_len = contig_len
        self.start = start
        self.end = end
        self.width = width
        self.height = height
        self.margin_left = 2
        self.margin_right = 2
        self.margin_top = 2
        self.margin_bottom = 2

    def update_pixels(self, node):
        from PyQt5.QtGui import QPainter, QColor
        self.image = self._get_image(self.width, self.height)
        qp = QPainter(self.image)

        # Draw black background (entire contig)
        qp.setBrush(QColor("black"))
        qp.drawRect(0, 0, self.width, self.height)

        # Draw red region (prophage)
        rel_start = int((self.start / self.contig_len) * self.width)
        rel_end = int((self.end / self.contig_len) * self.width)
        red_width = max(1, rel_end - rel_start)

        qp.setBrush(QColor("red"))
        qp.drawRect(rel_start, 0, red_width, self.height)

        qp.end()


def get_custom_layout(family_to_color: Dict[str, str]):
    def custom_layout(node):
        if node.is_leaf():
            # Family color box
            if hasattr(node, "family") and node.family:
                color = family_to_color.get(node.family, "#AAAAAA")
                node.add_face(RectFace(10, 10, color, color), column=0, position="aligned")

            # Draw prophage bar using standard RectFaces
            if getattr(node, "is_prophage", False):
                try:
                    contig_len = node.contig_length
                    start = node.start
                    end = node.end

                    # Single overlay bar
                    node.add_face(ProphageBarFace(contig_len, start, end, width=100, height=10), column=1, position="aligned")

                except Exception as e:
                    print(f"Failed to draw prophage for {node.name}: {e}")
    return custom_layout



def render_tree(tree: Tree, output_file: Path, family_to_color: Dict[str, str]) -> None:
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.force_topology = True
    ts.scale = 120
    ts.title.add_face(TextFace("TerL Tree with Prophage Annotations", fsize=12), column=0)
    ts.layout_fn = get_custom_layout(family_to_color)  # use closure
    tree.render(str(output_file), w=2000, tree_style=ts)

def main():
    wd = "/mnt/c/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages/blast_prophages_vs_ncbi_and_gtdb"
    crassus_phylogeny = "/mnt/c/crassvirales/CrassUS_old/CrassUS/results/prophage_analysis/5_phylogenies/3_iToL"
    tree_file = Path(f"{crassus_phylogeny}/TerL_iToL_renamed.nwk")
    itol_annotation = Path(f"{crassus_phylogeny}/TerL_iToL.txt")
    prophage_table = Path(f"{wd}/prophage_alignments_ncbi_and_gtdb_samtools_coordinates_edited_with_lengths.tsv")

    output_svg = Path(f"{wd}/TerL_with_prophages_annotated_tree.svg")
    prophage_plot = Path(f"{wd}/prophage_positions_plot.png")


    tree = read_tree(tree_file)
    leaf_to_family, family_to_color, _ = parse_itol_annotation(itol_annotation)

    prophage_dict, prophage_df = parse_prophage_table(prophage_table)

    # leaf_order = [leaf.name for leaf in tree.iter_leaves()]
    # ordered_prophage_dict = OrderedDict()
    # for leaf_name in leaf_order:
    #     genome_id = match_genome_id(leaf_name)
    #     for contig_id, prophage in prophage_dict.items():
    #         prophage_id = prophage.get("prophage_id", "")
    #         if prophage_id.startswith(genome_id):
    #             ordered_prophage_dict[contig_id] = prophage


    # prophage_dict = ordered_prophage_dict
    print(prophage_dict)

    leaf_ids = {match_genome_id(leaf.name) for leaf in tree.iter_leaves()}
    # matched_ids = set(prophage_dict.keys()) & leaf_ids
    # print(f"[INFO] Matching prophage IDs: {len(matched_ids)} out of {len(prophage_dict)} total prophages")

    matched = sum(
        any(leaf.name.startswith(prophage.get("prophage_id", "")) for leaf in tree.iter_leaves())
        for prophage in prophage_dict.values()
    )
    print(f"[INFO] Matching prophage IDs: {matched} out of {len(prophage_dict)} total prophages")


    annotate_tree_features(tree, leaf_to_family, family_to_color, prophage_dict)
    render_tree(tree, output_svg, family_to_color)
    print(f"Annotated tree saved to: {output_svg}")

    # New: plot genomic positions
    # plot_prophages(prophage_df, output_file=str(prophage_plot))
    # plot_prophage_positions(prophage_dict, output_file=str(prophage_plot))

    plot_prophage_positions(
        tree,
        prophage_dict,
        leaf_to_family,
        family_to_color,
        output_file=str(prophage_plot)
    )

    print(f"Prophage barplot saved to: {prophage_plot}")


if __name__ == "__main__":
    main()
