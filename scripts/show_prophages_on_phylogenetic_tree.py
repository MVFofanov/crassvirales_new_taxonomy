#!/usr/bin/env python3

from collections import OrderedDict
import os
from pathlib import Path
import re
from typing import Dict, Tuple, List, Set, Optional
from ete3 import Tree, TreeStyle, NodeStyle, RectFace, faces, TextFace, Face, NCBITaxa
import matplotlib
from matplotlib import gridspec
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

from PyQt5.QtGui import QPainter, QColor

from dna_features_viewer import GraphicFeature, GraphicRecord

import textwrap

matplotlib.use('Agg')
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering

ORDER_LIKE = re.compile(r".*ales$", re.IGNORECASE)

# top-level (after imports)
try:
    _NCBI = NCBITaxa()
except Exception:
    _NCBI = None


def read_tree(tree_file: Path) -> Tree:
    return Tree(str(tree_file), format=1)

def load_functional_annotation(path: str) -> Dict[str, List[GraphicFeature]]:
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    features_by_genome: Dict[str, List[GraphicFeature]] = {}

    allowed = {"TerL", "portal", "MCP", "Ttub", "Tstab", "primase", "PolB", "PolA", "DnaB"}

    for _, row in df.iterrows():
        genome = str(row["genome"])
        start = int(row["start"])
        end = int(row["end"])
        strand = 1 if str(row["strand"]).strip() == "+" else -1

        yutin_raw = row.get("yutin", None)
        yutin = (str(yutin_raw).strip() if pd.notna(yutin_raw) else "")

        # color: red only if yutin present, else lightgrey
        color = "red" if yutin else "lightgrey"
        # label: only for allowed markers
        label = yutin if yutin in allowed else None

        features_by_genome.setdefault(genome, []).append(
            GraphicFeature(
                start=start,
                end=end,
                strand=strand,
                color=color,
                label=label
            )
        )
    return features_by_genome


def _split_taxonomy_string(tax: str) -> List[str]:
    if not isinstance(tax, str) or not tax.strip():
        return []
    toks = [t.strip() for t in tax.split(";") if t.strip()]
    # drop “environmental samples”
    return [t for t in toks if t.lower() != "environmental samples"]

def _try_ncbi_order(tokens: List[str]) -> Optional[str]:
    if _NCBI is None:
        return None
    try:
        name2tax = _NCBI.get_name_translator(tokens)  # {name: [taxid]}
        if not name2tax:
            return None
        # build ranks in one go (fewer DB hits)
        taxids = [name2tax[t][0] for t in tokens if t in name2tax]
        ranks = _NCBI.get_rank(taxids)  # {taxid: rank}
        rev = {v: k for k, v in name2tax.items()}  # not reliable; better iterate
        for tok in tokens:
            tids = name2tax.get(tok)
            if not tids:
                continue
            taxid = tids[0]
            if ranks.get(taxid) == "order":
                return tok
    except Exception:
        return None
    return None

def _heuristic_order(tokens: List[str]) -> Optional[str]:
    # pick first token that looks like an order (…ales)
    for tok in tokens:
        if ORDER_LIKE.match(tok):
            return tok
    return None

def _pick_suffix(tokens: List[str], suffixes: Tuple[str, ...]) -> Optional[str]:
    for tok in tokens:
        if tok.endswith(suffixes):
            return tok
    return None

def extract_order_with_fallback(taxonomy: str) -> Tuple[str, str]:
    """
    Returns (used_rank, name) where used_rank ∈ {'order','class','phylum','family','unknown'}
    Priority: NCBITaxa (order) → heuristic order → class → phylum → family → Unknown
    """
    tokens = _split_taxonomy_string(taxonomy)
    if not tokens:
        return "unknown", "Unknown"

    # 1) NCBITaxa first
    ord_name = _try_ncbi_order(tokens)
    if ord_name:
        return "order", ord_name

    # 2) Heuristic (…ales)
    ord_name = _heuristic_order(tokens)
    if ord_name:
        return "order", ord_name

    # 3) Fallbacks by common suffixes
    cls = _pick_suffix(tokens, ("ia", "bia", "illi", "idia", "eria"))  # includes Bacilli       # Clostridia, Halobacteria
    if cls: return "class", cls
    phyl = _pick_suffix(tokens, ("ota", "otae"))      # Bacillota, Bacteroidota
    if phyl: return "phylum", phyl
    fam = _pick_suffix(tokens, ("aceae",))            # Lachnospiraceae
    if fam: return "family", fam

    return "unknown", "Unknown"

# def build_order_color_map(unique_orders: list) -> dict:
#     """Assign distinct colors to orders, keeping 'unknown' grey."""
#     import itertools
#     color_cycle = itertools.cycle(plt.cm.tab20.colors)
#     order_to_color = {}
#     for order in sorted(unique_orders):
#         if order.lower() == "unknown":
#             order_to_color[order] = "lightgrey"
#         else:
#             order_to_color[order] = next(color_cycle)
#     return order_to_color


def build_order_maps_for_contigs(
    taxonomy_df: pd.DataFrame,
    contig_ids_to_plot: List[str],
    accession_col: str = "accession",
    taxonomy_col: str = "taxonomy"
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Returns:
      contig_order_map: contig_id -> order_name_or_fallback
      order_color_map:  order_name -> color  (only for orders PRESENT in contig_ids_to_plot)
    """

    # === fixed order → color mapping ===
    fixed_order_color_map: Dict[str, str] = {
        'Acholeplasmatales': '#1f77b4',
        'Bacillales': '#aec7e8',
        'Bacteria': '#ff7f0e',
        'Bdellovibrionales': '#ffbb78',
        'Candidatus Borrarchaeota': '#2ca02c',
        'Candidatus Pacearchaeota': '#98df8a',
        'Chitinophagales': '#d62728',
        'Cytophagales': '#ff9896',
        'Erysipelotrichales': '#9467bd',
        'Eubacteriales': '#c5b0d5',
        'Halobacteriales': '#8c564b',
        'Lachnospirales': '#c49c94',
        'Methanobacteriales': '#e377c2',
        'Peptostreptococcales': '#f7b6d2',
        'Rickettsiales': '#7f7f00',
        'Vampirovibrionales': '#bcbd22',
    }

    contig_order_map: Dict[str, str] = {}
    seen_orders: List[str] = []

    tax_by_acc = dict(zip(taxonomy_df[accession_col].astype(str),
                          taxonomy_df[taxonomy_col].astype(str)))

    for contig_id in contig_ids_to_plot:
        tax = tax_by_acc.get(str(contig_id), "")
        _, ord_name = extract_order_with_fallback(tax)
        contig_order_map[str(contig_id)] = ord_name
        seen_orders.append(ord_name)

    unique_orders = sorted(set(seen_orders))

    # === build final color map ===
    order_color_map: Dict[str, str] = {}
    import itertools
    color_cycle = itertools.cycle(plt.cm.tab20.colors)

    for o in unique_orders:
        if o.lower() == "unknown":
            order_color_map[o] = "#999999"
        elif o in fixed_order_color_map:
            order_color_map[o] = fixed_order_color_map[o]
        else:
            # new order not in your mapping → assign next distinct color
            order_color_map[o] = next(color_cycle)

    print("\n🔶 Order → Color mapping (shown on plot):")
    for o in unique_orders:
        print(f"  '{o}': '{order_color_map[o]}',")
    print()

    return contig_order_map, order_color_map


def plot_genomic_map(ax, features_df, prophage_start, prophage_end, allowed_genes):
    prophage_len = max(1, prophage_end - prophage_start)
    features = []
    for _, row in features_df.iterrows():
        gene = row.get("yutin", "")
        if gene in allowed_genes:
            color = "blue"
        elif pd.notna(row.get("yutin")) and str(row["yutin"]).strip() != "":
            color = "red"
        else:
            color = "lightgrey"

        features.append(
            GraphicFeature(
                start=int(row["start"]) - prophage_start,
                end=int(row["end"]) - prophage_start,
                strand=1 if str(row["strand"]).strip() == "+" else -1,
                color=color,
                label=gene if gene in allowed_genes else None
            )
        )

    record = GraphicRecord(sequence_length=prophage_len, features=features)
    record.plot(ax=ax, with_ruler=False, strand_in_label_threshold=12)
    ax.set_xlim(0, prophage_len)
    ax.set_xticks([]); ax.set_yticks([])



# def plot_genomic_maps_on_subplot(
#     ax: Axes,
#     y_positions: List[int],
#     prophage_ids: List[str],
#     feature_dict: Dict[str, List[GraphicFeature]],
#     height: float = 0.8
# ) -> None:
#     matched = 0
#     for y, prophage_id in zip(y_positions, prophage_ids):
#         features = feature_dict.get(prophage_id)
#         if not features:
#             continue
#         matched += 1
#         sequence_length = max(f.end for f in features) + 100
#         record = GraphicRecord(sequence_length=sequence_length, features=features)
#         record.plot(ax=ax, with_ruler=False, strand_in_label_threshold=7)

#     ax.set_xlim(0, 100000)
#     ax.set_xticks([])
#     ax.set_yticks([])
#     ax.set_title("Genomic map", fontsize=10)
#     ax.invert_yaxis()

#     if matched == 0:
#         print("⚠️ No genomic maps were drawn. Check if feature_dict keys match the prophage IDs.")


def build_leaf_metadata(tree: Tree) -> Tuple[Dict[str, Dict[str, str]], Dict[str, str], Dict[str, str]]:
    metadata = {}
    contig_to_protein = {}
    genome_to_protein = {}

    for leaf in tree.iter_leaves():
        protein_id = leaf.name
        if "|" in protein_id:
            contig_id = "|".join(protein_id.split("|")[:-2])
        else:
            contig_id = protein_id

        # New logic: genome_id is just contig_id now
        genome_id = contig_id

        metadata[protein_id] = {
            "protein_id": protein_id,
            "contig_id": contig_id,
            "genome_id": genome_id  # 👈 they are now equal
        }

        contig_to_protein[contig_id] = protein_id
        genome_to_protein[genome_id] = protein_id

        leaf.add_features(
            protein_id=protein_id,
            contig_id=contig_id,
            genome_id=genome_id
        )

    return metadata, contig_to_protein, genome_to_protein

def propagate_family_annotations_by_mrca(
    tree: Tree,
    leaf_to_family: Dict[str, Optional[str]],
    family_color_map: Dict[str, str]
) -> None:
    """
    Annotate unannotated leaves within pure-family clades by MRCA.
    Updates `leaf_to_family` in-place.
    """
    for family in family_color_map:
        # Find leaves already annotated with this family
        family_leaves = [leaf for leaf, fam in leaf_to_family.items() if fam == family]
        if not family_leaves:
            continue

        # Get MRCA
        mrca = tree.get_common_ancestor(family_leaves)
        mrca_leaves = mrca.get_leaves()

        # Check if all leaves are either this family or unannotated
        invalid = [
            leaf.name for leaf in mrca_leaves
            if leaf_to_family.get(leaf.name) not in {None, family}
        ]
        if invalid:
            print(f"[SKIP] Cannot annotate {family} MRCA due to mixed leaves: {invalid[:3]}...")
            continue

        # Annotate all unannotated leaves in this subtree
        for leaf in mrca_leaves:
            if leaf_to_family.get(leaf.name) is None:
                leaf_to_family[leaf.name] = family
                # print(f"[INFO] Propagated {family} to {leaf.name}")

def compare_genomes_between_prophage_and_functional_tables(
    prophage_df: pd.DataFrame,
    functional_df: pd.DataFrame
) -> None:
    """
    Compare genome names between prophage summary table and functional annotation table.
    """
    prophage_genomes = set(prophage_df["prophage_id"].unique())
    functional_genomes = set(functional_df["genome"].unique())

    # Direct matches
    common = prophage_genomes & functional_genomes
    missing = prophage_genomes - functional_genomes

    print("🧬 Prophage table:")
    print(f"  Total prophage_id entries: {len(prophage_genomes)}")
    print("🔬 Functional annotation table:")
    print(f"  Total genome entries: {len(functional_genomes)}")
    print("📊 Match summary:")
    print(f"  ✅ Matched entries: {len(common)}")
    print(f"  ❌ Missing in functional table: {len(missing)}")

    if missing:
        print("⚠️ Missing entries:")
        for m in sorted(missing):
            print(f"   - {m}")


def merge_prophage_with_taxonomy(
    prophage_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """Merge taxonomy info into prophage DataFrame using contig_id ↔ accession."""
    merged_df = prophage_df.merge(
        taxonomy_df,
        how='left',
        left_on='contig_id',
        right_on='accession'
    )
    return merged_df


def create_gridspec_axes(n_rows: int) -> Tuple[plt.Figure, Dict[str, List[plt.Axes]]]:
    fig = plt.figure(figsize=(20, n_rows * 0.45))
    # columns: barplot | family | order | taxonomy | genemap
    gs = fig.add_gridspec(n_rows, 5, width_ratios=[10, 1, 1, 12, 15], wspace=0.05)

    axes = {"barplot": [], "family": [], "order": [], "taxonomy": [], "genemap": []}
    for i in range(n_rows):
        axes["barplot"].append(fig.add_subplot(gs[i, 0]))
        axes["family"].append(fig.add_subplot(gs[i, 1], sharey=axes["barplot"][i]))
        axes["order"].append(fig.add_subplot(gs[i, 2], sharey=axes["barplot"][i]))
        axes["taxonomy"].append(fig.add_subplot(gs[i, 3], sharey=axes["barplot"][i]))
        axes["genemap"].append(fig.add_subplot(gs[i, 4], sharey=axes["barplot"][i]))
    return fig, axes


def plot_prophage_positions(
    tree: Tree,
    prophage_dict: Dict[str, Dict[str, str]],
    leaf_metadata: Dict[str, Dict[str, str]],
    leaf_to_family: Dict[str, str],
    family_to_color: Dict[str, str],
    taxonomy_df: pd.DataFrame,
    output_file: str,
    feature_dict: Optional[Dict[str, List[GraphicFeature]]] = None
) -> None:
    # === order rows by tree leaves ===
    leaf_order = [leaf.name for leaf in tree.iter_leaves()]
    ordered_prophage_dict = OrderedDict()
    taxonomy_labels = []
    contigs_in_plot: List[str] = []

    for leaf_name in leaf_order:
        genome_id = match_genome_id(leaf_name)
        for contig_id, prophage in prophage_dict.items():
            if prophage.get("prophage_id", "").startswith(genome_id):
                ordered_prophage_dict[contig_id] = prophage
                contigs_in_plot.append(contig_id)
                tax_row = taxonomy_df[taxonomy_df["accession"].astype(str) == str(contig_id)]
                taxonomy = tax_row["taxonomy"].values[0] if not tax_row.empty else "NA"
                taxonomy_labels.append(taxonomy)
                break

    n = len(ordered_prophage_dict)
    fig, axes_dict = create_gridspec_axes(n)

    # Build order colors only for what we actually plot
    contig_order_map, order_color_map = build_order_maps_for_contigs(
        taxonomy_df, contigs_in_plot
    )

    # === global x-limits so widths vary across rows ===
    contig_lengths = [int(p["contig_length"]) for p in ordered_prophage_dict.values()]
    max_contig_len = max(contig_lengths) if contig_lengths else 1
    prophage_lens = [int(p["end"]) - int(p["start"]) for p in ordered_prophage_dict.values()]
    max_prophage_len = max(prophage_lens) if prophage_lens else 1

    y_labels = []
    for i, (contig_id, prophage_info) in enumerate(ordered_prophage_dict.items()):
        prophage_id = prophage_info["prophage_id"]
        start = int(prophage_info["start"])
        end = int(prophage_info["end"])
        contig_length = int(prophage_info["contig_length"])
        prophage_len = max(1, end - start)

        ax_bar = axes_dict["barplot"][i]
        ax_fam = axes_dict["family"][i]
        ax_ord = axes_dict["order"][i]
        ax_tax = axes_dict["taxonomy"][i]
        ax_genemap = axes_dict["genemap"][i]

        # label
        matching_leaf = leaf_metadata.get(prophage_id, {}).get("protein_id")
        base_label = f"{contig_id}\n{matching_leaf}" if matching_leaf else contig_id
        full_label = f"{prophage_id} && {base_label}"
        y_labels.append(full_label)

        # 1) Contig bar (shared max scale so widths differ across rows)
        # ax_bar.broken_barh([(0, contig_length)], (0, 0.8), facecolors='lightgrey')
        # ax_bar.broken_barh([(start, end - start)], (0, 0.8), facecolors='red')
        # --- Contig bar (fill full subplot height; width is true scale) ---
        # Lock vertical axis to [0, 1] so a height=1 bar fills the subplot.
        # --- Contig bar (full subplot height in Axes coordinates) ---
        ax_bar.set_xlim(0, max_contig_len)
        ax_bar.set_ylim(0, 1)  # Axes coords in Y

        # Clear ticks
        ax_bar.set_yticks([])
        ax_bar.set_xticks([])

        # Background bar (full contig length, full subplot height)
        ax_bar.add_patch(mpatches.Rectangle(
            (0, 0),                # x, y
            contig_length,         # width in data coords
            1,                     # height in Axes coords
            transform=ax_bar.get_xaxis_transform(),  # lock height to Axes
            facecolor='lightgrey',
            edgecolor='none'
        ))

        # Prophage interval
        ax_bar.add_patch(mpatches.Rectangle(
            (start, 0),
            end - start,
            1,
            transform=ax_bar.get_xaxis_transform(),
            facecolor='red',
            edgecolor='none'
        ))

        # Label
        ax_bar.set_ylabel(full_label, fontsize=6, ha='right', rotation=0)

        # 2) Family box  (fills entire mini-axis)
        # figure out the family for this row
        family = None
        if matching_leaf:
            family = leaf_to_family.get(matching_leaf)

        if family is None:
            # fallback: match by genome prefix
            genome_prefix = prophage_id.split("|")[0]
            family = next(
                (leaf_to_family[l] for l in leaf_to_family if l.startswith(genome_prefix)),
                None
            )

        fam_color = family_to_color.get(family, "#999999") if family else "#FFFFFF"
        ax_fam.set_xlim(0, 1); ax_fam.set_ylim(0, 1)
        ax_fam.set_xticks([]); ax_fam.set_yticks([])
        ax_fam.patch.set_facecolor(fam_color)  # full fill
        ax_fam.add_patch(mpatches.Rectangle((0, 0), 1, 1,
                                            transform=ax_fam.transAxes,
                                            fill=False, edgecolor="black", linewidth=0.4))

        # 3) Order box  (fills entire mini-axis)
        order_name = contig_order_map.get(contig_id, "unknown")
        ord_color = order_color_map.get(order_name, "#999999")
        ax_ord.set_xlim(0, 1); ax_ord.set_ylim(0, 1)
        ax_ord.set_xticks([]); ax_ord.set_yticks([])
        ax_ord.patch.set_facecolor(ord_color)
        ax_ord.add_patch(mpatches.Rectangle((0, 0), 1, 1,
                                            transform=ax_ord.transAxes,
                                            fill=False, edgecolor="black", linewidth=0.4))

        # 4) Taxonomy text
        ax_tax.text(0, 0.4, taxonomy_labels[i], fontsize=7, va='center', ha='left')
        ax_tax.set_xlim(0, 1); ax_tax.set_xticks([]); ax_tax.set_yticks([])

        # 5) Genomic map (dna_features_viewer) — shared max prophage scale
        # Genomic map (prophage-only, true scale)
        if feature_dict is not None:
            # you’re using feature_dict (pre-built GraphicFeatures), which is fine;
            # but if you want to use the raw table instead, filter the DF and call plot_genomic_map().
            # Example using your pre-built dict:
            feats = feature_dict.get(prophage_id)
            if feats:
                filtered = [f for f in feats if (start <= f.start <= f.end <= end)]
                if filtered:
                    shifted = [
                        GraphicFeature(
                            start=f.start - start, end=f.end - start,
                            strand=f.strand, color=f.color, label=f.label
                        ) for f in filtered
                    ]
                    record = GraphicRecord(sequence_length=prophage_len, features=shifted)
                    record.plot(ax=ax_genemap, with_ruler=False, strand_in_label_threshold=12)
                    ax_genemap.set_xlim(0, max_prophage_len)  # or prophage_len if you prefer per-row scaling
                else:
                    ax_genemap.set_xticks([]); ax_genemap.set_yticks([])
        else:
            # Example if you want to build from the functional TSV directly:
            func_df = taxonomy_df  # (placeholder; use the actual functional-annotation DataFrame)
            func_df = func_df[func_df["genome"] == prophage_id]
            allowed_genes = {"TerL","portal","MCP","Ttub","Tstab","primase","PolB","PolA","DnaB"}
            plot_genomic_map(ax_genemap, func_df, start, end, allowed_genes)
            ax_genemap.set_xlim(0, max_prophage_len)


    # Titles
    axes_dict["barplot"][0].set_title("Prophage positions in contigs")
    axes_dict["family"][0].set_title("Crassvirales family")
    axes_dict["order"][0].set_title("Bacterial order")
    axes_dict["taxonomy"][0].set_title("Taxonomy")
    axes_dict["genemap"][0].set_title("Genomic map (prophage)")

    # Order legend (ONLY orders present on this plot)
    handles = [mpatches.Patch(color=order_color_map[o], label=o) for o in sorted(order_color_map)]
    fig.legend(handles=handles, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.01), frameon=False, fontsize=8)

    fig.tight_layout(rect=[0, 0.05, 1, 1])  # leave space for legend
    plt.savefig(output_file, dpi=300)
    plt.savefig(Path(output_file).with_suffix(".svg"))
    plt.close()


# def add_taxonomy_subplot(
#     fig: plt.Figure,
#     ax_main: plt.Axes,
#     genome_ids: List[str],
#     taxonomies: List[str],
#     subplot_width: float = 0.25
# ) -> plt.Axes:
#     """
#     Add a new subplot to the right showing taxonomy info for each genome.

#     Args:
#         fig: The main matplotlib figure.
#         ax_main: The main axis (e.g., barplot).
#         genome_ids: List of genome IDs (y-axis positions).
#         taxonomies: List of corresponding taxonomy strings.
#         subplot_width: Fractional width of the subplot (e.g. 0.25 = 25% of total).
#     Returns:
#         ax_taxonomy: New subplot axis.
#     """
#     # from matplotlib import gridspec

#     # Get position of main axis
#     bbox = ax_main.get_position()
#     total_width = bbox.width + subplot_width
#     fig.clear()

#     # Setup new grid
#     gs = fig.add_gridspec(1, 2, width_ratios=[bbox.width, subplot_width], wspace=0.05)
#     ax_main_new = fig.add_subplot(gs[0, 0])
#     ax_taxonomy = fig.add_subplot(gs[0, 1], sharey=ax_main_new)

#     # Re-plot the main plot if needed
#     # (optional: you might want to pass drawing code here again)

#     # Plot taxonomy text
#     ax_taxonomy.set_xlim(0, 1)
#     ax_taxonomy.set_xticks([])
#     ax_taxonomy.set_yticks([])

#     for i, (gid, tax) in enumerate(zip(genome_ids, taxonomies)):
#         if pd.isna(tax):
#             tax = "NA"
#         ax_taxonomy.text(0, i, tax, va='center', ha='left', fontsize=8, clip_on=True)

#     ax_taxonomy.set_title("Taxonomy", fontsize=10)
#     ax_taxonomy.invert_yaxis()  # Align with main plot
#     return ax_taxonomy



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
    prophage_dict: Dict[str, Dict[str, str]],
    leaf_metadata: Dict[str, Dict[str, str]]
) -> None:
    matched = 0

    for leaf in tree.iter_leaves():
        meta = leaf_metadata[leaf.name]
        genome_id = meta["genome_id"]

        matched_prophage = None
        for prophage in prophage_dict.values():
            if prophage.get("prophage_id", "").startswith(genome_id):
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
            # print(f"[MATCHED] {leaf.name} ← {matched_prophage['prophage_id']}")

    print(f"[INFO] Matched {matched} prophage leaves out of {len(tree.get_leaves())}")

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
        # from PyQt5.QtGui import QPainter, QColor
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


def print_unmatched_prophages(prophage_dict, tree: Tree):
    tree_leaf_prefixes = {match_genome_id(leaf.name) for leaf in tree.iter_leaves()}

    unmatched = []
    for prophage in prophage_dict.values():
        prophage_id = prophage.get("prophage_id", "")
        prefix = prophage_id.split("|")[0]
        if prefix not in tree_leaf_prefixes:
            unmatched.append(prophage_id)

    print("\n❌ Prophages not matched to tree leaves:")
    for pid in unmatched:
        print(f"   - {pid}")



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
    crassus_output = "/mnt/c/crassvirales/CrassUS_old/CrassUS/results/prophage_analysis"
    crassus_phylogeny = f"{crassus_output}/5_phylogenies/3_iToL"
    tree_file = Path(f"{crassus_phylogeny}/TerL_iToL_renamed.nwk")
    itol_annotation = Path(f"{crassus_phylogeny}/TerL_iToL.txt")
    prophage_table = Path(f"{wd}/prophage_alignments_ncbi_and_gtdb_samtools_coordinates_edited_with_lengths.tsv")
    taxonomy_file = Path(f"{wd}/prophage_alignments_ncbi_and_gtdb_samtools_coordinates_edited_with_lengths_bacterial_contig_ids_taxonomy.tsv")
    functional_annotation_file = Path(f"{crassus_output}/4_ORF/3_functional_annot_table_renamed.tsv")

    output_svg = Path(f"{wd}/TerL_with_prophages_annotated_tree.svg")
    prophage_plot = Path(f"{wd}/prophage_positions_plot.png")


    tree = read_tree(tree_file)
    leaf_to_family, family_to_color, _ = parse_itol_annotation(itol_annotation)

    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", dtype=str)

    feature_dict = load_functional_annotation(str(functional_annotation_file))

    # print(f"{feature_dict=}")

    # NEW: propagate family annotation by MRCA
    propagate_family_annotations_by_mrca(tree, leaf_to_family, family_to_color)

    # === Build consistent metadata mapping for each leaf ===
    #leaf_metadata = build_leaf_metadata(tree)
    leaf_metadata, contig_to_protein, genome_to_protein = build_leaf_metadata(tree)

    # print(f'{leaf_metadata=}\n')

    # print(f'{contig_to_protein=}\n')

    # print(f'{genome_to_protein=}\n')

    prophage_dict, prophage_df = parse_prophage_table(prophage_table)

    merged_df = merge_prophage_with_taxonomy(prophage_df, taxonomy_df)


    # leaf_order = [leaf.name for leaf in tree.iter_leaves()]
    # ordered_prophage_dict = OrderedDict()
    # for leaf_name in leaf_order:
    #     genome_id = match_genome_id(leaf_name)
    #     for contig_id, prophage in prophage_dict.items():
    #         prophage_id = prophage.get("prophage_id", "")
    #         if prophage_id.startswith(genome_id):
    #             ordered_prophage_dict[contig_id] = prophage


    # prophage_dict = ordered_prophage_dict
    # print(prophage_dict)

    leaf_ids = {match_genome_id(leaf.name) for leaf in tree.iter_leaves()}
    # matched_ids = set(prophage_dict.keys()) & leaf_ids
    # print(f"[INFO] Matching prophage IDs: {len(matched_ids)} out of {len(prophage_dict)} total prophages")

    matched = sum(
        any(leaf.name.startswith(prophage.get("prophage_id", "")) for leaf in tree.iter_leaves())
        for prophage in prophage_dict.values()
    )
    print(f"[INFO] Matching prophage IDs: {matched} out of {len(prophage_dict)} total prophages")


    # annotate_tree_features(tree, leaf_to_family, family_to_color, prophage_dict)
    annotate_tree_features(tree, leaf_to_family, family_to_color, prophage_dict, leaf_metadata)
    render_tree(tree, output_svg, family_to_color)
    print(f"Annotated tree saved to: {output_svg}")

    print("🔍 Example keys in feature_dict:")
    print(list(feature_dict.keys())[:5])

    #print("🔍 First 5 prophage_ids:")
    #print(prophage_dict.keys()[:5])

    # New: plot genomic positions
    # plot_prophages(prophage_df, output_file=str(prophage_plot))
    # plot_prophage_positions(prophage_dict, output_file=str(prophage_plot))

    compare_genomes_between_prophage_and_functional_tables(
        prophage_df=prophage_df,  # ✅ Already read as DataFrame
        functional_df=pd.read_csv(functional_annotation_file, sep="\t")  # ✅ Explicit load
    )

    print_unmatched_prophages(prophage_dict, tree)

    plot_prophage_positions(
        tree,
        prophage_dict,
        leaf_metadata,
        leaf_to_family,
        family_to_color,
        merged_df,
        output_file=str(prophage_plot),
        feature_dict=feature_dict
    )

    print(f"Prophage plot saved to: {prophage_plot}")


if __name__ == "__main__":
    main()
