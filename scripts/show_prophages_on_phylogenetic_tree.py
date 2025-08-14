#!/usr/bin/env python3

import os
import re
from collections import OrderedDict
from pathlib import Path

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
from ete3 import Face, NCBITaxa, RectFace, TextFace, Tree, TreeStyle
from matplotlib import gridspec
from PyQt5.QtGui import QColor, QPainter

matplotlib.use("Agg")
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Ensure Qt offscreen rendering

DEFAULT_COLOR = "#bfbfbf"  # neutral grey

TITLE_FONTSIZE = 20  # or 11 if you want them same as y-axis labels

TEXT_FONTSIZE = 12

LEGEND_FONTSIZE = 20  # bigger text
LEGEND_TITLE_FONTSIZE = 30  # bigger text
LEGEND_HANDLE_HEIGHT = 5.0  # taller color boxes
LEGEND_HANDLE_LENGTH = 5.0  # wider color boxes

ORDER_LIKE = re.compile(r".*ales$", re.IGNORECASE)

# top-level (after imports)
try:
    _NCBI = NCBITaxa()
except Exception:
    _NCBI = None


def read_tree(tree_file: Path) -> Tree:
    return Tree(str(tree_file), format=1)


def adjust_left_margin_for_labels(fig, ax_bar_list):
    renderer = fig.canvas.get_renderer()
    max_label_width = 0
    for ax in ax_bar_list:
        labels = [lab.get_text() for lab in ax.get_yticklabels()] + [ax.get_ylabel()]
        for label in labels:
            if label:
                t = ax.text(0, 0, label)
                bbox = t.get_window_extent(renderer=renderer)
                max_label_width = max(max_label_width, bbox.width)
                t.remove()
    # Convert width in pixels to figure fraction
    fig_width_inch = fig.get_size_inches()[0]
    dpi = fig.dpi
    margin_frac = max_label_width / (fig_width_inch * dpi)
    # Adjust left margin
    fig.subplots_adjust(left=margin_frac + 0.02)  # small padding


def load_functional_annotation(path: str) -> dict[str, list[GraphicFeature]]:
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    features_by_genome: dict[str, list[GraphicFeature]] = {}

    allowed = {"TerL", "portal", "MCP", "Ttub", "Tstab", "primase", "PolB", "PolA", "DnaB"}

    for _, row in df.iterrows():
        genome = str(row["genome"]).strip()  # <- strip to avoid key mismatches
        s = int(row["start"])
        e = int(row["end"])
        # normalize
        start, end = (s, e) if s <= e else (e, s)

        strand = 1 if str(row["strand"]).strip() == "+" else -1
        yutin_raw = row.get("yutin", None)
        yutin = str(yutin_raw).strip() if pd.notna(yutin_raw) else ""

        color = "red" if yutin else "lightgrey"
        label = yutin if yutin in allowed else None

        features_by_genome.setdefault(genome, []).append(
            GraphicFeature(start=start, end=end, strand=strand, color=color, label=label)
        )
    return features_by_genome


def _split_taxonomy_string(tax: str) -> list[str]:
    if not isinstance(tax, str) or not tax.strip():
        return []
    toks = [t.strip() for t in tax.split(";") if t.strip()]
    # drop “environmental samples”
    return [t for t in toks if t.lower() != "environmental samples"]


def _try_ncbi_order(tokens: list[str]) -> str | None:
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


def _try_ncbi_rank(tokens: list[str], rank: str) -> str | None:
    """Return the first token that is that exact NCBI rank, if possible."""
    if _NCBI is None:
        return None
    try:
        name2tax = _NCBI.get_name_translator(tokens)
        if not name2tax:
            return None
        taxids = [name2tax[t][0] for t in tokens if t in name2tax]
        ranks = _NCBI.get_rank(taxids)  # {taxid: rank}
        for tok in tokens:
            tids = name2tax.get(tok)
            if not tids:
                continue
            if ranks.get(tids[0]) == rank:
                return tok
    except Exception:
        return None
    return None


def _heuristic_order(tokens: list[str]) -> str | None:
    # pick first token that looks like an order (…ales)
    for tok in tokens:
        if ORDER_LIKE.match(tok):
            return tok
    return None


def _pick_suffix(tokens: list[str], suffixes: tuple[str, ...]) -> str | None:
    for tok in tokens:
        if tok.endswith(suffixes):
            return tok
    return None


def extract_order_with_fallback(taxonomy: str) -> tuple[str, str]:
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
    cls = _pick_suffix(
        tokens, ("ia", "bia", "illi", "idia", "eria")
    )  # includes Bacilli       # Clostridia, Halobacteria
    if cls:
        return "class", cls
    phyl = _pick_suffix(tokens, ("ota", "otae"))  # Bacillota, Bacteroidota
    if phyl:
        return "phylum", phyl
    fam = _pick_suffix(tokens, ("aceae",))  # Lachnospiraceae
    if fam:
        return "family", fam

    return "unknown", "Unknown"


def extract_class_with_fallback(taxonomy: str) -> tuple[str, str]:
    """
    Returns (used_rank, name) where used_rank ∈ {'class','phylum','order','family','unknown'} for bacterial class display.
    Priority: NCBITaxa(class) → heuristic class by suffix → order → phylum → family → Unknown.
    """
    tokens = _split_taxonomy_string(taxonomy)
    if not tokens:
        return "unknown", "Unknown"

    # 1) NCBI
    cls = _try_ncbi_rank(tokens, "class")
    if cls:
        return "class", cls

    # 2) Heuristic for class-like suffixes
    cls = _pick_suffix(tokens, ("ia", "bia", "illi", "idia", "eria"))  # e.g., Bacilli, Clostridia, Halobacteria
    if cls:
        return "class", cls

    # 3) Loose fallbacks so we still color something if present
    ord_name = _pick_suffix(tokens, ("ales",))
    if ord_name:
        return "order", ord_name
    phyl = _pick_suffix(tokens, ("ota", "otae"))
    if phyl:
        return "phylum", phyl
    fam = _pick_suffix(tokens, ("aceae",))
    if fam:
        return "family", fam

    return "unknown", "Unknown"


def extract_phylum_with_fallback(taxonomy: str) -> tuple[str, str]:
    """
    Returns (used_rank, name) for phylum:
    Priority: NCBITaxa('phylum') → heuristic suffix (…ota/…otae/…ote) → 'Unknown'
    """
    tokens = _split_taxonomy_string(taxonomy)
    if not tokens:
        return "unknown", "Unknown"

    # 1) Exact NCBI rank if available
    phy = _try_ncbi_rank(tokens, "phylum")
    if phy:
        return "phylum", phy

    # 2) Heuristic by common phylum suffixes
    phy = _pick_suffix(tokens, ("ota", "otae", "ote"))
    if phy:
        return "phylum", phy

    return "unknown", "Unknown"


def extract_order_strict(taxonomy: str) -> str:
    """
    Return only the ORDER name if present, else 'Unknown'.
    Priority: exact NCBI rank 'order' → heuristic '...ales' → Unknown.
    """
    tokens = _split_taxonomy_string(taxonomy)
    if not tokens:
        return "Unknown"

    # exact NCBI 'order'
    order = _try_ncbi_rank(tokens, "order")
    if order:
        return order

    # heuristic suffix
    for tok in tokens:
        if ORDER_LIKE.match(tok):
            return tok

    return "Unknown"


def extract_family_strict(taxonomy: str) -> str:
    """
    Return only the FAMILY name if present, else 'Unknown'.
    Priority: exact NCBI rank 'family' → heuristic '...aceae' → Unknown.
    """
    tokens = _split_taxonomy_string(taxonomy)
    if not tokens:
        return "Unknown"

    fam = _try_ncbi_rank(tokens, "family")
    if fam:
        return fam

    for tok in tokens:
        if tok.endswith("aceae"):
            return tok

    return "Unknown"


def format_order_family_label(taxonomy: str, multiline: bool = True) -> str:
    """
    Build the text shown in the 'Taxonomy lineage' column: only Order and Family.
    If both missing → 'Unknown'. If multiline=True, put them on separate lines.
    """
    order = extract_order_strict(taxonomy)
    family = extract_family_strict(taxonomy)

    parts = []
    if order != "Unknown":
        parts.append(order)
    if family != "Unknown":
        parts.append(family)

    if not parts:
        return "Unknown"

    return ";".join(parts)


def format_order_family_or_full(taxonomy: str, multiline: bool = True) -> str:
    """
    Show only Order and Family; if neither can be extracted, return the original taxonomy string.
    """
    # Try extracting
    order = extract_order_strict(taxonomy)
    family = extract_family_strict(taxonomy)

    parts = []
    if order != "Unknown":
        parts.append(order)
    if family != "Unknown":
        parts.append(family)

    if parts:  # we found at least one of order/family
        return ("\n" if multiline else "; ").join(parts)

    # fallback: original taxonomy (could be empty if missing in table)
    return taxonomy if isinstance(taxonomy, str) and taxonomy.strip() else "Unknown"


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


# def build_order_maps_for_contigs(
#     taxonomy_df: pd.DataFrame,
#     contig_ids_to_plot: List[str],
#     accession_col: str = "accession",
#     taxonomy_col: str = "taxonomy"
# ) -> Tuple[Dict[str, str], Dict[str, str]]:
#     """
#     Returns:
#       contig_order_map: contig_id -> order_name_or_fallback
#       order_color_map:  order_name -> color  (only for orders PRESENT in contig_ids_to_plot)
#     """

#     # === fixed order → color mapping ===
#     fixed_order_color_map = {
#         # Mollicutes
#         'Acholeplasmatales':     '#9467bd',  # Mollicutes
#         # Bacilli
#         'Bacillales':            '#1f77b4',  # Bacilli
#         # use 'Unknown' instead of 'Bacteria' as an “order”
#         'Unknown':               '#999999',

#         # Bdellovibrionia
#         'Bdellovibrionales':     '#aec7e8',  # Bdellovibrionia

#         # Candidatus groups (keep same as class)
#         'Candidatus Borrarchaeota': '#ff7f0e',
#         'Candidatus Pacearchaeota': '#fbc15e',

#         # Bacteroidetes/CFB group
#         'Chitinophagales':       '#2ca02c',  # Chitinophagia
#         'Cytophagales':          '#98df8a',  # Cytophagia

#         # Firmicutes (Clostridia)
#         'Erysipelotrichales':    '#e377c2',  # Erysipelotrichia
#         'Eubacteriales':         '#ff7f0e',  # Clostridia (deprecated order but maps here)
#         'Lachnospirales':        '#ff7f0e',  # Clostridia
#         'Peptostreptococcales':  '#ff7f0e',  # Clostridia

#         # Archaea
#         'Halobacteriales':       '#7f7f7f',  # Halobacteria (see note above)
#         'Methanobacteriales':    '#ff9896',  # Methanobacteria

#         # Proteobacteria (Alpha)
#         'Rickettsiales':         '#d62728',  # Alphaproteobacteria

#         # Candidate phylum/class Vampirovibriophyceae
#         'Vampirovibrionales':    '#c5b0d5',  # Vampirovibriophyceae
#     }

#     contig_order_map: Dict[str, str] = {}
#     seen_orders: List[str] = []

#     tax_by_acc = dict(zip(taxonomy_df[accession_col].astype(str),
#                           taxonomy_df[taxonomy_col].astype(str)))

#     for contig_id in contig_ids_to_plot:
#         tax = tax_by_acc.get(str(contig_id), "")
#         _, ord_name = extract_order_with_fallback(tax)
#         contig_order_map[str(contig_id)] = ord_name
#         seen_orders.append(ord_name)

#     unique_orders = sorted(set(seen_orders))

#     # === build final color map ===
#     order_color_map: Dict[str, str] = {}
#     import itertools
#     color_cycle = itertools.cycle(plt.cm.tab20.colors)

#     for o in unique_orders:
#         if o.lower() == "unknown":
#             order_color_map[o] = "#999999"
#         elif o in fixed_order_color_map:
#             order_color_map[o] = fixed_order_color_map[o]
#         else:
#             # new order not in your mapping → assign next distinct color
#             order_color_map[o] = next(color_cycle)

#     print("\n🔶 Order → Color mapping (shown on plot):")
#     for o in unique_orders:
#         print(f"  '{o}': '{order_color_map[o]}',")
#     print()

#     return contig_order_map, order_color_map


def build_class_maps_for_contigs(
    taxonomy_df: pd.DataFrame,
    contig_ids_to_plot: list[str],
    accession_col: str = "accession",
    taxonomy_col: str = "taxonomy",
) -> tuple[dict[str, str], dict[str, str]]:
    """
    Returns:
      contig_class_map: contig_id -> class_name_or_fallback
      class_color_map:  class_name -> color (only for classes present)
    """

    # Fixed mapping (pick distinct tab20-ish colors; edit if you want specific classes)
    class_color_map = {
        # Bacterial classes
        # Bacillota
        "Bacilli": "#e6194B",
        "Clostridia": "#f58231",
        "Erysipelotrichia": "#911eb4",
        # Bacteroidota
        "Chitinophagia": "#f032e6",
        "Cytophagia": "#000075",
        "Ignavibacteria": "#fabed4",
        "Bdellovibrionia": "#bfef45",
        "Mollicutes": "#911eb4",
        "Vampirovibriophyceae": "#42d4f4",
        "Unknown": "#999999",  # fallback for unassigned bacterial classes
        # Archaeal classes (black)
        "Candidatus Borrarchaeota": "#000000",
        "Candidatus Pacearchaeota": "#000000",
        "Halobacteria": "#000000",
        "Methanobacteria": "#000000",
    }

    contig_class_map: dict[str, str] = {}
    seen_classes: list[str] = []

    tax_by_acc = dict(zip(taxonomy_df[accession_col].astype(str), taxonomy_df[taxonomy_col].astype(str), strict=False))

    for contig_id in contig_ids_to_plot:
        tax = tax_by_acc.get(str(contig_id), "")
        _, cls_name = extract_class_with_fallback(tax)
        contig_class_map[str(contig_id)] = cls_name
        seen_classes.append(cls_name)

    unique_classes = sorted(set(seen_classes))

    # Final color map
    # class_color_map: Dict[str, str] = {}
    # import itertools
    # color_cycle = itertools.cycle(plt.cm.tab20.colors)

    # for c in unique_classes:
    #     if c.lower() == "unknown":
    #         class_color_map[c] = "#999999"
    #     elif c in class_color_map:
    #         class_color_map[c] = class_color_map[c]
    #     else:
    #         class_color_map[c] = next(color_cycle)

    # (optional) print what’s used
    # print("\n🔷 Class → Color mapping (shown on plot):")
    # for c in unique_classes:
    #     print(f"  '{c}': '{class_color_map[c]}',")
    # print()

    return contig_class_map, class_color_map


def build_phylum_maps_for_contigs(
    taxonomy_df: pd.DataFrame,
    contig_ids_to_plot: list[str],
    accession_col: str = "accession",
    taxonomy_col: str = "taxonomy",
    color_map: dict[str, str] = None,
) -> tuple[dict[str, str], dict[str, str]]:
    """
    Returns:
      contig_phylum_map: contig_id -> phylum_name_or_fallback
      used_phylum_color_map: only colors for phyla present in the plot (no new colors generated)
    """

    # ===== Phylum colors (customizable) =====
    phylum_color_map = {
        # Bacteria phyla
        "Bacillota": "#1f77b4",  # from Bacilli
        "Firmicutes": "#1f77b4",
        "Bdellovibrionota": "#bfef45",  # from Bdellovibrionia
        "Bacteroidota": "#2ca02c",  # from Chitinophagia
        "Bacteroidetes": "#2ca02c",
        "Proteobacteria": "#9A6324",  # covers Alpha/others
        "Pseudomonadota": "#9A6324",
        "Mycoplasmatota": "#911eb4",
        "Chloroflexota": "#ffe119",
        "Ignavibacteriota": "#fabed4",
        "Candidatus Melainabacteria": "#42d4f4",
        # Unused colors
        # 'Thermotogota':           '#ff7f0e',
        # 'Cyanobacteriota':        '#98df8a',  # paired with Cytophagia tone
        # 'Actinomycetota':         '#e377c2',  # paired with Erysipelotrichia tone
        # 'Planctomycetota':        '#c5b0d5',  # paired with Vampirovibriophyceae tone
        # Archaea (always black)
        "Candidatus Borrarchaeota": "#000000",
        "Candidatus Pacearchaeota": "#000000",
        # 'Euryarchaeota':            '#000000',
        # 'Halobacteriota':           '#000000',
        "Methanobacteriota": "#000000",
        # 'Nanoarchaeota':            '#000000',
        # 'Thermoproteota':           '#000000',
        # Fallback
        "Unknown": DEFAULT_COLOR,
    }

    if color_map is None:
        color_map = phylum_color_map  # use the customizable global dict

    contig_phylum_map: dict[str, str] = {}
    seen_phyla: list[str] = []

    tax_by_acc = dict(zip(taxonomy_df[accession_col].astype(str), taxonomy_df[taxonomy_col].astype(str), strict=False))

    for contig_id in contig_ids_to_plot:
        tax = tax_by_acc.get(str(contig_id), "")
        _, phy_name = extract_phylum_with_fallback(tax)
        contig_phylum_map[str(contig_id)] = phy_name
        seen_phyla.append(phy_name)

    unique_phyla = sorted(set(seen_phyla))

    # Only use provided dict; unknowns get DEFAULT_COLOR
    used_phylum_color_map: dict[str, str] = {}
    for p in unique_phyla:
        used_phylum_color_map[p] = color_map.get(p, DEFAULT_COLOR)

    # Optional: print what’s used
    print("\n🔶 Phylum → Color mapping (shown on plot):")
    for p in unique_phyla:
        print(f"  '{p}': '{used_phylum_color_map[p]}',")
    print()

    return contig_phylum_map, used_phylum_color_map


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
                label=gene if gene in allowed_genes else None,
            )
        )

    record = GraphicRecord(sequence_length=prophage_len, features=features)
    record.plot(ax=ax, with_ruler=False, strand_in_label_threshold=12)
    ax.set_xlim(0, prophage_len)
    ax.set_xticks([])
    ax.set_yticks([])


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


def build_leaf_metadata(tree: Tree) -> tuple[dict[str, dict[str, str]], dict[str, str], dict[str, str]]:
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
            "genome_id": genome_id,  # 👈 they are now equal
        }

        contig_to_protein[contig_id] = protein_id
        genome_to_protein[genome_id] = protein_id

        leaf.add_features(protein_id=protein_id, contig_id=contig_id, genome_id=genome_id)

    return metadata, contig_to_protein, genome_to_protein


def propagate_family_annotations_by_mrca(
    tree: Tree, leaf_to_family: dict[str, str | None], family_color_map: dict[str, str]
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
        invalid = [leaf.name for leaf in mrca_leaves if leaf_to_family.get(leaf.name) not in {None, family}]
        if invalid:
            print(f"[SKIP] Cannot annotate {family} MRCA due to mixed leaves: {invalid[:3]}...")
            continue

        # Annotate all unannotated leaves in this subtree
        for leaf in mrca_leaves:
            if leaf_to_family.get(leaf.name) is None:
                leaf_to_family[leaf.name] = family
                # print(f"[INFO] Propagated {family} to {leaf.name}")


def compare_genomes_between_prophage_and_functional_tables(
    prophage_df: pd.DataFrame, functional_df: pd.DataFrame
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


def merge_prophage_with_taxonomy(prophage_df: pd.DataFrame, taxonomy_df: pd.DataFrame) -> pd.DataFrame:
    """Merge taxonomy info into prophage DataFrame using contig_id ↔ accession."""
    merged_df = prophage_df.merge(taxonomy_df, how="left", left_on="contig_id", right_on="accession")
    return merged_df


def create_gridspec_axes(n_rows: int) -> tuple[plt.Figure, dict[str, list[plt.Axes]]]:
    # Increase figure height to accommodate legends
    fig = plt.figure(figsize=(60, n_rows * 0.7 + 2))

    # columns: barplot | family | origin | phylum | class | taxonomy | genemap
    gs = gridspec.GridSpec(
        n_rows + 1,
        7,  # +1 row for legends; columns are 0..6
        width_ratios=[20, 1.2, 1.2, 1.2, 1.2, 20, 90],
        wspace=0.05,
        height_ratios=[1] * n_rows + [0.6],
        hspace=0.2,
    )

    axes = {
        "barplot": [],
        "family": [],
        "is_mag": [],
        "phylum": [],
        "class": [],
        "taxonomy": [],
        "genemap": [],
        "legend": None,
    }

    for i in range(n_rows):
        axes["barplot"].append(fig.add_subplot(gs[i, 0]))
        axes["family"].append(fig.add_subplot(gs[i, 1], sharey=axes["barplot"][i]))
        axes["is_mag"].append(fig.add_subplot(gs[i, 2], sharey=axes["barplot"][i]))
        axes["phylum"].append(fig.add_subplot(gs[i, 3], sharey=axes["barplot"][i]))
        axes["class"].append(fig.add_subplot(gs[i, 4], sharey=axes["barplot"][i]))
        axes["taxonomy"].append(fig.add_subplot(gs[i, 5], sharey=axes["barplot"][i]))  # was 6
        axes["genemap"].append(fig.add_subplot(gs[i, 6], sharey=axes["barplot"][i]))  # was 7

    axes["legend"] = fig.add_subplot(gs[n_rows, :])
    axes["legend"].axis("off")

    return fig, axes


def plot_prophage_positions(
    tree: Tree,
    prophage_dict: dict[str, dict[str, str]],
    leaf_metadata: dict[str, dict[str, str]],
    leaf_to_family: dict[str, str],
    family_to_color: dict[str, str],
    taxonomy_df: pd.DataFrame,
    output_file: str,
    feature_dict: dict[str, list[GraphicFeature]] | None = None,
) -> None:
    # === order rows by tree leaves ===
    leaf_order = [leaf.name for leaf in tree.iter_leaves()]
    ordered_prophage_dict = OrderedDict()
    taxonomy_labels = []
    contigs_in_plot: list[str] = []

    for leaf_name in leaf_order:
        genome_id = match_genome_id(leaf_name)
        for contig_id, prophage in prophage_dict.items():
            if prophage.get("prophage_id", "").startswith(genome_id):
                ordered_prophage_dict[contig_id] = prophage
                contigs_in_plot.append(contig_id)
                # tax_row = taxonomy_df[taxonomy_df["accession"].astype(str) == str(contig_id)]
                # taxonomy = tax_row["taxonomy"].values[0] if not tax_row.empty else "NA"
                # taxonomy_labels.append(taxonomy)
                tax_row = taxonomy_df[taxonomy_df["accession"].astype(str) == str(contig_id)]
                taxonomy = tax_row["taxonomy"].values[0] if not tax_row.empty else ""
                # taxonomy_labels.append(format_order_family_label(taxonomy))
                taxonomy_labels.append(format_order_family_or_full(taxonomy, multiline=False))
                break

    n = len(ordered_prophage_dict)
    fig, axes_dict = create_gridspec_axes(n)

    # Color map for origin (Isolate vs MAG)
    origin_color_map = {
        "Isolate": "#d62728",  # red
        "MAG": "#1f77b4",  # blue
    }

    contig_phylum_map, used_phylum_color_map = build_phylum_maps_for_contigs(taxonomy_df, contigs_in_plot)

    # NEW: build class colors only for what we actually plot
    contig_class_map, class_color_map = build_class_maps_for_contigs(taxonomy_df, contigs_in_plot)

    # Build order colors only for what we actually plot
    # contig_order_map, order_color_map = build_order_maps_for_contigs(
    #     taxonomy_df, contigs_in_plot
    # )

    # === global x-limits so widths vary across rows ===
    contig_lengths = [int(p["contig_length"]) for p in ordered_prophage_dict.values()]
    max_contig_len = max(contig_lengths) if contig_lengths else 1
    prophage_lens = [int(p["end"]) - int(p["start"]) for p in ordered_prophage_dict.values()]
    max_prophage_len = max(prophage_lens) if prophage_lens else 1

    y_labels = []

    families_used = set()

    for i, (contig_id, prophage_info) in enumerate(ordered_prophage_dict.items()):
        prophage_id = prophage_info["prophage_id"]
        start = int(prophage_info["start"])
        end = int(prophage_info["end"])
        contig_length = int(prophage_info["contig_length"])
        prophage_len = max(1, end - start)

        ax_bar = axes_dict["barplot"][i]
        ax_fam = axes_dict["family"][i]
        # ax_ord = axes_dict["order"][i]
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
        ax_bar.add_patch(
            mpatches.Rectangle(
                (0, 0),  # x, y
                contig_length,  # width in data coords
                1,  # height in Axes coords
                transform=ax_bar.get_xaxis_transform(),  # lock height to Axes
                facecolor="lightgrey",
                edgecolor="none",
            )
        )

        # Prophage interval
        ax_bar.add_patch(
            mpatches.Rectangle(
                (start, 0), end - start, 1, transform=ax_bar.get_xaxis_transform(), facecolor="red", edgecolor="none"
            )
        )

        # Label
        ax_bar.set_ylabel(full_label, fontsize=TITLE_FONTSIZE, ha="right", rotation=0)

        ax_cls = axes_dict["class"][i]  # NEW

        # 2) Family box  (fills entire mini-axis)
        # figure out the family for this row
        family = None
        if matching_leaf:
            family = leaf_to_family.get(matching_leaf)

        if family is None:
            # fallback: match by genome prefix
            genome_prefix = prophage_id.split("|")[0]
            family = next((leaf_to_family[l] for l in leaf_to_family if l.startswith(genome_prefix)), None)

        fam_color = family_to_color.get(family, "#999999") if family else "#FFFFFF"

        if family:
            families_used.add(family)  # NEW

        ax_fam.set_xlim(0, 1)
        ax_fam.set_ylim(0, 1)
        ax_fam.set_xticks([])
        ax_fam.set_yticks([])
        ax_fam.patch.set_facecolor(fam_color)  # full fill
        ax_fam.add_patch(
            mpatches.Rectangle((0, 0), 1, 1, transform=ax_fam.transAxes, fill=False, edgecolor="black", linewidth=0.4)
        )

        # 2.5) Origin box (Isolate/MAG)
        ax_mag = axes_dict["is_mag"][i]  # NEW origin axis

        origin_raw = prophage_info.get("is_mag", "")  # column name as you wrote it
        origin = str(origin_raw).strip()
        origin_color = origin_color_map.get(origin, DEFAULT_COLOR)

        ax_mag.set_xlim(0, 1)
        ax_mag.set_ylim(0, 1)
        ax_mag.set_xticks([])
        ax_mag.set_yticks([])
        ax_mag.patch.set_facecolor(origin_color)
        ax_mag.add_patch(
            mpatches.Rectangle((0, 0), 1, 1, transform=ax_mag.transAxes, fill=False, edgecolor="black", linewidth=0.4)
        )

        # 2.7 Phylum box

        ax_phy = axes_dict["phylum"][i]  # NEW
        phylum_name = contig_phylum_map.get(contig_id, "Unknown")
        phy_color = used_phylum_color_map.get(phylum_name, DEFAULT_COLOR)

        ax_phy.set_xlim(0, 1)
        ax_phy.set_ylim(0, 1)
        ax_phy.set_xticks([])
        ax_phy.set_yticks([])
        ax_phy.patch.set_facecolor(phy_color)
        ax_phy.add_patch(
            mpatches.Rectangle((0, 0), 1, 1, transform=ax_phy.transAxes, fill=False, edgecolor="black", linewidth=0.4)
        )

        # 3) Class box (NEW)
        class_name = contig_class_map.get(contig_id, "unknown")
        cls_color = class_color_map.get(class_name, DEFAULT_COLOR)
        ax_cls.set_xlim(0, 1)
        ax_cls.set_ylim(0, 1)
        ax_cls.set_xticks([])
        ax_cls.set_yticks([])
        ax_cls.patch.set_facecolor(cls_color)
        ax_cls.add_patch(
            mpatches.Rectangle((0, 0), 1, 1, transform=ax_cls.transAxes, fill=False, edgecolor="black", linewidth=0.4)
        )

        # # 4) Order box  (fills entire mini-axis)
        # order_name = contig_order_map.get(contig_id, "unknown")
        # ord_color = order_color_map.get(order_name, DEFAULT_COLOR)
        # ax_ord.set_xlim(0, 1); ax_ord.set_ylim(0, 1)
        # ax_ord.set_xticks([]); ax_ord.set_yticks([])
        # ax_ord.patch.set_facecolor(ord_color)
        # ax_ord.add_patch(mpatches.Rectangle((0, 0), 1, 1,
        #                                     transform=ax_ord.transAxes,
        #                                     fill=False, edgecolor="black", linewidth=0.4))

        # 5) Taxonomy text
        ax_tax.text(0, 0.4, taxonomy_labels[i], fontsize=20, va="center", ha="left")
        ax_tax.set_xlim(0, 1)
        ax_tax.set_xticks([])
        ax_tax.set_yticks([])

        # 6) Genomic map (handle absolute vs local coords + debug + fallbacks)
        DEBUG_GENEMAP = False  # set True to print diagnostics

        if feature_dict is not None:
            # --- try several keys ---
            lookups = [
                prophage_id.strip(),
                contig_id.strip(),
                contig_id.split("|")[0].strip(),
            ]
            feats = None
            used_key = None
            for k in lookups:
                if k in feature_dict:
                    feats = feature_dict[k]
                    used_key = k
                    break

            if feats:
                # normalize coordinates to integers and compute span
                fs_list, fe_list = [], []
                norm = []
                for f in feats:
                    fs = int(min(int(f.start), int(f.end)))
                    fe = int(max(int(f.start), int(f.end)))
                    fs_list.append(fs)
                    fe_list.append(fe)
                    norm.append((fs, fe, f))

                fs_min = min(fs_list)
                fe_max = max(fe_list)

                # Heuristic: treat as local if features look like 1..prophage_len
                is_local = (fs_min <= 5) and (fe_max <= prophage_len * 1.2)

                if is_local:
                    if DEBUG_GENEMAP:
                        print(
                            f"[GMAP] Using LOCAL coords for {prophage_id!r} (key={used_key!r}); "
                            f"window=({start},{end}), feat span=({fs_min},{fe_max}), len={prophage_len}"
                        )
                    local_feats = []
                    for fs, fe, f in norm:
                        # clip to [0, prophage_len] and (optionally) convert 1-based → 0-based
                        cs = max(0, fs - 1)  # if your local table is 1-based
                        ce = min(prophage_len, fe)
                        if ce > cs:
                            local_feats.append(
                                GraphicFeature(start=cs, end=ce, strand=f.strand, color=f.color, label=f.label)
                            )
                    if local_feats:
                        record = GraphicRecord(sequence_length=prophage_len, features=local_feats)
                        record.plot(ax=ax_genemap, with_ruler=False, strand_in_label_threshold=12)
                        ax_genemap.set_xlim(0, max_prophage_len)
                    else:
                        ax_genemap.set_xticks([])
                        ax_genemap.set_yticks([])

                else:
                    if DEBUG_GENEMAP and not (fe_max > start and fs_min < end):
                        print(
                            f"[GMAP] ABS coords but no overlap for {prophage_id!r} (key={used_key!r}); "
                            f"window=({start},{end}), feat span=({fs_min},{fe_max}), n={len(feats)}"
                        )

                    # absolute (contig) → clip to window and shift
                    clipped = []
                    for fs, fe, f in norm:
                        if fe > start and fs < end:
                            cs = max(fs, start)
                            ce = min(fe, end)
                            if ce > cs:
                                clipped.append(
                                    GraphicFeature(
                                        start=cs - start, end=ce - start, strand=f.strand, color=f.color, label=f.label
                                    )
                                )
                    if clipped:
                        record = GraphicRecord(sequence_length=prophage_len, features=clipped)
                        record.plot(ax=ax_genemap, with_ruler=False, strand_in_label_threshold=12)
                        ax_genemap.set_xlim(0, max_prophage_len)
                    else:
                        ax_genemap.set_xticks([])
                        ax_genemap.set_yticks([])

            else:
                if DEBUG_GENEMAP:
                    print(f"[GMAP] No features for keys={lookups!r}")
                ax_genemap.set_xticks([])
                ax_genemap.set_yticks([])

        else:
            ax_genemap.set_xticks([])
            ax_genemap.set_yticks([])

    # Titles (adjust pad if needed)
    axes_dict["barplot"][0].set_title(
        "Prophage coordinates in bacterial contigs, bp",
        fontsize=TITLE_FONTSIZE,
        pad=20,  # Increased pad
    )
    axes_dict["taxonomy"][0].set_title(
        "Bacterial order and family",
        fontsize=TITLE_FONTSIZE,
        pad=20,  # Increased pad
    )
    axes_dict["genemap"][0].set_title(
        "Prophage genomic map",
        fontsize=TITLE_FONTSIZE,
        pad=20,  # Increased pad
    )

    # Adjust the figure size first (if needed)
    fig.set_size_inches(fig.get_size_inches()[0], fig.get_size_inches()[1] * 1.2)  # Increase height

    # === LEGENDS ===
    # === Crassvirales family legend ===
    # Use only families that are actually present in the plot; fall back gracefully if empty
    families_used_sorted = sorted(families_used) if families_used else sorted(family_to_color.keys())

    family_handles = [
        mpatches.Patch(color=family_to_color.get(f, DEFAULT_COLOR), label=f) for f in families_used_sorted
    ]

    # Create legend handles
    phylum_handles = [mpatches.Patch(color=used_phylum_color_map[p], label=p) for p in sorted(used_phylum_color_map)]
    class_handles = [mpatches.Patch(color=class_color_map[c], label=c) for c in sorted(class_color_map)]
    # order_handles = [mpatches.Patch(color=order_color_map[o], label=o) for o in sorted(order_color_map)]
    # === Add origin legend (Isolate/MAG) ===
    origin_handles = [
        mpatches.Patch(color=origin_color_map["Isolate"], label="Isolate"),
        mpatches.Patch(color=origin_color_map["MAG"], label="MAG"),
    ]

    # Position legends in the dedicated space

    # leg1 = fig.legend(
    #     order_handles,
    #     [h.get_label() for h in order_handles],
    #     ncol=4,
    #     loc="lower center",
    #     bbox_to_anchor=(0.25, -0.02),  # further down
    #     frameon=False,
    #     fontsize=TEXT_FONTSIZE,
    #     title="Bacterial order",
    #     title_fontsize=TEXT_FONTSIZE
    # )

    # Place it above the phylum/class legends so it doesn’t collide; tweak ncol/bbox as you like
    leg_crassvirales_family = fig.legend(
        family_handles,
        [h.get_label() for h in family_handles],
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.1, -0.06),
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        title="Crassvirales family",
        title_fontsize=LEGEND_TITLE_FONTSIZE,
        handleheight=LEGEND_HANDLE_HEIGHT,
        handlelength=LEGEND_HANDLE_LENGTH,
    )

    leg_origin = fig.legend(
        origin_handles,
        [h.get_label() for h in origin_handles],
        ncol=1,
        loc="lower left",
        bbox_to_anchor=(0.1, 0.03),
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        title="Genome origin",
        title_fontsize=LEGEND_TITLE_FONTSIZE,
        handleheight=LEGEND_HANDLE_HEIGHT,
        handlelength=LEGEND_HANDLE_LENGTH,
    )

    leg_phylum = fig.legend(
        phylum_handles,
        [h.get_label() for h in phylum_handles],
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.35, -0.1),
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        title="Bacterial phylum",
        title_fontsize=LEGEND_TITLE_FONTSIZE,
        handleheight=LEGEND_HANDLE_HEIGHT,
        handlelength=LEGEND_HANDLE_LENGTH,
    )

    leg_class = fig.legend(
        class_handles,
        [h.get_label() for h in class_handles],
        ncol=2,
        loc="lower right",
        bbox_to_anchor=(0.75, -0.1),
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        title="Bacterial class",
        title_fontsize=LEGEND_TITLE_FONTSIZE,
        handleheight=LEGEND_HANDLE_HEIGHT,
        handlelength=LEGEND_HANDLE_LENGTH,
    )

    # Add both legends to the axis
    # Register on the legend axis so it survives tight_layout
    axes_dict["legend"].add_artist(leg_crassvirales_family)
    axes_dict["legend"].add_artist(leg_origin)
    axes_dict["legend"].add_artist(leg_class)
    axes_dict["legend"].add_artist(leg_phylum)
    # axes_dict["legend"].add_artist(leg1)

    # Adjust layout with more space at bottom
    plt.subplots_adjust(
        bottom=0.15,  # Increased bottom space
        top=0.95,
        hspace=0.2,
    )

    fig.tight_layout(rect=[0.2, 0.15, 1, 0.85])
    fig.canvas.draw()

    adjust_left_margin_for_labels(fig, axes_dict["barplot"])

    for ax, label in [
        (axes_dict["family"][0], "Crassvirales family"),
        (axes_dict["is_mag"][0], "Genome origin"),
        (axes_dict["phylum"][0], "Bacterial phylum"),  # NEW
        (axes_dict["class"][0], "Bacterial class"),
        # (axes_dict["order"][0],  "Bacterial order"),
    ]:
        bbox = ax.get_position()
        x_mid = bbox.x0 + bbox.width / 2.0
        y_top = bbox.y1 + 0.02
        fig.text(x_mid, y_top, label, rotation=90, ha="center", va="bottom", fontsize=TITLE_FONTSIZE, clip_on=False)

    plt.savefig(output_file, dpi=300, bbox_inches="tight", pad_inches=0.5)  # Added pad_inches
    plt.savefig(Path(output_file).with_suffix(".svg"), bbox_inches="tight", pad_inches=0.5)
    plt.close()


def parse_itol_annotation(file_path: Path) -> tuple[dict[str, str], dict[str, str], set[str]]:
    leaf_to_family: dict[str, str] = {}
    color_to_family: dict[str, str] = {}
    family_to_color: dict[str, str] = {}
    crassvirales_families: set[str] = set()

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

    for color, label in zip(legend_colors, legend_labels, strict=False):
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


def parse_prophage_table(path: Path) -> tuple[dict[str, dict[str, str]], pd.DataFrame]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    prophage_dict: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        contig_id = row["contig_id"]
        prophage_dict[contig_id] = row.to_dict()  # 👈 ключ = contig_id
    return prophage_dict, df


def match_genome_id(leaf_name: str) -> str:
    return leaf_name.split("|")[0]


def annotate_tree_features(
    tree: Tree,
    leaf_to_family: dict[str, str],
    family_to_color: dict[str, str],
    prophage_dict: dict[str, dict[str, str]],
    leaf_metadata: dict[str, dict[str, str]],
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
                contig_length=int(matched_prophage["contig_length"]),
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


def get_custom_layout(family_to_color: dict[str, str]):
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
                    node.add_face(
                        ProphageBarFace(contig_len, start, end, width=100, height=10), column=1, position="aligned"
                    )

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


def render_tree(tree: Tree, output_file: Path, family_to_color: dict[str, str]) -> None:
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
    prophage_table = Path(
        f"{wd}/prophage_alignments_ncbi_and_gtdb_samtools_coordinates_edited_with_lengths_and_mags.tsv"
    )
    taxonomy_file = Path(
        f"{wd}/prophage_alignments_ncbi_and_gtdb_samtools_coordinates_edited_with_lengths_bacterial_contig_ids_taxonomy.tsv"
    )
    functional_annotation_file = Path(f"{crassus_output}/4_ORF/3_functional_annot_table_renamed.tsv")

    output_svg = Path(f"{wd}/TerL_with_prophages_annotated_tree.svg")
    prophage_plot = Path(f"{wd}/Crassphage_prophage_analysis_plot.png")

    tree = read_tree(tree_file)
    leaf_to_family, family_to_color, _ = parse_itol_annotation(itol_annotation)

    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", dtype=str)

    feature_dict = load_functional_annotation(str(functional_annotation_file))

    # print(f"{feature_dict=}")

    # NEW: propagate family annotation by MRCA
    propagate_family_annotations_by_mrca(tree, leaf_to_family, family_to_color)

    # === Build consistent metadata mapping for each leaf ===
    # leaf_metadata = build_leaf_metadata(tree)
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

    # print("🔍 First 5 prophage_ids:")
    # print(prophage_dict.keys()[:5])

    # New: plot genomic positions
    # plot_prophages(prophage_df, output_file=str(prophage_plot))
    # plot_prophage_positions(prophage_dict, output_file=str(prophage_plot))

    compare_genomes_between_prophage_and_functional_tables(
        prophage_df=prophage_df,  # ✅ Already read as DataFrame
        functional_df=pd.read_csv(functional_annotation_file, sep="\t"),  # ✅ Explicit load
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
        feature_dict=feature_dict,
    )

    print(f"Prophage plot saved to: {prophage_plot}")


if __name__ == "__main__":
    main()
