#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

from ete3 import Tree
import pandas as pd


# ================================
# Color schemes
# ================================
CRASSVIRALES_FAMILY_COLORS = {
    "Intestiviridae": "#EE3B3B",
    "Crevaviridae":   "#EE9A00",
    "Suoliviridae":   "#4169E1",
    "Steigviridae":   "#00CED1",
    "Epsilon":        "#CD2990",
    "Zeta":           "#006400",
    "Outgroup":       "#EE82EE",
    "Other":          "#000000",
}

PHYLUM_COLORS_DIFFERENT = {
    "Bacteroidota":   "#33a02c",
    "Bacillota":      "#1f78b4",
    "Pseudomonadota": "#ff7f00",
    "Actinomycetota": "#6A3D9A",
    "Bacteria":       "#FFD700",
    "Other":          "#B3B3B3",
}

CLASS_COLORS_DIFFERENT = {
    "Flavobacteriia":      "#b2df8a",
    "Bacteroidia":         "#a6cee3",
    "Cytophagia":          "#ffff99",
    "Saprospiria":         "#cab2d6",
    "Chitinophagia":       "#6a3d9a",
    "Clostridia":          "#b15928",
    "Alphaproteobacteria": "#E31A1C",
    "Gammaproteobacteria": "#FB9A99",
    "Betaproteobacteria":  "#fdbf6f",
    "Other":               "#B3B3B3",
}

PHYLUM_COLORS_SIMILAR = {
    "Bacteroidota":   "#1B9E77",  # green
    "Bacillota":      "#7570B3",  # purple-blue
    "Pseudomonadota": "#D95F02",  # orange
    "Actinomycetota": "#E7298A",  # magenta
    "Bacteria":       "#666666",  # dark grey
    "Other":          "#BDBDBD",  # light grey
}

CLASS_COLORS_SIMILAR = {
    # Bacteroidota — green/teal shades
    "Bacteroidia":    "#006D2C",  # dark green
    "Flavobacteriia": "#31A354",  # medium green
    "Cytophagia":     "#74C476",  # light green
    "Saprospiria":    "#238B8D",  # teal
    "Chitinophagia":  "#00441B",  # very dark green

    # Bacillota — purple/blue shades
    "Clostridia":     "#54278F",  # dark purple

    # Pseudomonadota — orange/red shades
    "Alphaproteobacteria": "#D95F02",  # orange
    "Gammaproteobacteria": "#E6550D",  # red-orange
    "Betaproteobacteria":  "#Fdae6b",  # light orange

    # fallback
    "Other": "#BDBDBD",
}

INTEGRATION_BAR_COLORS = {
    "integrated": "#E31A1C",      # red
    "non_integrated": "#1F78B4",  # blue
}

SOURCE_TYPE_COLORS = {
    "isolate": "#000000",
    "MAG": "#000000",
}

COMPLETENESS_COLORS = {
    "75-100": "#33A02C",  # green
    "50-74":  "#FF7F00",  # orange
    "25-49":  "#E31A1C",  # red
    "0-24":   "#E31A1C",  # red
}

COMMON_MARGIN = 5
SYMBOL_SIZE = 15

HEX_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate multiple iTOL annotation layers from TSV."
    )
    parser.add_argument("--input_tsv", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--show_strip_labels", action="store_true")
    parser.add_argument("--tree_file", required=False,
                    help="Tree file used to identify MRCA nodes for Crassvirales families.")
    parser.add_argument("--reference_taxonomy_tsv", required=False,
                        help="TSV with reference Crassvirales leaf taxonomy.")
    return parser.parse_args()


def validate_colors(color_dict):
    bad = {k: v for k, v in color_dict.items() if not HEX_RE.match(v)}
    if bad:
        raise ValueError(
            "Invalid colors:\n" +
            "\n".join(f"{k}: {v}" for k, v in bad.items())
        )

def normalize_source_type(v):
    if pd.isna(v):
        return None
    v = str(v).strip()
    if v in {"MAG", "isolate"}:
        return v
    return None


def normalize_completeness_group(v):
    if pd.isna(v):
        return None
    v = str(v).strip()
    if v in {"75-100", "50-74", "25-49", "0-24"}:
        return v
    return None

# ================================
# Normalization
# ================================

def normalize_phylum(v):
    if pd.isna(v):
        return "Other"

    v = str(v).strip()

    synonym_map = {
        "Proteobacteria": "Pseudomonadota",
        "Firmicutes": "Bacillota",
        "Actinobacteria": "Actinomycetota",
        "Actinobacteriota": "Actinomycetota",
    }

    v = synonym_map.get(v, v)
    allowed = set(PHYLUM_COLORS_DIFFERENT) | set(PHYLUM_COLORS_SIMILAR)

    return v if v in allowed else "Other"


def normalize_class(v):
    if pd.isna(v):
        return "Other"

    v = str(v).strip()
    allowed = set(CLASS_COLORS_DIFFERENT) | set(CLASS_COLORS_SIMILAR)

    return v if v in allowed else "Other"


def normalize_integration_status(v):
    if pd.isna(v):
        return None
    v = str(v).strip()
    return v if v in INTEGRATION_BAR_COLORS else None


def parse_length_kb(v):
    if pd.isna(v):
        return None
    v = str(v).strip()
    if v == "":
        return None
    try:
        return float(v) / 1000.0   # convert nt to kb
    except ValueError:
        return None


# ================================
# Builders
# ================================

def build_itol_colorstrip(df, column, config, show_labels):
    colors = config["colors"]
    legend_order = config["legend_order"]

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{config['dataset_label']}",
        "COLOR\t#999999",
        "STRIP_WIDTH\t25",
        f"MARGIN\t{COMMON_MARGIN}",
        "BORDER_WIDTH\t1",
        "BORDER_COLOR\t#FFFFFF",
        "SHOW_INTERNAL\t0",
        f"SHOW_STRIP_LABELS\t{1 if show_labels else 0}",
        "STRIP_LABEL_AUTO_COLOR\t1",
        f"LEGEND_TITLE\t{config['legend_title']}",
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_order)),
        "LEGEND_COLORS\t" + "\t".join(colors[x] for x in legend_order),
        "LEGEND_LABELS\t" + "\t".join(legend_order),
        "DATA",
    ]

    for _, row in df.iterrows():
        node_id = row["label"]
        category = row[f"{column}_norm"]
        color = colors[category]
        lines.append(f"{node_id}\t{color}\t{category}")

    return "\n".join(lines) + "\n"


def build_itol_multibar_length(df):
    """
    Build one MULTIBAR dataset with two series:
    integrated (red) and non_integrated (blue).
    Each leaf gets value in only one of the two columns.
    """

    lines = [
        "DATASET_MULTIBAR",
        "SEPARATOR TAB",
        "DATASET_LABEL\tProphage length (kb)",
        "COLOR\t#999999",
        "FIELD_LABELS\tintegrated\tnon_integrated",
        "FIELD_COLORS\t#E31A1C\t#1F78B4",
        "LEGEND_TITLE\tIntegration status",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#E31A1C\t#1F78B4",
        "LEGEND_LABELS\tintegrated\tnon_integrated",
        "WIDTH\t200",
        f"MARGIN\t{COMMON_MARGIN}",
        "BORDER_WIDTH\t0",
        "DATASET_SCALE\t0-0 kb-#999999-1-0-1\t50-50 kb-#999999-1-0-1\t100-100 kb-#999999-1-0-1\t150-150 kb-#999999-1-0-1\t200-200 kb-#999999-1-0-1",
        "DATA",
    ]

    for _, row in df.iterrows():
        node_id = row["label"]
        value_kb = row["bact_genomad_length_num_kb"]
        status = row["integration_status_norm"]

        if status == "integrated":
            integrated_val = value_kb
            non_integrated_val = 0
        elif status == "non_integrated":
            integrated_val = 0
            non_integrated_val = value_kb
        else:
            continue

        lines.append(f"{node_id}\t{integrated_val:.3f}\t{non_integrated_val:.3f}")

    return "\n".join(lines) + "\n"


# ================================
# Config definitions
# ================================

CONFIGS = {
    "bact_phylum_different": {
        "source_column": "bact_phylum",
        "colors": PHYLUM_COLORS_DIFFERENT,
        "normalize": normalize_phylum,
        "dataset_label": "Bacterial phylum different",
        "legend_title": "Bacterial phylum",
        "legend_order": [
            "Bacteroidota",
            "Bacillota",
            "Pseudomonadota",
            "Actinomycetota",
            "Bacteria",
            "Other",
        ],
        "suffix": "bacterial_phylum_different",
    },
    "bact_phylum_similar": {
        "source_column": "bact_phylum",
        "colors": PHYLUM_COLORS_SIMILAR,
        "normalize": normalize_phylum,
        "dataset_label": "Bacterial phylum similar",
        "legend_title": "Bacterial phylum",
        "legend_order": [
            "Bacteroidota",
            "Bacillota",
            "Pseudomonadota",
            "Actinomycetota",
            "Bacteria",
            "Other",
        ],
        "suffix": "bacterial_phylum_similar",
    },
    "bact_class_different": {
        "source_column": "bact_class",
        "colors": CLASS_COLORS_DIFFERENT,
        "normalize": normalize_class,
        "dataset_label": "Bacterial class different",
        "legend_title": "Bacterial class",
        "legend_order": [
            "Flavobacteriia",
            "Bacteroidia",
            "Cytophagia",
            "Saprospiria",
            "Chitinophagia",
            "Clostridia",
            "Alphaproteobacteria",
            "Gammaproteobacteria",
            "Betaproteobacteria",
            "Other",
        ],
        "suffix": "bacterial_class_different",
    },
    "bact_class_similar": {
        "source_column": "bact_class",
        "colors": CLASS_COLORS_SIMILAR,
        "normalize": normalize_class,
        "dataset_label": "Bacterial class similar",
        "legend_title": "Bacterial class",
        "legend_order": [
            "Flavobacteriia",
            "Bacteroidia",
            "Cytophagia",
            "Saprospiria",
            "Chitinophagia",
            "Clostridia",
            "Alphaproteobacteria",
            "Gammaproteobacteria",
            "Betaproteobacteria",
            "Other",
        ],
        "suffix": "bacterial_class_similar",
    },
}

def normalize_crass_family(v):
    if pd.isna(v):
        return "Other"
    v = str(v).strip()
    return v if v in CRASSVIRALES_FAMILY_COLORS else "Other"


def build_itol_tree_colors(tree, ref_df, mode="clade"):
    """
    Build TREE_COLORS dataset for Crassvirales family highlights.

    mode:
        - "clade" : color actual tree branches
        - "range" : add outer clade/range highlight
    """

    if mode not in {"clade", "range"}:
        raise ValueError(f"Unsupported TREE_COLORS mode: {mode}")

    lines = [
        "TREE_COLORS",
        "SEPARATOR TAB",
        "DATA",
    ]

    family_counts = {}

    for family, sub in ref_df.groupby("family_norm"):
        leaves = sub["leaf_label"].dropna().astype(str).str.strip().tolist()
        leaves = [x for x in leaves if x]

        if len(leaves) == 0:
            continue

        present = [x for x in leaves if tree.search_nodes(name=x)]
        missing_n = len(leaves) - len(present)

        if len(present) == 0:
            print(f"[WARN] No leaves from family '{family}' found in tree")
            continue

        if len(present) == 1:
            node = tree.search_nodes(name=present[0])[0]
        else:
            node = tree.get_common_ancestor(present)

        color = CRASSVIRALES_FAMILY_COLORS.get(
            family,
            CRASSVIRALES_FAMILY_COLORS["Other"]
        )

        if mode == "clade":
            # node_id  type   color   style   width
            lines.append(f"{node.name}\tclade\t{color}\tnormal\t2")
        elif mode == "range":
            # node_id  type   color   label
            lines.append(f"{node.name}\trange\t{color}\t{family}")

        family_counts[family] = {
            "total_leaves": len(leaves),
            "present_in_tree": len(present),
            "missing_from_tree": missing_n,
            "node_id": node.name,
        }

    return "\n".join(lines) + "\n", family_counts

def assign_unique_internal_node_ids(tree, prefix="NODE"):
    """
    Assign unique names to ALL internal nodes.
    Existing internal labels (e.g. bootstrap strings like 100/100) are replaced.
    """
    counter = 1
    for node in tree.traverse("preorder"):
        if not node.is_leaf():
            node.name = f"{prefix}_{counter}"
            counter += 1
    return tree


def build_itol_binary_source_type(all_tree_leaves, prophage_df):
    """
    Build one DATASET_BINARY file for source_type.

    Mapping:
      isolate -> 1  (filled square)
      MAG     -> 0  (empty square)
      missing -> -1 (omitted completely)
    """

    # map label -> normalized source type from prophage table
    source_map = (
        prophage_df[["label", "source_type_norm"]]
        .drop_duplicates(subset=["label"])
        .set_index("label")["source_type_norm"]
        .to_dict()
    )

    lines = [
        "DATASET_BINARY",
        "SEPARATOR TAB",
        "DATASET_LABEL\tSource type",
        "COLOR\t#000000",
        "FIELD_SHAPES\t1",
        "FIELD_LABELS\tsource_type",
        "FIELD_COLORS\t#000000",
        "LEGEND_TITLE\tSource type",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#000000\t#000000",
        "LEGEND_LABELS\tisolate (filled)\tMAG (empty)",
        "DATA",
    ]

    counts = {"isolate": 0, "MAG": 0, "omitted": 0}

    for leaf in all_tree_leaves:
        src = source_map.get(leaf, None)

        if src == "isolate":
            value = 1
            counts["isolate"] += 1
        elif src == "MAG":
            value = 0
            counts["MAG"] += 1
        else:
            value = -1
            counts["omitted"] += 1

        lines.append(f"{leaf}\t{value}")

    return "\n".join(lines) + "\n", counts


def build_itol_source_type_symbols(all_tree_leaves, df):
    source_map = (
        df[["label", "source_type_norm"]]
        .drop_duplicates("label")
        .set_index("label")["source_type_norm"]
        .to_dict()
    )

    lines = [
        "DATASET_SYMBOL",
        "SEPARATOR TAB",
        "DATASET_LABEL\tSource type",
        "COLOR\t#000000",
        "LEGEND_TITLE\tSource type",
        "LEGEND_SHAPES\t1\t1",
        "LEGEND_COLORS\t#000000\t#000000",
        "LEGEND_LABELS\tisolate\tMAG",
        "LEGEND_SHAPE_INVERT\t0\t1",
        "DATA",
    ]

    counts = {"isolate": 0, "MAG": 0, "omitted": 0}

    for leaf in all_tree_leaves:
        src = source_map.get(leaf, None)

        if src == "isolate":
            # filled square
            lines.append(f"{leaf}\t1\t{SYMBOL_SIZE}\t#000000\t1\t-1\t")
            counts["isolate"] += 1
        elif src == "MAG":
            # empty square
            lines.append(f"{leaf}\t1\t{SYMBOL_SIZE}\t#000000\t0\t-1\t")
            counts["MAG"] += 1
        else:
            counts["omitted"] += 1
            # no entry -> omitted

    return "\n".join(lines) + "\n", counts

def build_itol_completeness_symbols(df, external_position=-2, symbol_size=SYMBOL_SIZE):
    """
    Build one DATASET_SYMBOL annotation for completeness_group.

    Circle colors:
      75-100 -> green
      50-74  -> orange
      25-49  -> red
      0-24   -> red
    """

    lines = [
        "DATASET_SYMBOL",
        "SEPARATOR TAB",
        "DATASET_LABEL\tCompleteness group",
        "COLOR\t#999999",
        "LEGEND_TITLE\tCompleteness group",
        "LEGEND_SHAPES\t2\t2\t2",
        "LEGEND_COLORS\t#33A02C\t#FF7F00\t#E31A1C",
        "LEGEND_LABELS\t75-100\t50-74\t0-49",
        "DATA",
    ]

    counts = {"75-100": 0, "50-74": 0, "25-49": 0, "0-24": 0}

    for _, row in df.iterrows():
        node_id = row["label"]
        group = row["completeness_group_norm"]

        if group is None:
            continue

        color = COMPLETENESS_COLORS[group]
        counts[group] += 1

        # format:
        # node_id  shape  size  color  fill  position  label
        # shape 2 = circle
        # fill 1 = filled
        lines.append(f"{node_id}\t2\t{SYMBOL_SIZE}\t{color}\t1\t{external_position}\t{group}")

    return "\n".join(lines) + "\n", counts


def parse_completeness_percent(v):
    if pd.isna(v):
        return None
    v = str(v).strip().replace("%", "")
    if v == "":
        return None
    try:
        x = float(v)
    except ValueError:
        return None
    return max(0.0, min(100.0, x))


def interpolate_hex(c0, c1, t):
    c0 = c0.lstrip("#")
    c1 = c1.lstrip("#")

    r0, g0, b0 = int(c0[0:2], 16), int(c0[2:4], 16), int(c0[4:6], 16)
    r1, g1, b1 = int(c1[0:2], 16), int(c1[2:4], 16), int(c1[4:6], 16)

    r = round(r0 + (r1 - r0) * t)
    g = round(g0 + (g1 - g0) * t)
    b = round(b0 + (b1 - b0) * t)

    return f"#{r:02X}{g:02X}{b:02X}"


def completeness_to_color(x):
    # viridis-like gradient (very distinguishable)
    stops = [
        (0,   "#440154"),
        (25,  "#3B528B"),
        (50,  "#21918C"),
        (75,  "#5EC962"),
        (100, "#FDE725"),
    ]

    for (x0, c0), (x1, c1) in zip(stops[:-1], stops[1:]):
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0)
            return interpolate_hex(c0, c1, t)

    return stops[-1][1]

def build_itol_completeness_tip_symbols(df, position=1, symbol_size=SYMBOL_SIZE):
    """
    Draw completeness as filled circles at leaf tips using continuous colors.
    Legend shows reference points along the gradient.
    """

    lines = [
        "DATASET_SYMBOL",
        "SEPARATOR TAB",
        "DATASET_LABEL\tCompleteness",
        "COLOR\t#999999",

        # Manual reference legend for continuous gradient
        "LEGEND_TITLE\tCompleteness (%)",
        "LEGEND_SHAPES\t2\t2\t2\t2\t2",
        "LEGEND_COLORS\t#440154\t#3B528B\t#21918C\t#5EC962\t#FDE725",
        "LEGEND_LABELS\t0\t25\t50\t75\t100",

        "DATA",
    ]

    counts = {"written": 0, "skipped": 0}

    for _, row in df.iterrows():
        node_id = row["label"]
        completeness = row["completeness_percent"]

        if pd.isna(completeness):
            counts["skipped"] += 1
            continue

        color = completeness_to_color(completeness)

        lines.append(
            f"{node_id}\t2\t{symbol_size}\t{color}\t1\t{position}\t{completeness:.1f}%"
        )

        counts["written"] += 1

    return "\n".join(lines) + "\n", counts

def main():
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = out_dir.name
    print(f"[INFO] Using prefix: {prefix}")

    validate_colors(PHYLUM_COLORS_DIFFERENT)
    validate_colors(PHYLUM_COLORS_SIMILAR)
    validate_colors(CLASS_COLORS_DIFFERENT)
    validate_colors(CLASS_COLORS_SIMILAR)
    validate_colors(INTEGRATION_BAR_COLORS)
    validate_colors(CRASSVIRALES_FAMILY_COLORS)
    validate_colors(SOURCE_TYPE_COLORS)
    validate_colors(COMPLETENESS_COLORS)

    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str)

    df = df.dropna(subset=["label"])
    df["label"] = df["label"].astype(str).str.strip()
    df = df[df["label"] != ""]

    print(f"[INFO] Loaded rows: {len(df)}")

    if df["label"].duplicated().any():
        print("[WARN] Duplicate labels found, keeping first")
        df = df.drop_duplicates(subset=["label"])

    # --------------------------------
    # 1) First layer: prophage length
    # --------------------------------
    required_bar_cols = {"label", "bact_genomad_length_num", "integration_status"}
    if required_bar_cols.issubset(df.columns):
        bar_df = df[["label", "bact_genomad_length_num", "integration_status"]].copy()
        bar_df["bact_genomad_length_num_kb"] = bar_df["bact_genomad_length_num"].apply(parse_length_kb)
        bar_df["integration_status_norm"] = bar_df["integration_status"].apply(normalize_integration_status)

        before = len(bar_df)
        bar_df = bar_df.dropna(subset=["bact_genomad_length_num_kb", "integration_status_norm"])
        after = len(bar_df)

        if after < before:
            print(f"[WARN] Skipped {before - after} rows for prophage length layer due to missing/invalid values")

        out_text = build_itol_multibar_length(bar_df)
        out_file = out_dir / f"{prefix}_prophage_length.txt"
        out_file.write_text(out_text, encoding="utf-8")
        print(f"[OK] Wrote: {out_file}")

        print("[INFO] Integration status counts:")
        print(bar_df["integration_status_norm"].value_counts().to_string(), "\n")
    else:
        missing = sorted(required_bar_cols - set(df.columns))
        print(f"[SKIP] Prophage length layer skipped; missing columns: {missing}")

    # --------------------------------
    # 2+) Color strip layers
    # --------------------------------
    for config_name, cfg in CONFIGS.items():
        source_col = cfg["source_column"]

        if source_col not in df.columns:
            print(f"[SKIP] Column {source_col} not found for {config_name}")
            continue

        print(f"[INFO] Processing {config_name} from column {source_col}")

        sub = df[["label", source_col]].copy()
        sub[f"{source_col}_norm"] = sub[source_col].apply(cfg["normalize"])

        out_text = build_itol_colorstrip(
            sub,
            column=source_col,
            config=cfg,
            show_labels=args.show_strip_labels
        )

        out_file = out_dir / f"{prefix}_{cfg['suffix']}.txt"
        out_file.write_text(out_text, encoding="utf-8")

        print(f"[OK] Wrote: {out_file}")
        print(sub[f"{source_col}_norm"].value_counts().to_string(), "\n")


    # --------------------------------
    # 3) Crassvirales family clade highlights
    # --------------------------------
    if args.tree_file and args.reference_taxonomy_tsv:
        print("[INFO] Processing Crassvirales family clade highlights")

        tree = Tree(args.tree_file, format=1)
        tree = assign_unique_internal_node_ids(tree)

        named_tree_file = out_dir / f"{prefix}_with_internal_ids.treefile"
        tree.write(outfile=str(named_tree_file), format=1)
        print(f"[OK] Wrote tree with internal node IDs: {named_tree_file}")

        ref_df = pd.read_csv(args.reference_taxonomy_tsv, sep="\t", dtype=str)

        required_ref_cols = {"leaf_label", "family"}
        missing_ref = required_ref_cols - set(ref_df.columns)
        if missing_ref:
            raise ValueError(f"Missing required columns in reference taxonomy TSV: {sorted(missing_ref)}")

        ref_df = ref_df[["leaf_label", "family"]].copy()
        ref_df["leaf_label"] = ref_df["leaf_label"].astype(str).str.strip()
        ref_df = ref_df[ref_df["leaf_label"] != ""]
        ref_df["family_norm"] = ref_df["family"].apply(normalize_crass_family)

                # --- clade coloring file ---
        tree_colors_clade_text, family_counts = build_itol_tree_colors(
            tree,
            ref_df,
            mode="clade"
        )

        out_file = out_dir / f"{prefix}_crassvirales_family_clades.txt"
        out_file.write_text(tree_colors_clade_text, encoding="utf-8")
        print(f"[OK] Wrote: {out_file}")

        # --- range highlighting file ---
        tree_colors_range_text, _ = build_itol_tree_colors(
            tree,
            ref_df,
            mode="range"
        )

        out_file = out_dir / f"{prefix}_crassvirales_family_ranges.txt"
        out_file.write_text(tree_colors_range_text, encoding="utf-8")
        print(f"[OK] Wrote: {out_file}")

        print("[INFO] Family clade summary:")
        for fam, info in family_counts.items():
            print(
                f"  {fam}: total={info['total_leaves']}, "
                f"present={info['present_in_tree']}, missing={info['missing_from_tree']}, "
                f"node_id={info['node_id']}"
            )
        print()
    
        # --------------------------------
    # Source type binary layer
    # --------------------------------
    if "source_type" in df.columns and args.tree_file:
        print("[INFO] Processing source_type symbols layer")

        tree_for_leaves = Tree(args.tree_file, format=1)
        all_tree_leaves = [leaf.name for leaf in tree_for_leaves.iter_leaves()]

        source_df = df[["label", "source_type"]].copy()
        source_df["source_type_norm"] = source_df["source_type"].apply(normalize_source_type)

        out_text, source_counts = build_itol_source_type_symbols(all_tree_leaves, source_df)

        out_file = out_dir / f"{prefix}_source_type_symbols.txt"
        out_file.write_text(out_text, encoding="utf-8")
        print(f"[OK] Wrote: {out_file}")

        print("[INFO] Source type counts:")
        for k, v in source_counts.items():
            print(f"  {k}: {v}")
        print()
    elif "source_type" not in df.columns:
        print("[SKIP] source_type binary layer skipped; column 'source_type' not found")
    elif not args.tree_file:
        print("[SKIP] source_type binary layer skipped; --tree_file is required to include -1 for missing/reference leaves")

    # --------------------------------
    # Completeness tip circles
    # --------------------------------
    completeness_col = "completeness"

    if completeness_col in df.columns:
        print("[INFO] Processing completeness tip-circle layer")

        comp_df = df[["label", completeness_col]].copy()
        comp_df["completeness_percent"] = comp_df[completeness_col].apply(parse_completeness_percent)

        before = len(comp_df)
        comp_df = comp_df.dropna(subset=["completeness_percent"])
        after = len(comp_df)

        if after < before:
            print(f"[WARN] Skipped {before - after} rows for completeness due to missing/invalid values")

        out_text, comp_counts = build_itol_completeness_tip_symbols(
            comp_df,
            position=1,
            symbol_size=SYMBOL_SIZE
        )

        out_file = out_dir / f"{prefix}_completeness_tip_circles.txt"
        out_file.write_text(out_text, encoding="utf-8")
        print(f"[OK] Wrote: {out_file}")

        print("[INFO] Completeness symbols:")
        for k, v in comp_counts.items():
            print(f"  {k}: {v}")
        print()

    else:
        print(f"[SKIP] completeness tip-circle layer skipped; column '{completeness_col}' not found")

if __name__ == "__main__":
    main()