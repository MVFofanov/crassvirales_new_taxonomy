#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
import pandas as pd


# ================================
# Color schemes
# ================================

PHYLUM_COLORS = {
    "Bacteroidota":   "#33a02c",
    "Bacillota":      "#1f78b4",
    "Pseudomonadota": "#ff7f00",
    "Actinomycetota": "#6A3D9A",
    "Bacteria":       "#FFD700",
    "Other":          "#B3B3B3",
}

CLASS_COLORS = {
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

INTEGRATION_BAR_COLORS = {
    "integrated": "#E31A1C",      # red
    "non_integrated": "#1F78B4",  # blue
}

HEX_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate multiple iTOL annotation layers from TSV."
    )
    parser.add_argument("--input_tsv", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--show_strip_labels", action="store_true")
    return parser.parse_args()


def validate_colors(color_dict):
    bad = {k: v for k, v in color_dict.items() if not HEX_RE.match(v)}
    if bad:
        raise ValueError(
            "Invalid colors:\n" +
            "\n".join(f"{k}: {v}" for k, v in bad.items())
        )


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
    return v if v in PHYLUM_COLORS else "Other"


def normalize_class(v):
    if pd.isna(v):
        return "Other"

    v = str(v).strip()
    return v if v in CLASS_COLORS else "Other"


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
        "MARGIN\t5",
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
        "MARGIN\t5",
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
    "bact_phylum": {
        "colors": PHYLUM_COLORS,
        "normalize": normalize_phylum,
        "dataset_label": "Bacterial phylum",
        "legend_title": "Bacterial phylum",
        "legend_order": [
            "Bacteroidota",
            "Bacillota",
            "Pseudomonadota",
            "Actinomycetota",
            "Bacteria",
            "Other",
        ],
        "suffix": "bacterial_phylum",
    },
    "bact_class": {
        "colors": CLASS_COLORS,
        "normalize": normalize_class,
        "dataset_label": "Bacterial class",
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
        "suffix": "bacterial_class",
    },
}


def main():
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = out_dir.name
    print(f"[INFO] Using prefix: {prefix}")

    validate_colors(PHYLUM_COLORS)
    validate_colors(CLASS_COLORS)
    validate_colors(INTEGRATION_BAR_COLORS)

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
    for col, cfg in CONFIGS.items():
        if col not in df.columns:
            print(f"[SKIP] Column {col} not found")
            continue

        print(f"[INFO] Processing {col}")

        sub = df[["label", col]].copy()
        sub[f"{col}_norm"] = sub[col].apply(cfg["normalize"])

        out_text = build_itol_colorstrip(
            sub,
            column=col,
            config=cfg,
            show_labels=args.show_strip_labels
        )

        out_file = out_dir / f"{prefix}_{cfg['suffix']}.txt"
        out_file.write_text(out_text, encoding="utf-8")

        print(f"[OK] Wrote: {out_file}")
        print(sub[f"{col}_norm"].value_counts().to_string(), "\n")


if __name__ == "__main__":
    main()