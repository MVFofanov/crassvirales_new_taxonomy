#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import fisher_exact


# =========================
# CLI
# =========================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate multipanel boxplots and stats for integration status comparison"
    )
    parser.add_argument(
        "--annotation_table",
        required=True,
        help="Input TSV table"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--prophages_to_exclude",
        default=None,
        help="Optional file containing labels of prophages to exclude (one per line)"
        )
    parser.add_argument(
        "--poster",
        action="store_true",
        help="Generate simplified poster version (2 panels)."
    )
    return parser.parse_args()


# =========================
# Stats
# =========================
def cliffs_delta(x, y):
    import numpy as np

    x = np.array(x)
    y = np.array(y)

    n = 0
    for xi in x:
        n += np.sum(xi > y) - np.sum(xi < y)

    delta = n / (len(x) * len(y))
    return delta


def compute_stats(df, col):
    g1 = df[df["integration_status"] == "integrated"][col]
    g2 = df[df["integration_status"] == "non_integrated"][col]

    stat, p = mannwhitneyu(g1, g2, alternative="two-sided")

    # Quartiles
    q1_int = g1.quantile(0.25)
    q3_int = g1.quantile(0.75)
    iqr_int = q3_int - q1_int

    q1_non = g2.quantile(0.25)
    q3_non = g2.quantile(0.75)
    iqr_non = q3_non - q1_non

    delta = cliffs_delta(g1, g2)

    return {
        "variable": col,

        "n_integrated": len(g1),
        "n_non_integrated": len(g2),

        "median_integrated": g1.median(),
        "q1_integrated": q1_int,
        "q3_integrated": q3_int,
        "iqr_integrated": iqr_int,

        "median_non_integrated": g2.median(),
        "q1_non_integrated": q1_non,
        "q3_non_integrated": q3_non,
        "iqr_non_integrated": iqr_non,

        "p_value": p,

        "cliffs_delta": delta
    }

def source_type_composition(df):
    tab = (
        df.groupby(["integration_status", "source_type"])
        .size()
        .reset_index(name="n")
    )

    tab["percent_within_integration_status"] = (
        tab["n"] /
        tab.groupby("integration_status")["n"].transform("sum") * 100
    )

    return tab

def above_median_summary(df, cols):
    results = []

    for col in cols:
        global_median = df[col].median()

        tmp = df.copy()
        tmp["above_global_median"] = tmp[col] > global_median

        summary = (
            tmp.groupby(["integration_status", "source_type", "above_global_median"])
            .size()
            .reset_index(name="n")
        )

        total = (
            tmp.groupby(["integration_status", "source_type"])
            .size()
            .reset_index(name="total")
        )

        summary = summary.merge(
            total,
            on=["integration_status", "source_type"],
            how="left"
        )

        summary["percent"] = summary["n"] / summary["total"] * 100
        summary["variable"] = col
        summary["global_median"] = global_median

        results.append(summary)

    return pd.concat(results, ignore_index=True)


def source_type_metric_tests(df, cols):
    results = []

    for status in ["integrated", "non_integrated"]:
        sub = df[df["integration_status"] == status].copy()

        for col in cols:
            iso = sub[sub["source_type"] == "isolate"][col].dropna()
            mag = sub[sub["source_type"] == "MAG"][col].dropna()

            if len(iso) >= 2 and len(mag) >= 2:
                stat, p = mannwhitneyu(iso, mag, alternative="two-sided")
                delta = cliffs_delta(iso, mag)
            else:
                stat, p, delta = np.nan, np.nan, np.nan

            results.append({
                "integration_status": status,
                "variable": col,

                "n_isolate": len(iso),
                "n_MAG": len(mag),

                "median_isolate": iso.median() if len(iso) else np.nan,
                "q1_isolate": iso.quantile(0.25) if len(iso) else np.nan,
                "q3_isolate": iso.quantile(0.75) if len(iso) else np.nan,

                "median_MAG": mag.median() if len(mag) else np.nan,
                "q1_MAG": mag.quantile(0.25) if len(mag) else np.nan,
                "q3_MAG": mag.quantile(0.75) if len(mag) else np.nan,

                "p_value": p,
                "cliffs_delta_isolate_vs_MAG": delta
            })

    return pd.DataFrame(results)

# =========================
# Plot
# =========================
def make_plot(df, outdir, poster=False):
    order = ["integrated", "non_integrated"]

    palette = {
        "integrated": "#EE3B3B",
        "non_integrated": "#4169E1"
    }

    sns.set_theme(style="white", context="talk")


    plt.rcParams.update({
        "font.size": 18,
        "axes.labelsize": 18,
        "axes.titlesize": 18,
        "xtick.labelsize": 16,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
    })

    if poster:
        fig, axes = plt.subplots(
            nrows=1,
            ncols=2,
            figsize=(15, 5)
        )
    else:
        fig, axes = plt.subplots(
            nrows=3,
            ncols=1,
            figsize=(5.5, 10.5)
        )

    def draw(ax, y_col, ylabel, label, log=False):
        # =========================
        # Boxplot (NO FILL)
        # =========================
        sns.boxplot(
            data=df,
            x="integration_status",
            y=y_col,
            order=order,
            showcaps=True,
            showfliers=False,
            width=0.5,
            boxprops=dict(facecolor="none", linewidth=1.5),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5),
            medianprops=dict(color="black", linewidth=1.5),
            ax=ax
        )

        # Manually recolor box edges
        for i, artist in enumerate(ax.artists):
            group = order[i]
            color = palette[group]
            artist.set_edgecolor(color)

        # =========================
        # Points
        # =========================
        import numpy as np
        rng = np.random.default_rng(42)

        x_positions = {group: i for i, group in enumerate(order)}

        for _, row in df.iterrows():
            group = row["integration_status"]
            source = row["source_type"]
            y = row[y_col]

            x = x_positions[group] + rng.uniform(-0.18, 0.18)
            color = palette[group]

            if source == "isolate":
                ax.scatter(
                    x, y,
                    s=30,
                    facecolors=color,
                    edgecolors="black",
                    linewidths=0.5,
                    alpha=0.9,
                    zorder=3
                )
            elif source == "MAG":
                ax.scatter(
                    x, y,
                    s=30,
                    facecolors="none",
                    edgecolors=color,
                    linewidths=1.2,
                    alpha=0.9,
                    zorder=3
                )

        # =========================
        # Axes
        # =========================
        if log:
            ax.set_yscale("log")

        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Integrated", "Non-integrated"])

        ax.text(
            -0.18, 1.05,
            label,
            transform=ax.transAxes,
            fontsize=20,
            fontweight="bold"
        )

        sns.despine(ax=ax)

    if poster:
        draw(
            axes[0],
            "completeness",
            "CheckV completeness (%)",
            "A",
            log=False
        )

        axes[0].set_ylim(-3, 103)

        draw(
            axes[1],
            "bact_prophage_ratio_num",
            "Prophage / contig length ratio",
            "B",
            log=False
        )


    else:
        draw(
            axes[0],
            "bact_genomad_length_num",
            "Prophage length (bp)",
            "A",
            log=True
        )

        draw(
            axes[1],
            "bact_prophage_ratio_num",
            "Prophage / contig length ratio",
            "B",
            log=False
        )

        draw(
            axes[2],
            "completeness",
            "CheckV completeness (%)",
            "C",
            log=False
        )

        axes[2].set_ylim(-3, 103)

    legend_elements = [
        Patch(facecolor=palette["integrated"], edgecolor="black", label="Integrated"),
        Patch(facecolor=palette["non_integrated"], edgecolor="black", label="Non-integrated"),
        Line2D(
            [0], [0],
            marker="o",
            color="black",
            label="Isolate",
            markerfacecolor="black",
            markeredgecolor="black",
            linestyle="None",
            markersize=6
        ),
        Line2D(
            [0], [0],
            marker="o",
            color="black",
            label="MAG",
            markerfacecolor="none",
            markeredgecolor="black",
            linestyle="None",
            markersize=6
        )
    ]

    if poster:
        axes[1].legend(
            handles=legend_elements,
            frameon=False,
            loc="lower left",
            bbox_to_anchor=(1.02, 0.0),
            borderaxespad=0.
        )
    else:
        axes[0].legend(
            handles=legend_elements,
            frameon=False,
            loc="best"
        )

    if poster:
        plt.tight_layout(rect=[0, 0, 0.86, 1])
    else:
        plt.tight_layout()

    png = outdir / "integration_boxplots.png"
    svg = outdir / "integration_boxplots.svg"

    plt.savefig(png, dpi=600, bbox_inches="tight")
    plt.savefig(svg, bbox_inches="tight")
    plt.close()

    return png, svg


# =========================
# Main
# =========================
def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.annotation_table, sep="\t")

    # Optionally exclude selected prophages
    if args.prophages_to_exclude:
        with open(args.prophages_to_exclude) as f:
            excluded = {line.strip() for line in f if line.strip()}

        n_before = len(df)
        df = df[~df["checkv_id"].isin(excluded)].copy()

        print(
            f"Excluded {n_before - len(df)} prophages "
            f"listed in {args.prophages_to_exclude}"
        )

    # Filter
    df = df[df["integration_status"].isin(["integrated", "non_integrated"])].copy()

    # Convert numeric
    cols = [
        "bact_genomad_length_num",
        "bact_prophage_ratio_num",
        "completeness"
    ]
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=cols)

    # =========================
    # Stats
    # =========================
    stats_results = []

    for col in cols:
        stats_results.append(compute_stats(df, col))

    stats_df = pd.DataFrame(stats_results)

    stats_tsv = outdir / "integration_stats.tsv"
    stats_txt = outdir / "integration_stats.txt"

    stats_df.to_csv(stats_tsv, sep="\t", index=False)

    # =========================
    # Source type composition
    # =========================
    composition_df = source_type_composition(df)
    composition_tsv = outdir / "source_type_composition_by_integration_status.tsv"
    composition_df.to_csv(composition_tsv, sep="\t", index=False)

    # =========================
    # Above-global-median summary
    # =========================
    above_median_df = above_median_summary(df, cols)
    above_median_tsv = outdir / "source_type_above_global_median_summary.tsv"
    above_median_df.to_csv(above_median_tsv, sep="\t", index=False)

    # =========================
    # Isolate vs MAG tests within integration groups
    # =========================
    source_tests_df = source_type_metric_tests(df, cols)
    source_tests_tsv = outdir / "source_type_metric_tests_within_integration_status.tsv"
    source_tests_df.to_csv(source_tests_tsv, sep="\t", index=False)

    print(f"Source type composition saved: {composition_tsv}")
    print(f"Above-median source type summary saved: {above_median_tsv}")
    print(f"Source type tests saved: {source_tests_tsv}")

    with open(stats_txt, "w") as f:
        for row in stats_results:
            f.write(f"{row['variable']}:\n")
            f.write(f"  n_integrated = {row['n_integrated']}\n")
            f.write(f"  n_non_integrated = {row['n_non_integrated']}\n")
            f.write(f"  median_integrated = {row['median_integrated']:.4g}\n")
            f.write(f"  median_non_integrated = {row['median_non_integrated']:.4g}\n")
            f.write(f"  p-value = {row['p_value']:.3e}\n\n")

    # =========================
    # Plot
    # =========================
    png, svg = make_plot(df, outdir, poster=args.poster)

    print(f"Figure saved: {png}")
    print(f"Figure saved: {svg}")
    print(f"Stats saved: {stats_tsv}")
    print(f"Stats saved: {stats_txt}")


if __name__ == "__main__":
    main()
