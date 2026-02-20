#!/usr/bin/env python3
"""Generate Nature Methods Figure 5: Benchmarking phylip-rs vs modern tools.

Three-panel figure:
  (a) Log-likelihood vs wall time (accuracy-speed tradeoff)
  (b) Normalized RF distance to true tree by dataset size
  (c) Scaling behavior: wall time vs number of taxa (log-log)
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

RESULTS_FILE = Path(__file__).parent / "results" / "benchmark_results.csv"
FIG_DIR = Path(__file__).parent / "figures"

# Okabe-Ito colorblind-safe palette
TOOL_COLORS = {
    "phylip-rs-ml": "#0072B2",      # blue
    "iqtree": "#E69F00",            # orange
    "raxml-ng": "#009E73",          # green
    "veryfasttree": "#CC79A7",      # pink
    "phylip-rs-nj": "#56B4E9",      # sky blue
}

TOOL_LABELS = {
    "phylip-rs-ml": "phylip-rs ML",
    "iqtree": "IQ-TREE 3",
    "raxml-ng": "RAxML-NG",
    "veryfasttree": "VeryFastTree",
    "phylip-rs-nj": "phylip-rs NJ",
}

TOOL_MARKERS = {
    "phylip-rs-ml": "o",
    "iqtree": "s",
    "raxml-ng": "^",
    "veryfasttree": "D",
    "phylip-rs-nj": "v",
}

ML_TOOLS = ["phylip-rs-ml", "iqtree", "raxml-ng", "veryfasttree"]
ALL_TOOLS = ML_TOOLS + ["phylip-rs-nj"]


def setup_style():
    plt.rcParams.update({
        "font.family": "Helvetica",
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 5.5,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def panel_a(ax, df: pd.DataFrame):
    """Log-likelihood vs wall time scatter."""
    ml_df = df[df["tool"].isin(ML_TOOLS) & ~df["timed_out"] & df["lnl_jc69"].notna()]

    for tool in ML_TOOLS:
        subset = ml_df[ml_df["tool"] == tool]
        if subset.empty:
            continue
        ax.scatter(subset["wall_time_s"], subset["lnl_jc69"],
                   c=TOOL_COLORS[tool], marker=TOOL_MARKERS[tool],
                   s=15, alpha=0.7, label=TOOL_LABELS[tool],
                   edgecolors="none", zorder=3)

    ax.set_xscale("log")
    ax.set_xlabel("Wall time (s)")
    ax.set_ylabel("Log-likelihood (JC69)")
    ax.set_title("Accuracy vs. speed", fontweight="bold", pad=6)
    ax.legend(frameon=False, loc="lower right", markerscale=1.2)


def panel_b(ax, df: pd.DataFrame):
    """Normalized RF distance to true tree by dataset size."""
    ml_df = df[df["tool"].isin(ML_TOOLS) & ~df["timed_out"] & df["nrf_distance"].notna()]

    taxa_sizes = sorted(ml_df["ntaxa"].unique())
    n_tools = len(ML_TOOLS)
    bar_width = 0.8 / n_tools
    x = np.arange(len(taxa_sizes))

    for i, tool in enumerate(ML_TOOLS):
        subset = ml_df[ml_df["tool"] == tool]
        means = []
        stds = []
        for nt in taxa_sizes:
            vals = subset[subset["ntaxa"] == nt]["nrf_distance"]
            means.append(vals.mean() if len(vals) > 0 else np.nan)
            stds.append(vals.std() if len(vals) > 1 else 0)

        offset = (i - n_tools / 2 + 0.5) * bar_width
        ax.bar(x + offset, means, bar_width, yerr=stds,
               color=TOOL_COLORS[tool], label=TOOL_LABELS[tool],
               edgecolor="white", linewidth=0.3, capsize=1.5,
               error_kw={"linewidth": 0.5})

    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in taxa_sizes])
    ax.set_xlabel("Number of taxa")
    ax.set_ylabel("Normalized RF distance")
    ax.set_title("Topology accuracy", fontweight="bold", pad=6)
    ax.legend(frameon=False, loc="upper left", fontsize=5)
    ax.set_ylim(bottom=0)


def panel_c(ax, df: pd.DataFrame):
    """Scaling: wall time vs number of taxa (log-log), fixed at 1000 sites."""
    scaling_df = df[df["nsites"] == 1000]
    # Also include NJ for comparison
    tools_to_plot = ALL_TOOLS

    for tool in tools_to_plot:
        subset = scaling_df[(scaling_df["tool"] == tool) & ~scaling_df["timed_out"]]
        if subset.empty:
            continue
        medians = subset.groupby("ntaxa")["wall_time_s"].median()
        ax.plot(medians.index, medians.values,
                color=TOOL_COLORS[tool], marker=TOOL_MARKERS[tool],
                markersize=3.5, linewidth=1, label=TOOL_LABELS[tool])

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of taxa")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Scaling behavior (1000 sites)", fontweight="bold", pad=6)
    ax.legend(frameon=False, loc="upper left", fontsize=5)

    # Reference lines for O(n^2) and O(n^3)
    taxa_range = np.array([10, 500])
    # Faint reference lines
    ax.plot(taxa_range, 0.001 * (taxa_range / 10) ** 2,
            "--", color="#CCCCCC", linewidth=0.5, zorder=1)
    ax.plot(taxa_range, 0.0001 * (taxa_range / 10) ** 3,
            ":", color="#CCCCCC", linewidth=0.5, zorder=1)
    ax.text(400, 0.001 * (400 / 10) ** 2 * 1.5, "O(n²)",
            fontsize=5, color="#AAAAAA")
    ax.text(400, 0.0001 * (400 / 10) ** 3 * 1.5, "O(n³)",
            fontsize=5, color="#AAAAAA")


def main():
    setup_style()

    if not RESULTS_FILE.exists():
        print(f"Error: {RESULTS_FILE} not found. Run run_benchmarks.py first.")
        return

    df = pd.read_csv(RESULTS_FILE)
    # Ensure boolean column
    if "timed_out" in df.columns:
        df["timed_out"] = df["timed_out"].fillna(False).astype(bool)
    else:
        df["timed_out"] = False

    print(f"Loaded {len(df)} benchmark results")
    print(f"Tools: {sorted(df['tool'].unique())}")
    print(f"Dataset sizes: {sorted(df['ntaxa'].unique())} taxa")

    # Vertical layout (single column, Nature Methods)
    fig, axes = plt.subplots(3, 1, figsize=(88 / 25.4, 200 / 25.4),
                             gridspec_kw={"hspace": 0.45})

    panel_a(axes[0], df)
    panel_b(axes[1], df)
    panel_c(axes[2], df)

    for ax, label in zip(axes, ["a", "b", "c"]):
        ax.text(-0.08, 1.08, label, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "figure5_combined.pdf", bbox_inches="tight", pad_inches=0.1)
    fig.savefig(FIG_DIR / "figure5_combined.png", bbox_inches="tight",
                pad_inches=0.1, dpi=300)
    plt.close()

    # Also wide version
    fig2, axes2 = plt.subplots(1, 3, figsize=(180 / 25.4, 55 / 25.4),
                               gridspec_kw={"wspace": 0.35})
    panel_a(axes2[0], df)
    panel_b(axes2[1], df)
    panel_c(axes2[2], df)

    for ax, label in zip(axes2, ["a", "b", "c"]):
        ax.text(-0.12, 1.12, label, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")

    fig2.savefig(FIG_DIR / "figure5_wide.pdf", bbox_inches="tight", pad_inches=0.1)
    fig2.savefig(FIG_DIR / "figure5_wide.png", bbox_inches="tight",
                 pad_inches=0.1, dpi=300)
    plt.close()

    print(f"\nFigures saved to {FIG_DIR}/")

    # Print summary statistics
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS FOR PAPER")
    print(f"{'='*60}")

    for tool in ALL_TOOLS:
        subset = df[(df["tool"] == tool) & ~df["timed_out"]]
        if subset.empty:
            continue
        print(f"\n{TOOL_LABELS[tool]}:")
        print(f"  Datasets completed: {len(subset)}")
        print(f"  Median wall time: {subset['wall_time_s'].median():.3f}s")
        if "lnl_jc69" in subset.columns and subset["lnl_jc69"].notna().any():
            print(f"  Median lnL: {subset['lnl_jc69'].median():.1f}")
        if "nrf_distance" in subset.columns and subset["nrf_distance"].notna().any():
            print(f"  Median nRF: {subset['nrf_distance'].median():.3f}")
        if "peak_mem_kb" in subset.columns and subset["peak_mem_kb"].notna().any():
            print(f"  Median peak mem: {subset['peak_mem_kb'].median():.0f} KB")


if __name__ == "__main__":
    main()
