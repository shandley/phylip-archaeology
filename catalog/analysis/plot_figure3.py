#!/usr/bin/env python3
"""Generate Nature Methods Figure 4: Software Catalog Preservation.

Three-panel figure:
  (a) Timeline of tool publications by 5-year bins, stacked by method category
  (b) Preservation status breakdown
  (c) Category evolution showing the parsimony → ML → Bayesian shift

Nature Methods specs: 180mm full-page width, 300+ DPI, Helvetica 7-8pt,
colorblind-safe palette (Okabe-Ito).
"""

import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

ENRICHED_FILE = Path(__file__).parent / "tools_enriched.json"
STATS_FILE = Path(__file__).parent / "stats" / "narrative_stats.json"
FIG_DIR = Path(__file__).parent / "figures"

# ── Okabe-Ito colorblind-safe palette ──────────────────────
COLORS = {
    # Status
    "archived": "#56B4E9",    # sky blue
    "dormant": "#E69F00",     # orange
    "dead": "#D55E00",        # vermillion
    "unknown": "#999999",     # gray

    # Categories
    "distance": "#0072B2",       # blue
    "parsimony": "#D55E00",      # vermillion
    "likelihood": "#009E73",     # green
    "bayesian": "#CC79A7",       # pink
    "bootstrap": "#56B4E9",      # sky blue
    "tree-visualization": "#E69F00",  # orange
    "alignment": "#F0E442",      # yellow

    # Languages
    "C": "#0072B2",
    "Java": "#D55E00",
    "Perl": "#009E73",
    "Python": "#E69F00",
    "R": "#CC79A7",
    "C++": "#56B4E9",
}


def setup_style():
    """Configure matplotlib for Nature Methods."""
    plt.rcParams.update({
        "font.family": "Helvetica",
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6,
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


def panel_a(ax, tools: list[dict]):
    """Timeline histogram: tools by 5-year bins, stacked by category."""
    # Define era bins (only include bins with data)
    era_range = range(1985, 2015, 5)  # Focus on 1985-2014
    eras = [f"{e}-{e+4}" for e in era_range]
    x = np.arange(len(eras))

    # Count tools by era and category
    top_cats = ["distance", "parsimony", "likelihood", "bayesian"]
    cat_labels = ["Distance", "Parsimony", "Max. Likelihood", "Bayesian"]

    # Build per-category counts per era
    era_cat_counts: dict[str, list[int]] = {c: [] for c in top_cats}
    era_other: list[int] = []

    for e in era_range:
        era_tools = [t for t in tools
                     if t.get("first_release") and (e <= t["first_release"] < e + 5)]
        for cat in top_cats:
            count = sum(1 for t in era_tools if cat in t.get("categories", []))
            era_cat_counts[cat].append(count)
        # "Other" = tools not in any top category
        other = sum(1 for t in era_tools
                    if not any(c in t.get("categories", []) for c in top_cats))
        era_other.append(other)

    # Stacked bar chart
    bar_width = 0.7
    bottom = np.zeros(len(eras))

    for cat, label in zip(top_cats, cat_labels):
        vals = np.array(era_cat_counts[cat])
        ax.bar(x, vals, bar_width, bottom=bottom, label=label,
               color=COLORS[cat], edgecolor="white", linewidth=0.3)
        bottom += vals

    ax.bar(x, era_other, bar_width, bottom=bottom, label="Other",
           color="#BBBBBB", edgecolor="white", linewidth=0.3)

    # Annotations for key milestones
    milestones = {
        1989: "PHYLIP 3.2",
        2001: "MrBayes",
        2004: "RAxML",
    }
    for year, label in milestones.items():
        idx = (year - 1985) / 5
        if 0 <= idx < len(eras):
            ax.annotate(label, xy=(idx, bottom[int(idx)] + 1),
                       fontsize=5, ha="center", va="bottom",
                       color="#444444", style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels([f"{e}\n–{e+4}" for e in era_range], fontsize=5.5)
    ax.set_ylabel("Number of tools")
    ax.set_title("Timeline of phylogenetics software development", fontweight="bold", pad=8)
    ax.legend(loc="upper left", frameon=False, ncol=2, fontsize=5.5)

    # Add total count annotation
    total_with_year = sum(1 for t in tools if t.get("first_release"))
    ax.text(0.98, 0.95, f"n = {total_with_year} tools\nwith known year",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5,
            color="#666666")


def panel_b(ax, tools: list[dict]):
    """Preservation status breakdown — horizontal stacked bar."""
    status_counts = Counter(t["status"] for t in tools)
    total = len(tools)

    categories = ["archived", "dormant", "dead", "unknown"]
    labels = ["Archived", "Dormant", "Dead", "Unknown"]
    counts = [status_counts.get(c, 0) for c in categories]
    pcts = [c / total * 100 for c in counts]

    # Single horizontal stacked bar
    left = 0
    for cat, label, count, pct in zip(categories, labels, counts, pcts):
        width = count
        bar = ax.barh(0, width, left=left, height=0.5,
                      color=COLORS[cat], edgecolor="white", linewidth=0.5)
        # Label inside bar if wide enough
        if pct > 8:
            ax.text(left + width / 2, 0, f"{label}\n{count} ({pct:.0f}%)",
                    ha="center", va="center", fontsize=6, fontweight="bold",
                    color="white" if cat in ("dead",) else "black")
        left += width

    ax.set_xlim(0, total)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel(f"Tools (n = {total})")
    ax.set_title("Preservation status of cataloged tools", fontweight="bold", pad=8)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(ticker.MultipleLocator(50))

    # Add context below
    with_archive = sum(1 for t in tools if t.get("url_archive"))
    dead_no_archive = sum(1 for t in tools
                          if t["status"] == "dead" and not t.get("url_archive"))
    ax.text(0.5, -0.6, f"Wayback Machine archived: {with_archive} tools  |  "
            f"Dead with no archive: {dead_no_archive} tools",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.5,
            color="#666666")


def panel_c(ax, tools: list[dict]):
    """Category evolution: stacked area chart showing method paradigm shifts."""
    era_range = range(1985, 2015, 5)
    eras = list(era_range)

    cats = ["parsimony", "distance", "likelihood", "bayesian"]
    cat_labels = ["Parsimony", "Distance", "Max. Likelihood", "Bayesian"]

    # Compute proportions per era
    data = {c: [] for c in cats}
    for e in era_range:
        era_tools = [t for t in tools
                     if t.get("first_release") and (e <= t["first_release"] < e + 5)]
        era_total = len(era_tools) if era_tools else 1
        for cat in cats:
            count = sum(1 for t in era_tools if cat in t.get("categories", []))
            data[cat].append(count / era_total * 100)

    # Stacked area
    x = np.arange(len(eras))
    bottom = np.zeros(len(eras))

    for cat, label in zip(cats, cat_labels):
        vals = np.array(data[cat])
        ax.fill_between(x, bottom, bottom + vals, alpha=0.8,
                        color=COLORS[cat], label=label, linewidth=0.5,
                        edgecolor="white")
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([f"{e}\n–{e+4}" for e in era_range], fontsize=5.5)
    ax.set_ylabel("% of tools in era")
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(25))
    ax.set_title("Methodological paradigm shift", fontweight="bold", pad=8)
    ax.legend(loc="upper left", frameon=False, fontsize=5.5)

    # Annotation
    ax.annotate("Bayesian\nrevolution", xy=(3, 70), fontsize=5,
                ha="center", color="#444444", style="italic")


def main():
    setup_style()
    tools = json.loads(ENRICHED_FILE.read_text())

    # Create 3-panel figure (vertical layout, single column for Nature Methods)
    fig, axes = plt.subplots(3, 1, figsize=(88 / 25.4, 200 / 25.4),
                             gridspec_kw={"hspace": 0.55})

    panel_a(axes[0], tools)
    panel_b(axes[1], tools)
    panel_c(axes[2], tools)

    # Panel labels
    for ax, label in zip(axes, ["a", "b", "c"]):
        ax.text(-0.08, 1.08, label, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Save in multiple formats
    fig.savefig(FIG_DIR / "figure3_combined.pdf", bbox_inches="tight",
                pad_inches=0.1)
    fig.savefig(FIG_DIR / "figure3_combined.png", bbox_inches="tight",
                pad_inches=0.1, dpi=300)
    plt.close()

    print(f"Figure 3 saved:")
    print(f"  {FIG_DIR / 'figure3_combined.pdf'}")
    print(f"  {FIG_DIR / 'figure3_combined.png'}")

    # Also generate full-width version
    fig2, axes2 = plt.subplots(1, 3, figsize=(180 / 25.4, 55 / 25.4),
                               gridspec_kw={"wspace": 0.35})

    panel_a(axes2[0], tools)
    panel_b(axes2[1], tools)
    panel_c(axes2[2], tools)

    for ax, label in zip(axes2, ["a", "b", "c"]):
        ax.text(-0.12, 1.12, label, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")

    fig2.savefig(FIG_DIR / "figure3_wide.pdf", bbox_inches="tight",
                 pad_inches=0.1)
    fig2.savefig(FIG_DIR / "figure3_wide.png", bbox_inches="tight",
                 pad_inches=0.1, dpi=300)
    plt.close()

    print(f"  {FIG_DIR / 'figure3_wide.pdf'}")
    print(f"  {FIG_DIR / 'figure3_wide.png'}")


if __name__ == "__main__":
    main()
