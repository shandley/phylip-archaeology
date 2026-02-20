#!/usr/bin/env python3
"""
Figure 2: Cross-disciplinary algorithmic connections found in PHYLIP.

Four-panel figure for Nature Methods manuscript:
  (a) Pruning algorithm = Belief propagation
  (b) Independent contrasts = Kirchhoff's circuit laws
  (c) Genetic code step matrix and optimization
  (d) LogDet compositional robustness

Usage:
    python manuscript/figures/plot_figure2.py
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, BoundaryNorm

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Genetic code data (used in panel c)
# ---------------------------------------------------------------------------
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_ORDER)}

POLAR_REQ = {
    "A": 7.0, "C": 4.8, "D": 13.0, "E": 12.5, "F": 5.0,
    "G": 7.9, "H": 8.4, "I": 4.9, "K": 10.1, "L": 4.9,
    "M": 5.3, "N": 10.0, "P": 6.6, "Q": 8.6, "R": 9.1,
    "S": 7.5, "T": 6.6, "V": 5.6, "W": 5.2, "Y": 5.4,
}

# Build codon table: index = pos1*16 + pos2*4 + pos3  (U=0 C=1 A=2 G=3)
CODON_TO_AA: list[str | None] = [None] * 64  # None = stop

_assignments: dict[str, list[tuple[int, int, int]]] = {
    "F": [(0,0,0),(0,0,1)],
    "L": [(0,0,2),(0,0,3),(1,0,0),(1,0,1),(1,0,2),(1,0,3)],
    "I": [(2,0,0),(2,0,1),(2,0,2)],
    "M": [(2,0,3)],
    "V": [(3,0,0),(3,0,1),(3,0,2),(3,0,3)],
    "S": [(0,1,0),(0,1,1),(0,1,2),(0,1,3),(2,3,0),(2,3,1)],
    "P": [(1,1,0),(1,1,1),(1,1,2),(1,1,3)],
    "T": [(2,1,0),(2,1,1),(2,1,2),(2,1,3)],
    "A": [(3,1,0),(3,1,1),(3,1,2),(3,1,3)],
    "Y": [(0,2,0),(0,2,1)],
    "H": [(1,2,0),(1,2,1)],
    "Q": [(1,2,2),(1,2,3)],
    "N": [(2,2,0),(2,2,1)],
    "K": [(2,2,2),(2,2,3)],
    "D": [(3,2,0),(3,2,1)],
    "E": [(3,2,2),(3,2,3)],
    "C": [(0,3,0),(0,3,1)],
    "W": [(0,3,3)],
    "R": [(1,3,0),(1,3,1),(1,3,2),(1,3,3),(2,3,2),(2,3,3)],
    "G": [(3,3,0),(3,3,1),(3,3,2),(3,3,3)],
    # Stops
    "*": [(0,2,2),(0,2,3),(0,3,2)],
}

for _aa, _codons in _assignments.items():
    for _p1, _p2, _p3 in _codons:
        _idx = _p1 * 16 + _p2 * 4 + _p3
        CODON_TO_AA[_idx] = _aa if _aa != "*" else None


def _codons_for_aa(aa: str) -> list[int]:
    """Return list of codon indices encoding amino acid *aa*."""
    return [i for i in range(64) if CODON_TO_AA[i] == aa]


def _hamming(c1: int, c2: int) -> int:
    """Hamming distance between two codons given as integer indices."""
    t1 = (c1 // 16, (c1 // 4) % 4, c1 % 4)
    t2 = (c2 // 16, (c2 // 4) % 4, c2 % 4)
    return sum(a != b for a, b in zip(t1, t2))


def build_step_matrix() -> np.ndarray:
    """20x20 minimum nucleotide substitution cost matrix."""
    n = len(AA_ORDER)
    mat = np.zeros((n, n), dtype=int)
    for i in range(n):
        codons_i = _codons_for_aa(AA_ORDER[i])
        for j in range(i + 1, n):
            codons_j = _codons_for_aa(AA_ORDER[j])
            best = 3
            for ci in codons_i:
                for cj in codons_j:
                    d = _hamming(ci, cj)
                    if d < best:
                        best = d
            mat[i, j] = best
            mat[j, i] = best
    return mat


def compute_ms_polar(codon_aa: list[str | None], pr_map: dict[str, float]) -> float:
    """
    Mean squared polar-requirement change for nonsynonymous point mutations.

    For every sense codon, for each of 3 positions, for each of 3 alternative
    nucleotides: if the mutation is nonsynonymous (different AA and not to a
    stop), accumulate (delta_polar)^2.  Return mean over all such mutations.
    """
    total = 0.0
    count = 0
    for idx in range(64):
        aa = codon_aa[idx]
        if aa is None:
            continue
        t = (idx // 16, (idx // 4) % 4, idx % 4)
        for pos in range(3):
            for nuc in range(4):
                if nuc == t[pos]:
                    continue
                new_t = list(t)
                new_t[pos] = nuc
                new_idx = new_t[0] * 16 + new_t[1] * 4 + new_t[2]
                new_aa = codon_aa[new_idx]
                if new_aa is None:
                    continue  # mutation to stop
                if new_aa == aa:
                    continue  # synonymous
                delta = pr_map[new_aa] - pr_map[aa]
                total += delta * delta
                count += 1
    return total / count if count > 0 else 0.0


def run_permutation_test(
    n_perm: int = 10_000, seed: int = 42
) -> tuple[np.ndarray, float]:
    """Permutation test: shuffle AA-to-PR mapping, return distribution & real value."""
    rng = np.random.default_rng(seed)

    # Real value
    real_val = compute_ms_polar(CODON_TO_AA, POLAR_REQ)

    # For permutations we shuffle which polar requirement value goes with which AA
    aa_list = sorted(POLAR_REQ.keys())
    pr_vals = np.array([POLAR_REQ[a] for a in aa_list])

    perm_vals = np.empty(n_perm)
    for k in range(n_perm):
        shuffled = rng.permutation(pr_vals)
        pr_map_perm = {aa_list[i]: float(shuffled[i]) for i in range(len(aa_list))}
        perm_vals[k] = compute_ms_polar(CODON_TO_AA, pr_map_perm)

    return perm_vals, real_val


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def _add_panel_label(ax, label: str):
    ax.text(
        -0.08, 1.08, label, transform=ax.transAxes,
        fontsize=10, fontweight="bold", va="top",
    )


# ===========================================================================
# Panel (a): Pruning = Belief Propagation
# ===========================================================================

def draw_panel_a(ax):
    ax.axis("off")
    _add_panel_label(ax, "a")
    ax.set_title("Pruning algorithm is belief propagation", pad=8)

    # --- Timeline (left side) ---
    tl_x = 0.10
    ax.plot([tl_x, tl_x], [0.15, 0.85], color="grey", linewidth=1.0,
            transform=ax.transAxes, clip_on=False)

    # Tick marks and labels
    for yr, ypos, col in [
        ("Felsenstein\n1981", 0.78, "#0072B2"),
        ("Pearl\n1988", 0.30, "#D55E00"),
    ]:
        ax.plot([tl_x - 0.02, tl_x + 0.02], [ypos, ypos],
                color=col, linewidth=1.5, transform=ax.transAxes, clip_on=False)
        ax.text(tl_x - 0.04, ypos, yr, ha="right", va="center",
                fontsize=6, color=col, fontweight="bold",
                transform=ax.transAxes)

    # Double-headed arrow with "7 years"
    ax.annotate(
        "", xy=(tl_x + 0.04, 0.30), xytext=(tl_x + 0.04, 0.78),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle="<->", color="grey", lw=0.8),
    )
    ax.text(tl_x + 0.07, 0.54, "7 years", ha="left", va="center",
            fontsize=6, color="grey", style="italic",
            transform=ax.transAxes, rotation=90)

    # --- Binary tree with message-passing arrows (right side) ---
    # Tree topology:  ((L0,L1),(L2,(L3,L4)))
    leaves_y = 0.12
    nodes = {
        "L0": (0.36, leaves_y),
        "L1": (0.48, leaves_y),
        "L2": (0.58, leaves_y),
        "L3": (0.72, leaves_y),
        "L4": (0.84, leaves_y),
        "N3": (0.78, 0.30),   # parent of L3,L4
        "N2": (0.65, 0.48),   # parent of L2, N3
        "N1": (0.42, 0.48),   # parent of L0, L1
        "R":  (0.54, 0.72),   # root
    }

    edges = [
        ("L0", "N1"), ("L1", "N1"),
        ("L2", "N2"), ("N3", "N2"),
        ("L3", "N3"), ("L4", "N3"),
        ("N1", "R"), ("N2", "R"),
    ]

    # Draw edges as L-shaped lines
    for child, parent in edges:
        cx, cy = nodes[child]
        px, py = nodes[parent]
        ax.plot([cx, cx], [cy, py], color="black", linewidth=0.7,
                transform=ax.transAxes, clip_on=False)
        ax.plot([cx, px], [py, py], color="black", linewidth=0.7,
                transform=ax.transAxes, clip_on=False)

    # Draw nodes
    for name, (x, y) in nodes.items():
        if name.startswith("L"):
            ax.plot(x, y, "o", color="#0072B2", markersize=4,
                    transform=ax.transAxes, clip_on=False, zorder=5)
        elif name == "R":
            ax.plot(x, y, "s", color="#D55E00", markersize=5,
                    transform=ax.transAxes, clip_on=False, zorder=5)
        else:
            ax.plot(x, y, "o", color="grey", markersize=3.5,
                    transform=ax.transAxes, clip_on=False, zorder=5)

    # Message-passing arrows (along branches, from child toward parent)
    arrow_color = "#009E73"
    for child, parent in edges:
        cx, cy = nodes[child]
        px, py = nodes[parent]
        mid_y = (cy + py) / 2
        dx = 0.013
        ax.annotate(
            "", xy=(cx + dx, mid_y + 0.06), xytext=(cx + dx, mid_y - 0.06),
            xycoords="axes fraction", textcoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", color=arrow_color, lw=0.7),
        )

    # Leaf labels
    for i, name in enumerate(["L0", "L1", "L2", "L3", "L4"]):
        x, y = nodes[name]
        ax.text(x, y - 0.06, f"$s_{i+1}$", ha="center", va="top",
                fontsize=6, color="#0072B2", transform=ax.transAxes)

    ax.text(nodes["R"][0], nodes["R"][1] + 0.04, "root", ha="center",
            va="bottom", fontsize=6, fontweight="bold", color="#D55E00",
            transform=ax.transAxes)

    # Recurrence formula
    ax.text(0.60, 0.93,
            r"$L(b) = \prod_c \sum_j P(b,j) \cdot L_c(j)$",
            ha="center", va="top", fontsize=7, transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey",
                      lw=0.5, alpha=0.9))


# ===========================================================================
# Panel (b): Independent Contrasts = Kirchhoff's Laws
# ===========================================================================

def draw_panel_b(ax):
    ax.axis("off")
    _add_panel_label(ax, "b")
    ax.set_title("Independent contrasts are Kirchhoff's laws", pad=8)

    # Tree: ((A:0.5,B:0.3):0.2,(C:0.4,(D:0.1,E:0.6):0.3):0.1)
    # Leaf values: A=2V, B=4V, C=1V, D=5V, E=3V
    # Use the upper-left quadrant of the panel for the tree,
    # the right side for the comparison table.
    nodes = {
        "A":   (0.03, 0.18),
        "B":   (0.13, 0.18),
        "AB":  (0.08, 0.38),
        "C":   (0.23, 0.18),
        "D":   (0.33, 0.18),
        "E":   (0.43, 0.18),
        "DE":  (0.38, 0.32),
        "CDE": (0.30, 0.48),
        "R":   (0.19, 0.62),
    }

    edges_info = [
        ("A", "AB", 0.5), ("B", "AB", 0.3),
        ("D", "DE", 0.1), ("E", "DE", 0.6),
        ("C", "CDE", 0.4), ("DE", "CDE", 0.3),
        ("AB", "R", 0.2), ("CDE", "R", 0.1),
    ]

    for child, parent, blen in edges_info:
        cx, cy = nodes[child]
        px, py = nodes[parent]
        ax.plot([cx, cx], [cy, py], color="black", linewidth=0.8,
                transform=ax.transAxes, clip_on=False)
        ax.plot([cx, px], [py, py], color="black", linewidth=0.8,
                transform=ax.transAxes, clip_on=False)
        # branch length label
        mid_y = (cy + py) / 2
        ax.text(cx + 0.018, mid_y, f"{blen}", ha="left", va="center",
                fontsize=5, color="grey", transform=ax.transAxes)

    # Leaf labels with voltages
    leaf_labels = {"A": "A=2V", "B": "B=4V", "C": "C=1V", "D": "D=5V", "E": "E=3V"}
    for name, label in leaf_labels.items():
        x, y = nodes[name]
        ax.plot(x, y, "o", color="#0072B2", markersize=3.5,
                transform=ax.transAxes, clip_on=False, zorder=5)
        ax.text(x, y - 0.045, label, ha="center", va="top",
                fontsize=5, fontweight="bold", transform=ax.transAxes)

    # Internal nodes
    for name in ["AB", "DE", "CDE", "R"]:
        x, y = nodes[name]
        ax.plot(x, y, "o", color="grey", markersize=3,
                transform=ax.transAxes, clip_on=False, zorder=5)

    ax.text(nodes["R"][0], nodes["R"][1] + 0.03, "root", ha="center",
            va="bottom", fontsize=5.5, color="#D55E00", fontweight="bold",
            transform=ax.transAxes)

    # --- Comparison table using matplotlib table ---
    # Place it on the right side of the panel.
    # Use plain text in a column layout, with monospace numbers right-aligned.
    # Use separate text calls per column to control alignment precisely.

    # Column positions (axes fraction x-coords for left edge of each column)
    c0_x = 0.52   # junction (left-aligned)
    c1_x = 0.72   # Felsenstein (left-aligned, monospace)
    c2_x = 0.88   # Kirchhoff (left-aligned, monospace)

    tbl_top = 0.65
    row_h = 0.075

    # Headers
    ax.text(c0_x, tbl_top, "Junction", ha="left", va="top",
            fontsize=5, fontweight="bold", transform=ax.transAxes)
    ax.text(c1_x, tbl_top, "Felsenstein", ha="center", va="top",
            fontsize=5, fontweight="bold", transform=ax.transAxes)
    ax.text(c2_x, tbl_top, "Kirchhoff", ha="center", va="top",
            fontsize=5, fontweight="bold", transform=ax.transAxes)
    # Header underline
    ax.plot([c0_x, 0.98], [tbl_top - 0.018, tbl_top - 0.018],
            color="grey", linewidth=0.4, transform=ax.transAxes, clip_on=False)

    # Data rows
    rows = [
        ("A-B",     "-2.236068", "-2.236068"),
        ("D-E",      "2.390457",  "2.390457"),
        ("C-DE",    "-4.190279", "-4.190279"),
        ("AB-CDE",   "0.434230",  "0.434230"),
    ]
    for r, (junc, f_val, k_val) in enumerate(rows):
        y = tbl_top - (r + 1) * row_h
        ax.text(c0_x, y, junc, ha="left", va="top",
                fontsize=4.5, family="monospace", transform=ax.transAxes)
        ax.text(c1_x, y, f_val, ha="center", va="top",
                fontsize=4.5, family="monospace", transform=ax.transAxes)
        ax.text(c2_x, y, k_val, ha="center", va="top",
                fontsize=4.5, family="monospace", transform=ax.transAxes)

    # "Match to 10 decimal places" below table
    ax.text(0.75, tbl_top - 5 * row_h + 0.015,
            "Match to 10 decimal places",
            ha="center", va="top", fontsize=5.5, style="italic",
            color="#009E73", transform=ax.transAxes)

    # Key formula at bottom of panel
    ax.text(0.50, 0.03,
            r"$v = \frac{v_L \cdot v_R}{v_L + v_R}$"
            r"$\;$= parallel resistor formula",
            ha="center", va="bottom", fontsize=6.5,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey",
                      lw=0.5, alpha=0.9))


# ===========================================================================
# Panel (c): Genetic Code Step Matrix & Optimization
# ===========================================================================

def draw_panel_c(fig, gs_cell):
    """Panel c occupies one gridspec cell but has two sub-panels."""
    inner = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs_cell, width_ratios=[1.1, 0.9], wspace=0.65,
    )
    ax_heat = fig.add_subplot(inner[0])
    ax_hist = fig.add_subplot(inner[1])

    # Place "c" label slightly higher and further left to avoid heatmap overlap
    ax_heat.text(-0.18, 1.08, "c", transform=ax_heat.transAxes,
                 fontsize=10, fontweight="bold", va="top")

    # --- Heatmap ---
    mat = build_step_matrix()
    cmap = ListedColormap(["white", "#a6cee3", "#1f78b4", "#08306b"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    im = ax_heat.imshow(mat, cmap=cmap, norm=norm, origin="upper", aspect="equal")

    ax_heat.set_xticks(range(20))
    ax_heat.set_yticks(range(20))
    ax_heat.set_xticklabels(AA_ORDER, fontsize=5)
    ax_heat.set_yticklabels(AA_ORDER, fontsize=5)
    ax_heat.tick_params(length=0, pad=2)
    ax_heat.set_title("Genetic code step matrix\nand optimization",
                       fontsize=8, pad=6)

    # Colorbar -- place below the heatmap to avoid overlapping the histogram
    cbar = fig.colorbar(im, ax=ax_heat, ticks=[0, 1, 2, 3],
                         orientation="horizontal",
                         fraction=0.06, pad=0.12, shrink=0.8)
    cbar.ax.set_xticklabels(["0", "1", "2", "3"], fontsize=5)
    cbar.ax.set_xlabel("Min. nucleotide substitutions", fontsize=5.5,
                         labelpad=2)
    cbar.outline.set_linewidth(0.4)

    # Restore spines for heatmap (overridden by rcParams)
    for sp in ax_heat.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.4)

    # --- Histogram ---
    perm_vals, real_val = run_permutation_test(n_perm=10_000, seed=42)
    perm_mean = np.mean(perm_vals)
    perm_std = np.std(perm_vals)
    z_score = (real_val - perm_mean) / perm_std

    ax_hist.hist(perm_vals, bins=50, color="#a6cee3", edgecolor="white",
                  linewidth=0.3, density=True)
    ax_hist.axvline(real_val, color="#e31a1c", linewidth=1.2, linestyle="--",
                     zorder=5)

    # Label for the real code line -- place in upper-left of histogram
    ax_hist.text(0.05, 0.95,
                  f"Real code\nz = {z_score:.2f}",
                  ha="left", va="top", fontsize=5.5, color="#e31a1c",
                  fontweight="bold", transform=ax_hist.transAxes)

    ax_hist.set_xlabel("Mean squared\npolar requirement change", fontsize=6)
    ax_hist.set_ylabel("Density", fontsize=6)
    ax_hist.set_title("Random code\npermutation test", fontsize=7, pad=4)

    return ax_heat, ax_hist


# ===========================================================================
# Panel (d): LogDet Compositional Robustness
# ===========================================================================

def draw_panel_d(ax):
    _add_panel_label(ax, "d")
    ax.set_title("LogDet resists compositional bias", pad=8)

    bias_rates = [20, 30, 40]
    jc69_correct = [84, 64, 74]
    logdet_correct = [100, 100, 80]

    x = np.arange(len(bias_rates))
    width = 0.32

    bars_jc = ax.bar(x - width / 2, jc69_correct, width,
                      color="#D55E00", label="JC69", edgecolor="white",
                      linewidth=0.4, zorder=3)
    bars_ld = ax.bar(x + width / 2, logdet_correct, width,
                      color="#0072B2", label="LogDet", edgecolor="white",
                      linewidth=0.4, zorder=3)

    ax.axhline(100, color="grey", linewidth=0.6, linestyle="--", zorder=2)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{r}%" for r in bias_rates])
    ax.set_xlabel("GC-bias rate")
    ax.set_ylabel("Correct topology (%)")
    ax.set_ylim(0, 115)
    ax.set_yticks([0, 20, 40, 60, 80, 100])

    # Value labels on bars
    for bar_set in [bars_jc, bars_ld]:
        for bar in bar_set:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 1.5,
                    f"{int(h)}%", ha="center", va="bottom", fontsize=5.5)

    ax.legend(loc="upper right", frameon=True, framealpha=0.9,
               edgecolor="grey", fancybox=False)


# ===========================================================================
# Main
# ===========================================================================

def main():
    fig = plt.figure(figsize=(180 / 25.4, 160 / 25.4))
    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        hspace=0.50, wspace=0.38,
        left=0.10, right=0.96, top=0.94, bottom=0.07,
    )

    # Panel (a) -- top-left
    ax_a = fig.add_subplot(gs[0, 0])
    draw_panel_a(ax_a)

    # Panel (b) -- top-right
    ax_b = fig.add_subplot(gs[0, 1])
    draw_panel_b(ax_b)

    # Panel (c) -- bottom-left (has sub-panels, pass gridspec cell)
    draw_panel_c(fig, gs[1, 0])

    # Panel (d) -- bottom-right
    ax_d = fig.add_subplot(gs[1, 1])
    draw_panel_d(ax_d)

    # --- Save ---
    out_dir = Path(__file__).resolve().parent
    for suffix in ("pdf", "png"):
        out_path = out_dir / f"figure2.{suffix}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
        print(f"Saved {out_path}")

    plt.close(fig)


if __name__ == "__main__":
    main()
