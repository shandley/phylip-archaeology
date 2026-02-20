#!/usr/bin/env python3
"""
Figure 1: LLM-assisted code archaeology workflow.

Generates a workflow schematic showing 5 phases of LLM-assisted code
archaeology, connected by directional arrows with a feedback loop and
a human expertise annotation bar.

Produces both a full-page-width (180 mm) and single-column (88 mm)
version for Nature Methods formatting.

Usage:
    python manuscript/figures/plot_figure1.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path as MplPath

# ---------------------------------------------------------------------------
# Style (Nature Methods)
# ---------------------------------------------------------------------------
plt.rcParams.update(
    {
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
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

# ---------------------------------------------------------------------------
# Output paths (relative to this script)
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR


# ---------------------------------------------------------------------------
# Okabe-Ito inspired soft fills + darker borders
# ---------------------------------------------------------------------------
PHASE_COLORS = [
    {"fill": "#E8F0FE", "border": "#4285F4"},  # Phase 1 - blue
    {"fill": "#FFF3E0", "border": "#E8910C"},  # Phase 2 - orange
    {"fill": "#E8F5E9", "border": "#34A853"},  # Phase 3 - green
    {"fill": "#FDE7F0", "border": "#CC0066"},  # Phase 4 - pink/magenta
    {"fill": "#EDE7F6", "border": "#7B1FA2"},  # Phase 5 - purple
]

ARROW_COLOR = "#444444"
HUMAN_BAR_COLOR = "#F5F0E1"
HUMAN_BAR_BORDER = "#B8A77E"

# ---------------------------------------------------------------------------
# Phase data
# ---------------------------------------------------------------------------
PHASES = [
    {
        "num": 1,
        "title": "Source Code\nReading",
        "desc": "LLM reads and annotates\nlegacy C source code",
    },
    {
        "num": 2,
        "title": "Algorithm\nExtraction",
        "desc": "Algorithms expressed in\nmathematical notation",
    },
    {
        "num": 3,
        "title": "Literature\nCross-referencing",
        "desc": "Cross-referenced against\npublished work",
    },
    {
        "num": 4,
        "title": "Reimplementation",
        "desc": "Modern implementation\nwith validation tests",
    },
    {
        "num": 5,
        "title": "Cross-disciplinary\nConnections",
        "desc": "Connections to other\nfields identified",
    },
]


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _draw_box(ax, x, y, w, h, phase, colors, fontsize_title=7, fontsize_desc=5.5):
    """Draw a single phase box with number badge, title, and description."""
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012",
        facecolor=colors["fill"],
        edgecolor=colors["border"],
        linewidth=1.0,
        zorder=2,
    )
    ax.add_patch(box)

    # Phase number badge
    badge_r = 0.016
    badge_cx = x + 0.028
    badge_cy = y + h - 0.028
    badge = plt.Circle(
        (badge_cx, badge_cy),
        badge_r,
        facecolor=colors["border"],
        edgecolor="none",
        zorder=3,
    )
    ax.add_patch(badge)
    ax.text(
        badge_cx,
        badge_cy,
        str(phase["num"]),
        ha="center",
        va="center",
        fontsize=5.5,
        fontweight="bold",
        color="white",
        zorder=4,
    )

    # Title
    ax.text(
        x + w / 2,
        y + h * 0.63,
        phase["title"],
        ha="center",
        va="center",
        fontsize=fontsize_title,
        fontweight="bold",
        color="#222222",
        zorder=3,
        linespacing=1.15,
    )

    # Description
    ax.text(
        x + w / 2,
        y + h * 0.23,
        phase["desc"],
        ha="center",
        va="center",
        fontsize=fontsize_desc,
        color="#555555",
        zorder=3,
        linespacing=1.15,
    )


def _draw_straight_arrow(ax, x1, y1, x2, y2):
    """Draw a straight forward arrow between two points."""
    arrow = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>,head_width=4,head_length=3",
        connectionstyle="arc3,rad=0.0",
        color=ARROW_COLOR,
        linewidth=1.2,
        zorder=5,
        mutation_scale=1,
    )
    ax.add_patch(arrow)


def _draw_bezier_arrow(ax, points, color=ARROW_COLOR, linewidth=1.0,
                       linestyle="-", zorder=5):
    """Draw an arrow along a Bezier curve defined by control points.

    *points*: list of (x, y) tuples.
      - length 4 = single cubic Bezier
      - length 7 = two joined cubic Bezier segments
    """
    n = len(points)
    if n == 4:
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    elif n == 7:
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
            MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        ]
    else:
        raise ValueError(f"Expected 4 or 7 control points, got {n}")

    path = MplPath(points, codes)
    patch = mpatches.FancyArrowPatch(
        path=path,
        arrowstyle="-|>,head_width=4,head_length=3",
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        zorder=zorder,
        mutation_scale=1,
    )
    ax.add_patch(patch)


def _draw_human_bar(ax, x, y, w, h):
    """Draw the 'Human domain expertise at every stage' annotation bar."""
    bar = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.008",
        facecolor=HUMAN_BAR_COLOR,
        edgecolor=HUMAN_BAR_BORDER,
        linewidth=0.7,
        linestyle="--",
        zorder=1,
    )
    ax.add_patch(bar)
    ax.text(
        x + w / 2,
        y + h / 2,
        "Human domain expertise at every stage",
        ha="center",
        va="center",
        fontsize=6,
        fontstyle="italic",
        color="#6B5D3E",
        zorder=2,
    )


# ===================================================================
# WIDE VERSION  (full-page width, 180 mm x 70 mm)
# ===================================================================
def make_wide_figure():
    fig, ax = plt.subplots(figsize=(180 / 25.4, 70 / 25.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    n = len(PHASES)
    box_w = 0.145
    box_h = 0.48
    gap = 0.035
    total_w = n * box_w + (n - 1) * gap
    x_start = (1 - total_w) / 2
    y_box = 0.25

    # Human expertise bar
    _draw_human_bar(ax, x_start - 0.01, 0.08, total_w + 0.02, 0.12)

    # Draw boxes
    box_xs = []
    for i, phase in enumerate(PHASES):
        x = x_start + i * (box_w + gap)
        _draw_box(ax, x, y_box, box_w, box_h, phase, PHASE_COLORS[i])
        box_xs.append(x)

    # Forward arrows between consecutive boxes
    y_mid = y_box + box_h / 2
    for i in range(n - 1):
        _draw_straight_arrow(
            ax,
            box_xs[i] + box_w + 0.004,
            y_mid,
            box_xs[i + 1] - 0.004,
            y_mid,
        )

    # Feedback arrow: smooth arc above all boxes from phase-5 to phase-1
    x5_mid = box_xs[-1] + box_w / 2
    x1_mid = box_xs[0] + box_w / 2
    y_top = y_box + box_h + 0.015
    y_apex = 0.96

    _draw_bezier_arrow(
        ax,
        [
            (x5_mid, y_top),
            (x5_mid, y_apex),
            (x1_mid, y_apex),
            (x1_mid, y_top),
        ],
        color=PHASE_COLORS[4]["border"],
        linewidth=0.9,
        linestyle="--",
    )

    ax.text(
        (x5_mid + x1_mid) / 2,
        y_apex + 0.005,
        "Iterative refinement",
        ha="center",
        va="bottom",
        fontsize=5.5,
        fontstyle="italic",
        color=PHASE_COLORS[4]["border"],
    )

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    return fig


# ===================================================================
# SINGLE-COLUMN VERSION  (88 mm x 120 mm), 3-top / 2-bottom layout
# ===================================================================
def make_narrow_figure():
    fig, ax = plt.subplots(figsize=(88 / 25.4, 110 / 25.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box_w = 0.27
    box_h = 0.155
    h_gap = 0.05

    # Row geometry
    row1_total = 3 * box_w + 2 * h_gap
    row2_total = 2 * box_w + h_gap
    x1_start = (1 - row1_total) / 2
    x2_start = (1 - row2_total) / 2

    row1_y = 0.60
    row2_y = 0.33

    positions = []  # (x, cx, cy) per box

    # Row 1: phases 1-3
    for i in range(3):
        x = x1_start + i * (box_w + h_gap)
        _draw_box(
            ax, x, row1_y, box_w, box_h, PHASES[i], PHASE_COLORS[i],
            fontsize_title=6.5, fontsize_desc=5,
        )
        positions.append((x, x + box_w / 2, row1_y + box_h / 2))

    # Row 2: phases 4-5
    for i in range(2):
        x = x2_start + i * (box_w + h_gap)
        _draw_box(
            ax, x, row2_y, box_w, box_h, PHASES[3 + i], PHASE_COLORS[3 + i],
            fontsize_title=6.5, fontsize_desc=5,
        )
        positions.append((x, x + box_w / 2, row2_y + box_h / 2))

    # Arrows 1->2, 2->3 (horizontal in row 1)
    for i in range(2):
        _draw_straight_arrow(
            ax,
            positions[i][0] + box_w + 0.004,
            positions[i][2],
            positions[i + 1][0] - 0.004,
            positions[i][2],
        )

    # Arrow 3->4: route from bottom-center of box 3 to top-center of box 4
    # with a gentle S-curve that stays in the gap between rows
    x3_bot = positions[2][1]            # center-x of box 3
    y3_bot = row1_y - 0.004            # just below box 3
    x4_top = positions[3][1]            # center-x of box 4
    y4_top = row2_y + box_h + 0.004    # just above box 4
    y_mid_gap = (y3_bot + y4_top) / 2  # midpoint of gap

    _draw_bezier_arrow(
        ax,
        [
            (x3_bot, y3_bot),
            (x3_bot, y_mid_gap),
            (x4_top, y_mid_gap),
            (x4_top, y4_top),
        ],
        color=ARROW_COLOR,
        linewidth=1.2,
    )

    # Arrow 4->5 (horizontal in row 2)
    _draw_straight_arrow(
        ax,
        positions[3][0] + box_w + 0.004,
        positions[3][2],
        positions[4][0] - 0.004,
        positions[4][2],
    )

    # Feedback arrow: from bottom of box 5, route down-left-up to top of box 1
    x5_bot = positions[4][1]
    y5_bot = row2_y - 0.005
    x1_mid = positions[0][1]
    y1_top = row1_y + box_h + 0.005

    x_route = x1_start - 0.06    # how far left the feedback arrow goes
    y_route = row2_y - 0.06      # how far down the feedback arrow goes

    _draw_bezier_arrow(
        ax,
        [
            (x5_bot, y5_bot),          # start: bottom of box 5
            (x5_bot, y_route),         # control: straight down
            (x_route, y_route),        # control: across bottom
            (x_route, (y_route + y1_top) / 2),  # midpoint: left side
            (x_route, y1_top + 0.04),  # control: up along left
            (x1_mid, y1_top + 0.04),   # control: curve toward box 1
            (x1_mid, y1_top),          # end: top of box 1
        ],
        color=PHASE_COLORS[4]["border"],
        linewidth=0.9,
        linestyle="--",
    )

    ax.text(
        x_route - 0.02,
        (y_route + y1_top) / 2,
        "Iterative\nrefinement",
        ha="center",
        va="center",
        fontsize=5,
        fontstyle="italic",
        color=PHASE_COLORS[4]["border"],
        rotation=90,
    )

    # Human expertise bar
    _draw_human_bar(ax, x1_start - 0.02, 0.15, row1_total + 0.04, 0.09)

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    return fig


# ===================================================================
# Main
# ===================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Wide version
    fig_wide = make_wide_figure()
    fig_wide.savefig(OUT_DIR / "figure1.pdf", bbox_inches="tight", pad_inches=0.1)
    fig_wide.savefig(OUT_DIR / "figure1.png", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig_wide)
    print(f"Saved: {OUT_DIR / 'figure1.pdf'}")
    print(f"Saved: {OUT_DIR / 'figure1.png'}")

    # Single-column version
    fig_narrow = make_narrow_figure()
    fig_narrow.savefig(
        OUT_DIR / "figure1_single_column.pdf", bbox_inches="tight", pad_inches=0.1
    )
    fig_narrow.savefig(
        OUT_DIR / "figure1_single_column.png", bbox_inches="tight", pad_inches=0.1
    )
    plt.close(fig_narrow)
    print(f"Saved: {OUT_DIR / 'figure1_single_column.pdf'}")
    print(f"Saved: {OUT_DIR / 'figure1_single_column.png'}")


if __name__ == "__main__":
    main()
