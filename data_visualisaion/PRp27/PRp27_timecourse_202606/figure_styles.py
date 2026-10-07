"""Shared plotting styles for the PRp27 analysis figures.

This module centralises visual choices that should stay consistent across the
individual analysis scripts and combined master figures. Add new plot-type
helpers here as the figure set grows.
"""

from matplotlib.collections import LineCollection, PathCollection
from matplotlib.lines import Line2D
import numpy as np

TEXT_COLOR = "#4A4A4A"
FONT_STACK = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]

# Reference-inspired palette: muted blue/orange/green first, then compatible
# extension colours for future datasets and secondary encodings.
PALETTE_SEQUENCE = [
    "#4E79A7",  # muted blue
    "#F28E2B",  # orange
    "#59A14F",  # green
    "#E15759",  # rose
    "#B07AA1",  # muted purple
    "#000101",  # teal
    "#EDC948",  # yellow
    "#FF9DA7",  # soft pink
    "#9C755F",  # brown
    "#BAB0AC",  # grey
]
COLORS = {"WT": PALETTE_SEQUENCE[0], "PRp27#1": PALETTE_SEQUENCE[1], "PRp27#2": PALETTE_SEQUENCE[2]}
MARKERS = {"WT": "o", "PRp27#1": "^", "PRp27#2": "s"}
LINESTYLES = {"WT": "-", "PRp27#1": "--", "PRp27#2": ":"}
MARKER_FACES = {line: color for line, color in COLORS.items()}
LINES = ["WT", "PRp27#1", "PRp27#2"]
DPI_COLORS = {0: PALETTE_SEQUENCE[0], 1: PALETTE_SEQUENCE[1], 3: PALETTE_SEQUENCE[2], 5: PALETTE_SEQUENCE[3]}

MASTER_RCPARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "axes.labelcolor": TEXT_COLOR,
    "axes.edgecolor": TEXT_COLOR,
    "axes.linewidth": 1.0,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "xtick.color": TEXT_COLOR,
    "ytick.color": TEXT_COLOR,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "legend.fontsize": 13,
    "legend.title_fontsize": 13,
    "text.color": TEXT_COLOR,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
}

ANALYSIS_RCPARAMS = {
    **MASTER_RCPARAMS,
    "font.size": 14,
    "axes.titlesize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "legend.title_fontsize": 12,
    "figure.dpi": 120,
    "savefig.dpi": 300,
}

INDIVIDUAL_POINT_ALPHA = 0.90
INDIVIDUAL_POINT_EDGE = "#000000"
INDIVIDUAL_POINT_LINEWIDTH = 0.45
INDIVIDUAL_POINT_STYLE = {
    "alpha": INDIVIDUAL_POINT_ALPHA,
    "edgecolor": INDIVIDUAL_POINT_EDGE,
    "linewidth": INDIVIDUAL_POINT_LINEWIDTH,
}
MEAN_MARKER_EDGE = TEXT_COLOR
MEAN_MARKER_LINEWIDTH = 0.8
ERRORBAR_COLOR = "#000000"
PANEL_LABEL_SIZE = 20
LETTER_FONT_SIZE = 14
LETTER_LABEL_GAP_FRACTION = 0.10
LETTER_LABEL_STEP_FRACTION = 0.14
LOG_LETTER_GAP = 1.70
LOG_LETTER_STEP = 1.90


def panel_label(ax, label, x=-0.18, y=1.12):
    """Place a bold panel label in a consistent figure position."""
    return ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=PANEL_LABEL_SIZE,
        fontweight="bold",
        color=TEXT_COLOR,
        fontfamily=FONT_STACK,
    )


def restyle_scatter_points(ax):
    """Make individual scatter points dark, opaque, and clearly outlined."""
    for collection in ax.collections:
        if not isinstance(collection, PathCollection):
            continue
        collection.set_edgecolor(INDIVIDUAL_POINT_EDGE)
        collection.set_linewidth(INDIVIDUAL_POINT_LINEWIDTH)
        collection.set_alpha(INDIVIDUAL_POINT_ALPHA)


def restyle_line_markers(ax):
    """Give only mean-line markers the outlined marker style."""
    for line in ax.lines:
        marker = line.get_marker()
        if marker in (None, "None", "", " "):
            continue
        line.set_markeredgecolor(MEAN_MARKER_EDGE)
        line.set_markeredgewidth(MEAN_MARKER_LINEWIDTH)


def make_errorbars_black(ax):
    """Set error-bar line collections/caps to black without touching points."""
    for collection in ax.collections:
        if isinstance(collection, LineCollection):
            collection.set_color(ERRORBAR_COLOR)

    for line in ax.lines:
        if line.get_linestyle() in ("None", "none", ""):
            line.set_color(ERRORBAR_COLOR)


def recolor_grouped_scatter_by_x(ax, colors=COLORS):
    """Re-apply filled group colours after helpers draw hollow grouped points."""
    for collection in ax.collections:
        if not isinstance(collection, PathCollection):
            continue
        offsets = collection.get_offsets()
        if len(offsets) == 0:
            continue
        x = float(np.ma.median(offsets[:, 0]))
        nearest = round(x)
        delta = x - nearest
        if delta < -0.09:
            collection.set_facecolor(colors["WT"])
        elif delta > 0.09:
            collection.set_facecolor(colors["PRp27#2"])
        else:
            collection.set_facecolor(colors["PRp27#1"])


def apply_axis_style(ax):
    """Apply shared master styling to axis text, ticks, spines, and markers."""
    ax.xaxis.label.set_color(TEXT_COLOR)
    ax.yaxis.label.set_color(TEXT_COLOR)
    ax.title.set_color(TEXT_COLOR)
    ax.title.set_fontfamily(FONT_STACK)
    ax.xaxis.label.set_fontfamily(FONT_STACK)
    ax.yaxis.label.set_fontfamily(FONT_STACK)
    ax.tick_params(colors=TEXT_COLOR)
    restyle_scatter_points(ax)
    restyle_line_markers(ax)

    for text in ax.texts:
        text.set_fontfamily(FONT_STACK)

    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)

    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_color(TEXT_COLOR)
        tick_label.set_fontfamily(FONT_STACK)


def legend_handles(lines=LINES, colors=COLORS, markers=MARKERS):
    """Return consistent sample-type legend handles."""
    display_labels = {
        "PRp27#1": "prp27-1",
        "PRp27#2": "prp27-2",
    }
    return [
        Line2D(
            [0],
            [0],
            color=colors[line],
            marker=markers[line],
            markerfacecolor=colors[line],
            markeredgecolor=MEAN_MARKER_EDGE,
            markeredgewidth=MEAN_MARKER_LINEWIDTH,
            linewidth=1.8,
            markersize=6,
            label=display_labels.get(line, line),
        )
        for line in lines
    ]
