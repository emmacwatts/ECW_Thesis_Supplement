"""Draw an abstract two-route schematic for cloning acdS into pBBR1MCS-3."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch, Wedge


HERE = Path(__file__).resolve().parent
OUT = HERE / "202608_acds_cloning_strategy_figures"
LOGS = HERE / "202608_acds_cloning_strategy_logs"

FONT_STACK = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]
TEXT = "#4A4A4A"
PANEL_BG = "#F3F6F9"
BACKBONE = "#4E79A7"
ACDS = "#59A14F"
PROMOTER_NATIVE = "#B07AA1"
PROMOTER_VIRE = "#7A5195"
TERMINATOR = "#6B6B6B"
TETR = "#F28E2B"
ARROW = "#555555"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 11,
    "axes.labelcolor": TEXT,
    "text.color": TEXT,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})


def rounded_box(ax, xy, width, height, *, facecolor="white", edgecolor="#B8C2CC",
                linewidth=1.0, radius=0.025, zorder=1):
    box = FancyBboxPatch(
        xy, width, height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, zorder=zorder,
    )
    ax.add_patch(box)
    return box


def arrow(ax, start, end, *, color=ARROW, mutation_scale=13, linewidth=1.3):
    patch = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=mutation_scale,
                            linewidth=linewidth, color=color, shrinkA=2, shrinkB=2)
    ax.add_patch(patch)
    return patch


def cassette(ax, x, y, width, *, promoter, promoter_label, terminator=False,
             no_terminator_label=False):
    """Draw promoter → acdS → optional terminator as a compact linear cassette."""
    promoter_width = width * 0.30
    gene_width = width * 0.43
    term_width = width * 0.17
    height = 0.075
    rounded_box(ax, (x, y), promoter_width, height, facecolor=promoter,
                edgecolor="white", linewidth=0.7, radius=0.012)
    ax.text(x + promoter_width / 2, y + height / 2, promoter_label,
            ha="center", va="center", color="white", fontsize=9.5, fontweight="bold")
    arrow(ax, (x + promoter_width + 0.008, y + height / 2),
          (x + promoter_width + 0.055, y + height / 2), color=ACDS, mutation_scale=10)
    gene_x = x + promoter_width + 0.06
    rounded_box(ax, (gene_x, y), gene_width, height, facecolor=ACDS,
                edgecolor="white", linewidth=0.7, radius=0.012)
    ax.text(gene_x + gene_width / 2, y + height / 2, r"$\it{acdS}$",
            ha="center", va="center", color="white", fontsize=10.5, fontweight="bold")
    right = gene_x + gene_width
    if terminator:
        arrow(ax, (right + 0.008, y + height / 2),
              (right + 0.046, y + height / 2), color=TERMINATOR, mutation_scale=9)
        term_x = right + 0.05
        rounded_box(ax, (term_x, y), term_width, height, facecolor=TERMINATOR,
                    edgecolor="white", linewidth=0.7, radius=0.012)
        ax.text(term_x + term_width / 2, y + height / 2, "rrnB\nT1/T2",
                ha="center", va="center", color="white", fontsize=8.5, fontweight="bold",
                linespacing=0.85)
    elif no_terminator_label:
        ax.text(right + 0.02, y + height / 2, "no added\nterminator",
                ha="left", va="center", fontsize=8.5, color="#777777", linespacing=0.9)


def plasmid(ax, center, radius, cassette_kind):
    """Draw a simple final plasmid with backbone, TetR, and expression cassette."""
    cx, cy = center
    ax.add_patch(Arc(center, 2 * radius, 2 * radius, theta1=0, theta2=360,
                     linewidth=8, color=BACKBONE, capstyle="round"))
    ax.add_patch(Wedge(center, radius + 0.015, 25, 77, width=0.065,
                       facecolor=TETR, edgecolor="white", linewidth=0.8))
    ax.text(cx + radius + 0.035, cy + radius * 0.48, "TetR", ha="left", va="center",
            fontsize=9, fontweight="bold", color=TETR)
    ax.text(cx, cy + 0.015, "pBBR1MCS-3", ha="center", va="center",
            fontsize=10.5, fontweight="bold", color=TEXT)
    ax.text(cx, cy - 0.040, "broad-host-range backbone", ha="center", va="center",
            fontsize=8.3, color="#6D7882")
    # A coloured cassette segment at the bottom communicates the construct
    # without implying exact feature lengths or orientation.
    ax.add_patch(Wedge(center, radius + 0.015, 215, 322, width=0.065,
                       facecolor=ACDS, edgecolor="white", linewidth=0.8))
    if cassette_kind == "native":
        label = r"native T3 | $\it{acdS}$"
    else:
        label = r"PvirE | $\it{acdS}$ | rrnB T1/T2"
    ax.text(cx, cy - radius - 0.055, label, ha="center", va="top",
            fontsize=9.2, color=TEXT)


def scissors(ax, x, y, size=0.035):
    """Small vector scissors icon without relying on a font glyph."""
    ax.add_patch(Circle((x - size * 0.35, y + size * 0.28), size * 0.20,
                        fill=False, edgecolor=TEXT, linewidth=1.2))
    ax.add_patch(Circle((x - size * 0.35, y - size * 0.28), size * 0.20,
                        fill=False, edgecolor=TEXT, linewidth=1.2))
    ax.plot([x - size * 0.18, x + size * 0.48],
            [y + size * 0.16, y - size * 0.38], color=TEXT, linewidth=1.2)
    ax.plot([x - size * 0.18, x + size * 0.48],
            [y - size * 0.16, y + size * 0.38], color=TEXT, linewidth=1.2)


def promoter_symbol(ax, x, y, color, label):
    """Standard bent-arrow promoter symbol."""
    ax.plot([x, x, x + 0.052], [y - 0.035, y + 0.018, y + 0.018],
            color=color, linewidth=2.2, solid_capstyle="round")
    arrow(ax, (x + 0.050, y + 0.018), (x + 0.085, y + 0.018),
          color=color, mutation_scale=10, linewidth=2.2)
    ax.text(x + 0.025, y - 0.060, label, ha="center", va="top",
            fontsize=8.3, color=color, fontweight="bold")


def terminator_symbol(ax, x, y, label="rrnB T1/T2"):
    """Standard T-shaped transcription terminator symbol."""
    ax.plot([x, x], [y - 0.035, y + 0.035], color=TERMINATOR, linewidth=2.2)
    ax.plot([x - 0.025, x + 0.025], [y + 0.035, y + 0.035],
            color=TERMINATOR, linewidth=2.2)
    ax.text(x, y - 0.058, label, ha="center", va="top", fontsize=8.0,
            color=TERMINATOR, fontweight="bold")


def simple_plasmid(ax, center, radius, *, final=None):
    """Minimal plasmid map with optional expression cassette symbols."""
    cx, cy = center
    ax.add_patch(Arc(center, 2 * radius, 2 * radius, theta1=0, theta2=360,
                     linewidth=9, color=BACKBONE, capstyle="round"))
    ax.add_patch(Wedge(center, radius + 0.017, 25, 78, width=0.075,
                       facecolor=TETR, edgecolor="white", linewidth=0.8))
    ax.text(cx + radius + 0.020, cy + radius * 0.48, "TetR", ha="left", va="center",
            fontsize=8.8, color=TETR, fontweight="bold")
    ax.text(cx, cy + 0.012, "pBBR1MCS-3", ha="center", va="center",
            fontsize=9.5, color=TEXT, fontweight="bold")
    if final is None:
        # MCS marker at the lower edge.
        ax.plot([cx - 0.018, cx + 0.018], [cy - radius, cy - radius],
                color="#222222", linewidth=3)
        return

    # acdS is shown as a single green feature; orientation is intentionally
    # abstracted because it is not part of the biological comparison.
    ax.add_patch(Wedge(center, radius + 0.017, 218, 310, width=0.075,
                       facecolor=ACDS, edgecolor="white", linewidth=0.8))
    ax.text(cx, cy - radius * 0.64, r"$\it{acdS}$", ha="center", va="center",
            fontsize=10, color="white", fontweight="bold")
    if final == "native":
        promoter_symbol(ax, cx - radius * 0.88, cy - radius * 0.48,
                        PROMOTER_NATIVE, "native T3")
    else:
        promoter_symbol(ax, cx - radius * 0.90, cy - radius * 0.48,
                        PROMOTER_VIRE, "PvirE")
        terminator_symbol(ax, cx + radius * 0.87, cy - radius * 0.47)


def inverse_pcr_primers(ax, center, radius):
    """Two inward-facing primer arrows flanking the region removed by inverse PCR."""
    cx, cy = center
    y = cy - radius * 0.86
    arrow(ax, (cx - radius * 0.95, y - 0.020), (cx - radius * 0.28, y + 0.005),
          color=PROMOTER_VIRE, mutation_scale=10, linewidth=1.8)
    arrow(ax, (cx + radius * 0.95, y + 0.020), (cx + radius * 0.28, y - 0.005),
          color=PROMOTER_VIRE, mutation_scale=10, linewidth=1.8)
    ax.text(cx, cy - radius - 0.060, "inverse-PCR primers", ha="center", va="top",
            fontsize=8.2, color=PROMOTER_VIRE, fontweight="bold")


def route_panel(ax, route):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    rounded_box(ax, (0.015, 0.06), 0.97, 0.86, facecolor=PANEL_BG,
                edgecolor="#D5DDE4", linewidth=0.9, radius=0.025, zorder=0)
    start = (0.27, 0.50)
    result = (0.74, 0.50)
    radius = 0.128
    simple_plasmid(ax, start, radius)
    simple_plasmid(ax, result, radius,
                   final="native" if route == "native" else "engineered")

    if route == "native":
        scissors(ax, start[0], start[1] - radius, size=0.050)
        ax.text(start[0], 0.285, "SmaI digest at MCS", ha="center", va="top",
                fontsize=9.2, fontweight="bold")
        process = "Gibson assembly"
    else:
        inverse_pcr_primers(ax, start, radius)
        process = "inverse PCR + Gibson"

    arrow(ax, (0.415, 0.50), (0.595, 0.50), mutation_scale=16, linewidth=1.7)
    ax.text(0.505, 0.655, process, ha="center", va="bottom",
            fontsize=9.2, fontweight="bold", color=TEXT)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 6.1), facecolor="white")
    route_panel(axes[0], "native")
    route_panel(axes[1], "engineered")
    fig.text(0.015, 0.955, "A", fontsize=16, fontweight="bold", ha="left", va="top")
    fig.text(0.505, 0.955, "B", fontsize=16, fontweight="bold", ha="left", va="top")
    fig.subplots_adjust(left=0.035, right=0.985, top=0.95, bottom=0.05, wspace=0.08)

    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"202608_acds_pBBR1MCS3_cloning_strategy.{extension}",
                    dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    LOGS.mkdir(parents=True, exist_ok=True)
    note = """acdS cloning-strategy schematic

Panel A represents the simple route described by the user: SmaI linearisation of pBBR1MCS-3 followed by Gibson assembly of acdS under the native T3 promoter, with no added terminator.

Panel B represents the engineered-expression route described by the user and shown in the supplied plasmid-map image: inverse-PCR shortening of the backbone/removal of the native promoter region, followed by Gibson assembly of an approximately 700-bp virE upstream/promoter region, acdS, and a strong bacterial rrnB T1/T2 terminator.

Both constructs retain tetracycline resistance (TetR). The diagram is intentionally abstract and not to scale. Insert direction is not used as a distinguishing feature. The source notebook was not visible as a separate workspace file at generation time; promoter/terminator wording therefore follows the user's description and should be checked against the final sequence record before publication.
"""
    (LOGS / "acds_cloning_strategy_note.txt").write_text(note, encoding="utf-8")


if __name__ == "__main__":
    main()
