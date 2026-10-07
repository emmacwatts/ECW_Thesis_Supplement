"""Rebuild the four ethylene thesis figures in the shared house style.
"""

from __future__ import annotations

import csv
import io
import math
import tempfile
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from matplotlib import font_manager
from matplotlib.gridspec import GridSpec
from PIL import Image
from scipy import stats
from statsmodels.formula.api import ols


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "originalData"
ORIGINAL = DATA / "OriginalFigures"
OUT = ROOT / "202608_ethylene_figures"
LOGS = ROOT / "202608_ethylene_logs"

PANEL_BG = "#E9EFF6"
GRID = "#FFFFFF"
TEXT = "#303030"
WT = "#4E79A7"
MUTANT = "#F28E2B"
EDGE = "#303030"
FONT = "Helvetica Neue"
BAR_ALPHA = 0.76
POINT_ALPHA = 0.90


def register_helvetica_neue_light() -> None:
    """Expose Light TTC faces to Matplotlib, which otherwise selects Regular."""
    collection_path = Path("/System/Library/Fonts/HelveticaNeue.ttc")
    if not collection_path.exists():
        return
    try:
        from fontTools.ttLib import TTCollection

        collection = TTCollection(collection_path)
        cache = Path(tempfile.gettempdir()) / "ethylene_helvetica_neue"
        cache.mkdir(exist_ok=True)
        for index, filename in (
            (7, "HelveticaNeue-Light.ttf"),
            (8, "HelveticaNeue-LightItalic.ttf"),
            (1, "HelveticaNeue-Bold.ttf"),
            (3, "HelveticaNeue-BoldItalic.ttf"),
        ):
            target = cache / filename
            if not target.exists():
                collection.fonts[index].save(target)
            font_manager.fontManager.addfont(target)
        # The unsplit TTC entries point Matplotlib back to face 0 (Regular),
        # even when their metadata advertises weight 300. Prefer the two
        # explicitly extracted faces for this process.
        font_manager.fontManager.ttflist[:] = [
            entry for entry in font_manager.fontManager.ttflist
            if not (entry.name == "Helvetica Neue" and Path(entry.fname) == collection_path)
        ]
    except (ImportError, OSError, IndexError):
        # Helvetica Neue Regular remains a close fallback on systems where the
        # TTC cannot be split or fontTools is unavailable.
        return


register_helvetica_neue_light()

plt.rcParams.update(
    {
        "font.family": FONT,
        "font.weight": 300,
        "font.size": 12,
        "axes.labelsize": 13,
        "axes.labelweight": 300,
        "axes.titleweight": 300,
        "axes.edgecolor": EDGE,
        "axes.linewidth": 0.6,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "axes.labelcolor": TEXT,
        "text.color": TEXT,
        "savefig.dpi": 300,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def prism_csv(path: Path, table_uid: str) -> list[list[str]]:
    member = f"data/tables/{table_uid}/data.csv"
    with zipfile.ZipFile(path) as archive:
        text = archive.read(member).decode("utf-8-sig")
    return list(csv.reader(io.StringIO(text)))


def numeric_columns(rows: list[list[str]]) -> list[np.ndarray]:
    width = max(map(len, rows))
    columns = []
    for index in range(width):
        values = []
        for row in rows:
            if index < len(row) and row[index].strip():
                try:
                    values.append(float(row[index]))
                except ValueError:
                    pass
        columns.append(np.asarray(values, dtype=float))
    return columns


def image_crop(figure: int, box: tuple[int, int, int, int]) -> Image.Image:
    return Image.open(ORIGINAL / f"Fig{figure}.png").convert("RGB").crop(box)


def image_axis(ax, image: Image.Image) -> None:
    ax.imshow(image, interpolation="lanczos")
    ax.set_axis_off()


def quantitative_axis(ax) -> None:
    ax.set_facecolor(PANEL_BG)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, linewidth=1.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=0.6)


def figure_panel_label(fig, ax, label: str, *, x_pad: float = 0.035, y: float | None = None) -> None:
    """Place a panel letter in the surrounding gutter, never over panel content."""
    position = ax.get_position()
    fig.text(
        position.x0 - x_pad,
        position.y1 + 0.02 if y is None else y,
        label,
        ha="left",
        va="bottom",
        fontsize=15,
        fontweight="bold",
        color=TEXT,
    )


def jitter(n: int) -> np.ndarray:
    if n == 1:
        return np.zeros(1)
    return np.linspace(-0.085, 0.085, n)


def mean_sem_plot(
    ax,
    groups: list[np.ndarray],
    labels: list[str],
    ylabel: str,
    ylim: tuple[float, float] | None = None,
) -> None:
    quantitative_axis(ax)
    x = np.arange(len(groups), dtype=float)
    means = [np.mean(group) for group in groups]
    sems = [np.std(group, ddof=1) / math.sqrt(len(group)) for group in groups]
    ax.bar(
        x,
        means,
        yerr=sems,
        width=0.62,
        color=[WT, MUTANT] if len(groups) == 2 else [WT, MUTANT, WT, MUTANT],
        alpha=BAR_ALPHA,
        edgecolor=EDGE,
        linewidth=0.7,
        capsize=3,
        error_kw={"ecolor": EDGE, "elinewidth": 0.8, "capthick": 0.8},
        zorder=2,
    )
    for xpos, group in zip(x, groups):
        color = [WT, MUTANT] if len(groups) == 2 else [WT, MUTANT, WT, MUTANT]
        ax.scatter(xpos + jitter(len(group)), group, s=23, color=color[int(xpos)],
                   edgecolor="#000000", linewidth=0.45, alpha=POINT_ALPHA, zorder=3)
    ax.set_xticks(x, labels)
    for tick, label in zip(ax.get_xticklabels(), labels):
        if "40-1" in label:
            tick.set_fontstyle("italic")
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)


def significance(ax, left: float, right: float, y: float, text: str, height: float | None = None) -> None:
    ymin, ymax = ax.get_ylim()
    h = height if height is not None else (ymax - ymin) * 0.035
    ax.plot([left, left, right, right], [y, y + h, y + h, y], color=EDGE, linewidth=0.7)
    ax.text((left + right) / 2, y + h + (ymax - ymin) * 0.015, text, ha="center", va="bottom", fontweight=300)


def grouped_p19_plot(ax, groups: list[np.ndarray]) -> None:
    """Plot genotype by p19, using p19 on x and genotype as colour."""
    quantitative_axis(ax)
    centers = np.arange(2, dtype=float)
    offsets = {"WT": -0.18, "40-1": 0.18}
    colors = {"WT": WT, "40-1": MUTANT}
    for condition, center in enumerate(centers):
        for genotype_index, genotype in enumerate(("WT", "40-1")):
            group = groups[condition * 2 + genotype_index] / 1000
            xpos = center + offsets[genotype]
            mean = np.mean(group)
            sem = np.std(group, ddof=1) / math.sqrt(len(group))
            ax.bar(xpos, mean, width=0.32, color=colors[genotype], alpha=BAR_ALPHA,
                   edgecolor=EDGE, linewidth=0.7, zorder=2)
            ax.errorbar(xpos, mean, yerr=sem, fmt="none", ecolor=EDGE,
                        elinewidth=0.8, capsize=3, capthick=0.8, zorder=4)
            ax.scatter(xpos + jitter(len(group)) * 0.55, group, s=23,
                       color=colors[genotype], edgecolor="#000000", linewidth=0.45,
                       alpha=POINT_ALPHA, zorder=3)
    ax.set_xticks(centers, ["−p19", "+p19"])
    ax.set_ylabel("GFP intensity (AU x1000)")
    ax.set_ylim(0, 29)
    significance(ax, centers[0] - 0.18, centers[0] + 0.18, 25.2, "**")
    significance(ax, centers[1] - 0.18, centers[1] + 0.18, 25.2, "***")
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=colors[label], edgecolor=EDGE,
                      linewidth=0.8, alpha=BAR_ALPHA, label=label)
        for label in ("WT", "40-1")
    ]
    legend = ax.legend(handles=handles, frameon=False, loc="lower center",
                       bbox_to_anchor=(0.5, 1.01), ncol=2, handlelength=1.2,
                       columnspacing=1.2, borderaxespad=0)
    legend.get_texts()[1].set_fontstyle("italic")


def save_all(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure1() -> None:
    leaf = numeric_columns(prism_csv(DATA / "Leaf_count.prism", "7A956D5F-049C-4979-8718-FA5BA3E617E8"))
    # Approximate values digitised from Fig1.png; the underlying table was absent.
    fresh = [np.array([13.0, 18.0, 23.0, 39.0]), np.array([50.0, 67.0, 80.0, 83.0])]
    fig = plt.figure(figsize=(10.4, 5.25), facecolor="white")
    gs = GridSpec(2, 2, figure=fig, width_ratios=[3.1, 1], hspace=0.34, wspace=0.16)
    ax_a = fig.add_subplot(gs[:, 0])
    image_axis(ax_a, image_crop(1, (8, 35, 522, 470)))
    ax_b = fig.add_subplot(gs[0, 1])
    mean_sem_plot(ax_b, leaf[:2], ["WT", "40-1"], "Leaves per plant", (0, 210))
    significance(ax_b, 0, 1, 184, "***")
    ax_c = fig.add_subplot(gs[1, 1])
    mean_sem_plot(ax_c, fresh, ["WT", "40-1"], "Fresh weight (g)", (0, 105))
    significance(ax_c, 0, 1, 88, "**")
    fig.subplots_adjust(left=0.055, right=0.985, top=0.90, bottom=0.10)
    top_y = 0.94
    figure_panel_label(fig, ax_a, "A", x_pad=0.030, y=top_y)
    figure_panel_label(fig, ax_b, "B", x_pad=0.075, y=top_y)
    figure_panel_label(fig, ax_c, "C", x_pad=0.075, y=ax_c.get_position().y1 + 0.025)
    save_all(fig, "202608_ethylene_Fig1_reformatted")


def figure2() -> None:
    virb = numeric_columns(prism_csv(DATA / "PvirB_lux_chemiluminescence.prism", "853A02C8-0456-4D37-8536-642A4F52597B"))
    fig = plt.figure(figsize=(10.4, 4.7), facecolor="white")
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.75, 1], wspace=0.35)
    ax_a = fig.add_subplot(gs[0])
    image_axis(ax_a, image_crop(2, (36, 25, 520, 415)))
    ax_b = fig.add_subplot(gs[1])
    mean_sem_plot(ax_b, virb[:2], ["WT", "40-1"], "Chemiluminescence (AU)", (0, 1550))
    significance(ax_b, 0, 1, 1370, "*")
    fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.12)
    top_y = 0.93
    figure_panel_label(fig, ax_a, "A", x_pad=0.035, y=top_y)
    figure_panel_label(fig, ax_b, "B", x_pad=0.040, y=top_y)
    save_all(fig, "202608_ethylene_Fig2_reformatted")


def cova_elisa_plot(ax) -> None:
    """Draw the former Figure 4 quantitative panel on an existing axis."""
    rows = prism_csv(DATA / "CoVA_ELISA.prism", "4294EB22-6092-4860-93D1-C8A35680D009")
    dilutions = [row[0] for row in rows]
    wt = np.asarray([[float(value) for value in row[1:7]] for row in rows])
    mutant = np.asarray([[float(value) for value in row[7:13]] for row in rows])
    x = np.arange(len(dilutions))
    quantitative_axis(ax)
    for values, color, marker, label, offset in (
        (wt, WT, "o", "WT", -0.05),
        (mutant, MUTANT, "o", "40-1", 0.05),
    ):
        means = values.mean(axis=1)
        sems = values.std(axis=1, ddof=1) / np.sqrt(values.shape[1])
        ax.plot(x, means, color=color, linewidth=1.25, marker=marker,
                markersize=4.2, markeredgewidth=0, label=label, zorder=3)
        ax.fill_between(x, means - sems, means + sems, color=color,
                        alpha=0.20, linewidth=0, zorder=2)
        for xpos, replicates in zip(x, values):
            ax.scatter(
                np.full(len(replicates), xpos + offset) + jitter(len(replicates)) * 0.34,
                replicates, s=17, color=color, edgecolor="#000000", linewidth=0.45,
                alpha=POINT_ALPHA, zorder=4,
            )
    ax.set_yscale("log")
    ax.set_ylim(0.07, 1.35)
    ax.set_xticks(x, dilutions)
    ax.set_xlabel("Sample dilution")
    ax.set_ylabel("A₄₅₀")
    legend = ax.legend(frameon=False, loc="upper right", fontsize=11,
                       handlelength=1.5, borderaxespad=0.3)
    for text in legend.get_texts():
        if text.get_text() == "40-1":
            text.set_fontstyle("italic")


def figure3() -> None:
    gfp = numeric_columns(prism_csv(DATA / "GFP_fluorescence.prism", "E1628383-0656-4DB8-B611-52ECAC611246"))
    western = numeric_columns(prism_csv(DATA / "GFP_western.prism", "FC3AC1FA-9FE9-4094-B82A-E9C501C3BB23"))
    fig = plt.figure(figsize=(14.2, 7.6), facecolor="white")
    outer = GridSpec(2, 1, figure=fig, height_ratios=[1.34, 0.86], hspace=0.23)
    # B retains the taller top-row height but receives more width. This uses
    # space that would otherwise sit around aspect-preserved panel A.
    top = outer[0].subgridspec(1, 2, width_ratios=[2.05, 1], wspace=0.20)
    # D is squeezed horizontally, while retaining the full bottom-row height.
    bottom = outer[1].subgridspec(1, 3, width_ratios=[1.72, 0.72, 1.18], wspace=0.42)
    ax_a = fig.add_subplot(top[0, 0])
    image_axis(ax_a, image_crop(3, (26, 0, 530, 270)))
    ax_b = fig.add_subplot(top[0, 1])
    grouped_p19_plot(ax_b, gfp[:4])
    ax_c = fig.add_subplot(bottom[0, 0])
    image_axis(ax_c, image_crop(3, (28, 274, 530, 532)))
    ax_d = fig.add_subplot(bottom[0, 1])
    mean_sem_plot(ax_d, western[:2], ["WT", "40-1"], "Anti-GFP:Ponceau (AU)", (0, 4.6))
    significance(ax_d, 0, 1, 4.05, "P = 0.110")
    ax_e = fig.add_subplot(bottom[0, 2])
    cova_elisa_plot(ax_e)
    fig.subplots_adjust(left=0.050, right=0.985, top=0.89, bottom=0.105)
    top_y = 0.922
    bottom_y = max(ax_c.get_position().y1, ax_e.get_position().y1) + 0.018
    fig.text(0.014, top_y, "A", ha="left", va="bottom", fontsize=15,
             fontweight="bold", color=TEXT)
    figure_panel_label(fig, ax_b, "B", x_pad=0.060, y=top_y)
    fig.text(0.014, bottom_y, "C", ha="left", va="bottom", fontsize=15,
             fontweight="bold", color=TEXT)
    figure_panel_label(fig, ax_d, "D", x_pad=0.045, y=bottom_y)
    figure_panel_label(fig, ax_e, "E", x_pad=0.045, y=bottom_y)
    save_all(fig, "202608_ethylene_Fig3_reformatted")


def figure4() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.0), facecolor="white")
    cova_elisa_plot(ax)
    save_all(fig, "202608_ethylene_Fig4_reformatted")


def write_figure4_stats() -> None:
    """Re-run the original two-way ANOVA and exploratory per-dilution tests."""
    rows = prism_csv(DATA / "CoVA_ELISA.prism", "4294EB22-6092-4860-93D1-C8A35680D009")
    records = []
    for row in rows:
        for genotype, values in (("WT", row[1:7]), ("40-1", row[7:13])):
            for replicate, value in enumerate(values, start=1):
                records.append({"dilution": row[0], "genotype": genotype,
                                "replicate": replicate, "a450": float(value)})
    data = pd.DataFrame(records)
    model = ols("a450 ~ C(genotype) * C(dilution)", data=data).fit()
    anova = sm.stats.anova_lm(model, typ=2).reset_index().rename(columns={"index": "term"})
    LOGS.mkdir(parents=True, exist_ok=True)
    anova.to_csv(LOGS / "cova_elisa_two_way_anova.csv", index=False)
    comparisons = []
    for dilution, subset in data.groupby("dilution", sort=False):
        wt = subset.loc[subset["genotype"].eq("WT"), "a450"]
        mutant = subset.loc[subset["genotype"].eq("40-1"), "a450"]
        result = stats.ttest_ind(wt, mutant, equal_var=False)
        comparisons.append({"dilution": dilution, "test": "two-sided Welch t-test",
                            "wt_mean": wt.mean(), "mutant_mean": mutant.mean(),
                            "t": result.statistic, "df": result.df, "p_value": result.pvalue,
                            "significant_p_lt_0.05": result.pvalue < 0.05})
    pd.DataFrame(comparisons).to_csv(LOGS / "cova_elisa_per_dilution_welch_tests.csv", index=False)

def main() -> None:
    figure1()
    figure2()
    figure3()
    figure4()
    write_figure4_stats()
    write_methods_note()


if __name__ == "__main__":
    main()
