"""Rebuild the two six-week NPR1 GFP/COVA figure sets."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
LOG = HERE / "logs"
OUT.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "1.General/ReportsandPapers/IssyVIGS/NPR1_FollowOn_Exps/6wk_GFPCOVA"
)
STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, ERRORBAR_COLOR, FONT_STACK, TEXT_COLOR, apply_axis_style  # noqa: E402

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 14,
    "axes.labelsize": 15,
    "xtick.labelsize": 13.5,
    "ytick.labelsize": 13.5,
    "legend.fontsize": 11.5,
    "savefig.dpi": 300,
})

WT_COLOR = "#666666"
NPR1_COLOR = "#315F78"
POINT_ALPHA = 0.82


def read_gfp(path):
    data = pd.read_excel(path)
    data.columns = [str(c).strip() for c in data.columns]
    value_col = next(c for c in data.columns if c.replace(" ", "").lower() == "normalisedgfp")
    data = data.rename(columns={"Plant Type": "genotype", value_col: "relative_gfp"})
    data["genotype"] = data["genotype"].astype(str).str.lower().map(
        lambda x: "WT" if x == "wt" else "npr1−"
    )
    return data[["genotype", "Repeat", "Replicate", "relative_gfp"]]


def plot_gfp(ax, data, p_value, scale_thousands=False):
    order = ["WT", "npr1−"]
    colors = [WT_COLOR, NPR1_COLOR]
    rng = np.random.default_rng(202609)
    ymax = 0
    for x, (group, color) in enumerate(zip(order, colors)):
        values = data.loc[data["genotype"] == group, "relative_gfp"].to_numpy(float)
        mean, sd = values.mean(), values.std(ddof=1)
        ax.bar(x, mean, width=0.58, color=color, alpha=0.576, edgecolor="black", linewidth=1.0, zorder=2)
        ax.errorbar(x, mean, yerr=sd, fmt="none", ecolor=ERRORBAR_COLOR,
                    elinewidth=1.1, capsize=4, capthick=1.1, zorder=5)
        jitter = rng.uniform(-0.12, 0.12, len(values))
        ax.scatter(x + jitter, values, s=54, facecolor=color, edgecolor="#202020",
                   linewidth=0.65, alpha=POINT_ALPHA, zorder=6)
        ymax = max(ymax, values.max(), mean + sd)
    ax.set_xticks([0, 1], ["WT", r"$\it{npr1}$"])
    ax.set_ylabel(r"Relative GFP intensity ($\times 10^3$)" if scale_thousands else "Relative GFP intensity")
    if scale_thousands:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1000:g}"))
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, ymax * 1.35)
    label = "ns" if p_value >= 0.05 else "****"
    y = ymax * 1.10
    ax.plot([0, 0, 1, 1], [y * 0.96, y, y, y * 0.96], color=TEXT_COLOR, lw=1.0, clip_on=False)
    ax.text(0.5, y * 1.01, label, ha="center", va="bottom", fontsize=14, fontweight="bold" if label != "ns" else "normal")
    apply_axis_style(ax)


def representative_elisa():
    # Repeat 3 values from the workbook Pivot table (right-hand block, K5:R8).
    dilution = np.array([1 / 3000, 1 / 300, 1 / 30, 1 / 3])
    series = {
        ("npr1", "COVA"): [0.0762, 0.08625, 0.1878, 0.2593],
        ("WT", "COVA"): [0.06615, 0.0763, 0.18555, 0.2414],
        ("npr1", "EV"): [0.0689, 0.06425, 0.0618, 0.0678],
        ("WT", "EV"): [0.0557, 0.0554, 0.0720, 0.0581],
    }
    rows = []
    for (genotype, infiltration), values in series.items():
        for d, value in zip(dilution, values):
            rows.append({"repeat": 3, "genotype": genotype, "infiltration": infiltration,
                         "sample_dilution": d, "a450": value})
    return pd.DataFrame(rows)


def plot_elisa(ax, data, legend_above=False):
    x = np.arange(4)
    labels = ["1:3000", "1:300", "1:30", "1:3"]
    for genotype, color, marker in [("WT", WT_COLOR, "o"), ("npr1", NPR1_COLOR, "s")]:
        for infiltration, linestyle, alpha in [("COVA", "-", 1.0), ("EV", "--", 0.72)]:
            values = data.loc[(data.genotype == genotype) & (data.infiltration == infiltration), "a450"]
            ax.plot(x, values, color=color, marker=marker, markersize=5.5,
                    markerfacecolor=("none" if genotype == "npr1" else color),
                    markeredgecolor=color, markeredgewidth=1.3,
                    linewidth=1.8, linestyle=linestyle, alpha=alpha,
                    label=(f"{genotype}, {infiltration}" if genotype == "WT"
                           else rf"$\it{{npr1}}$, {infiltration}"))
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.set_xlabel("Sample dilution")
    ax.set_ylabel(r"Absorbance ($A_{450}$)")
    ax.set_ylim(0, 0.32)
    legend_kwargs = ({"bbox_to_anchor": (0, 1.04, 1, 0.1), "loc": "lower left",
                      "mode": "expand", "borderaxespad": 0}
                     if legend_above else {"loc": "upper left"})
    ax.legend(frameon=False, ncol=(4 if legend_above else 2), handlelength=2.2,
              columnspacing=1.1, handletextpad=0.5, **legend_kwargs)
    apply_axis_style(ax)


def crop_western(source, box):
    with Image.open(source) as image:
        return image.convert("RGB").crop(box)


def show_leaf_panel(ax, leaf):
    # Separate the source halves with black space and append a full-width footer.
    width, height = leaf.size
    half = width // 2
    expanded_width = int(round(width * 1.32))
    footer = int(round(height * 0.18))
    canvas = Image.new("RGB", (expanded_width, height + footer), "black")
    left_x = int(round(expanded_width * 0.25 - half * 0.5))
    right_width = width - half
    right_x = int(round(expanded_width * 0.75 - right_width * 0.5))
    canvas.paste(leaf.crop((0, 0, half, height)), (left_x, 0))
    canvas.paste(leaf.crop((half, 0, width, height)), (right_x, 0))
    ax.imshow(canvas)
    ax.text(0.25, 0.07, "WT", transform=ax.transAxes, color="white", fontsize=15,
            ha="center", va="center")
    ax.text(0.75, 0.07, r"$\it{npr1}$",
            transform=ax.transAxes, color="white", fontsize=15,
            ha="center", va="center")
    ax.axis("off")


def show_western_panel(ax, western, ponceau, ponceau_has_outline=False):
    ax.imshow(western, extent=(0, 8, 1.05, 2.35), aspect="auto")
    ax.imshow(ponceau, extent=(0, 8, 0.12, 0.82), aspect="auto")
    ax.add_patch(Rectangle((0, 1.05), 8, 1.30, fill=False, edgecolor=TEXT_COLOR,
                           linewidth=1.0, zorder=8, clip_on=False))
    if not ponceau_has_outline:
        ax.add_patch(Rectangle((0, 0.12), 8, 0.70, fill=False, edgecolor=TEXT_COLOR,
                               linewidth=1.0, zorder=8, clip_on=False))
    for x, label in zip([1, 3, 5, 7], ["WT", r"$\it{npr1}$", "WT", r"$\it{npr1}$"]):
        ax.text(x, 3.12, label, ha="center", va="bottom", fontsize=15)
    ax.text(2, 3.72, "EV", ha="center", va="bottom", fontsize=15)
    ax.text(6, 3.72, "GFP", ha="center", va="bottom", fontsize=15)
    for x, replicate in zip(np.arange(0.5, 8, 1), [1, 2] * 4):
        ax.text(x, 2.62, str(replicate), ha="center", va="bottom", fontsize=13)
    ax.text(-0.28, 1.70, "25", ha="right", va="center", fontsize=14)
    ax.text(-0.28, 0.47, "55", ha="right", va="center", fontsize=14)
    ax.text(8, -0.18, r"$\it{anti}$–GFP", ha="right", va="top", fontsize=13)
    ax.set_xlim(-0.65, 8.05)
    ax.set_ylim(-0.42, 4.05)
    ax.axis("off")


def panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=20, fontweight="bold", ha="left", va="top")


def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight", facecolor="white")


def make_no_p19(gfp):
    # The leaf image is extracted once from the PowerPoint source by main().
    leaf = Image.open(HERE / "assets" / "leaf_gfp_localisation.png").convert("RGB")
    western = crop_western(
        SOURCE / "20241217_npr1Nop19" / "Screenshot 2024-12-21 at 14.05.22.png",
        (455, 1080, 1538, 1242),
    )
    ponceau = crop_western(
        SOURCE / "20241217_npr1Nop19" / "Screenshot 2024-12-21 at 14.05.22.png",
        (455, 1306, 1538, 1405),
    )
    fig = plt.figure(figsize=(12.0, 8.2))
    outer = fig.add_gridspec(1, 2, width_ratios=[0.76, 1.24], left=0.08, right=0.985,
                             top=0.94, bottom=0.09, wspace=0.30)
    right = outer[0, 1].subgridspec(2, 1, height_ratios=[0.86, 1.14], hspace=0.42)
    ax_a = fig.add_subplot(outer[0, 0]); plot_gfp(ax_a, gfp, 4.064056439741591e-08, scale_thousands=True)
    ax_b = fig.add_subplot(right[0, 0]); show_leaf_panel(ax_b, leaf)
    ax_c = fig.add_subplot(right[1, 0]); show_western_panel(ax_c, western, ponceau)
    ax_a.text(-0.28, 1.04, "A", transform=ax_a.transAxes, fontsize=20, fontweight="bold")
    fig.text(0.45, 0.94, "B", fontsize=20, fontweight="bold")
    fig.text(0.45, 0.525, "C", fontsize=20, fontweight="bold")
    save(fig, "202609_npr1_6wk_no_p19")
    plt.close(fig)


def make_p19(gfp, elisa):
    western = crop_western(
        SOURCE / "20241217_npr1Nop19" / "Screenshot 2024-12-21 at 13.45.59.png",
        (330, 496, 2554, 881),
    )
    ponceau = crop_western(
        SOURCE / "20241217_npr1Nop19" / "Screenshot 2024-12-21 at 13.45.59.png",
        (330, 974, 2554, 1199),
    )
    fig = plt.figure(figsize=(12.0, 8.2))
    outer = fig.add_gridspec(1, 2, width_ratios=[0.76, 1.24], left=0.08, right=0.985,
                             top=0.94, bottom=0.09, wspace=0.42)
    right = outer[0, 1].subgridspec(2, 1, height_ratios=[0.90, 1.10], hspace=0.62)
    ax_a = fig.add_subplot(outer[0, 0]); plot_gfp(ax_a, gfp, 0.4030375557367415)
    ax_b = fig.add_subplot(right[0, 0]); show_western_panel(ax_b, western, ponceau)
    ax_c = fig.add_subplot(right[1, 0]); plot_elisa(ax_c, elisa, legend_above=True)
    c_position = ax_c.get_position()
    ax_c.set_position([c_position.x0, c_position.y0 + 0.065,
                       c_position.width, c_position.height])
    fig.text(0.0, 0.965, "A", fontsize=20, fontweight="bold")
    fig.text(0.45, 0.958, "B", fontsize=20, fontweight="bold")
    ax_c.text(-0.16, 1.04, "C", transform=ax_c.transAxes, fontsize=20, fontweight="bold")
    save(fig, "202609_npr1_6wk_p19")
    plt.close(fig)


def extract_leaf_asset():
    import zipfile
    assets = HERE / "assets"; assets.mkdir(exist_ok=True)
    target = assets / "leaf_gfp_localisation.png"
    with zipfile.ZipFile(SOURCE / "6wk_NPR1_GFP_COVA_FinalPlots.pptx") as archive:
        target.write_bytes(archive.read("ppt/media/image8.png"))


def main():
    extract_leaf_asset()
    no_p19 = read_gfp(SOURCE / "20241217_npr1Nop19" / "NPR1_nop19.xlsx")
    p19 = read_gfp(SOURCE / "20241211" / "npr1_6wk_p19GFPCOVA.xlsx")
    elisa = representative_elisa()
    no_p19.to_csv(LOG / "no_p19_gfp_cleaned.csv", index=False)
    p19.to_csv(LOG / "p19_gfp_cleaned.csv", index=False)
    elisa.to_csv(LOG / "p19_elisa_repeat3.csv", index=False)
    make_no_p19(no_p19)
    make_p19(p19, elisa)
print(f"Figures saved to {OUT}")


if __name__ == "__main__":
    main()
