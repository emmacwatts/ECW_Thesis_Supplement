"""Cas9 WT versus core-704-6 GFP signal five days after GUS infiltration."""

from pathlib import Path
import os
import sys

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

SOURCE_DIR = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/VIGSinCORE/cas9vCORE_GUS_5dpi"
)
DATA_FILE = SOURCE_DIR / "GUSinfil_core_wt.xlsx"
PHOTO_DIR = HERE / "source_photos"
PHOTOS = {
    "Cas9 WT": PHOTO_DIR / "cas9_wt_gus.png",
    "core-704-6": PHOTO_DIR / "core_704_6_gus.png",
}
OUT = HERE / "202607_corevigs_gus_5dpi_figures"
LOG = HERE / "202607_corevigs_gus_5dpi_logs"
for folder in (OUT, LOG):
    folder.mkdir(parents=True, exist_ok=True)

ORDER = ["Cas9 WT", "core-704-6"]
COLORS = {"Cas9 WT": "#A6B0BA", "core-704-6": "#287FB8"}
BAR_OUTLINE = "#6F7780"

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 17,
    "axes.labelsize": 18,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})


def load_data():
    source = pd.read_excel(DATA_FILE, sheet_name="Sheet1", header=1, usecols="B:E")
    source = source.rename(columns={
        "Sample": "source_group", "Background": "background",
        "GFP": "gfp", "Normalised GFP": "source_normalised_gfp",
    }).dropna(subset=["source_group", "background", "gfp"])
    source["genotype"] = source["source_group"].replace({
        "WT": "Cas9 WT", "704-6": "core-704-6",
    })
    source["normalised_gfp"] = source["gfp"] - source["background"]
    source["replicate"] = source.groupby("genotype").cumcount() + 1
    return source


def analyse(data):
    summary = (
        data.groupby("genotype")["normalised_gfp"]
        .agg(n="size", mean="mean", sd="std", sem="sem")
        .reindex(ORDER).reset_index()
    )
    wt = data.loc[data["genotype"].eq("Cas9 WT"), "normalised_gfp"]
    core = data.loc[data["genotype"].eq("core-704-6"), "normalised_gfp"]
    result = stats.ttest_ind(wt, core, equal_var=True, alternative="two-sided")
    test = pd.DataFrame([{
        "comparison": "Cas9 WT vs core-704-6",
        "test": "ordinary two-sided unpaired t-test",
        "t": result.statistic,
        "df": len(wt) + len(core) - 2,
        "p": result.pvalue,
    }])
    with pd.ExcelWriter(LOG / "corevigs_gus_5dpi_processed_data_and_statistics.xlsx") as writer:
        data.to_excel(writer, sheet_name="processed_data", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        test.to_excel(writer, sheet_name="unpaired_t_test", index=False)
    return summary, test.iloc[0]


def add_plot(ax, data, summary, test):
    x = np.arange(len(ORDER))
    means = summary["mean"].to_numpy() / 1000
    sems = summary["sem"].to_numpy() / 1000
    ax.bar(
        x, means, width=0.58, color=[COLORS[group] for group in ORDER], alpha=0.60,
        edgecolor="#303030", linewidth=1.2, yerr=sems, capsize=5,
        error_kw={"elinewidth": 1.5, "ecolor": TEXT_COLOR, "capthick": 1.5},
    )
    rng = np.random.default_rng(20260715)
    for index, group in enumerate(ORDER):
        values = data.loc[data["genotype"].eq(group), "normalised_gfp"].to_numpy() / 1000
        jitter = rng.uniform(-0.10, 0.10, len(values))
        point_x = np.full(len(values), index) + jitter
        ax.scatter(point_x, values, s=42, color=COLORS[group], alpha=0.82,
                   edgecolor="#202020", linewidth=0.65, zorder=4)

    ax.set_xticks(x, ORDER)
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_ylabel("Normalised GFP signal (×10³ a.u.)")
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    ax.spines[["top", "right"]].set_visible(False)

    ymax = max((data["normalised_gfp"] / 1000).max(), (means + sems).max())
    bracket_bottom = ymax + 0.10
    bracket_top = ymax + 0.25
    ax.plot([0, 0, 1, 1], [bracket_bottom, bracket_top, bracket_top, bracket_bottom],
            color=TEXT_COLOR, lw=1.4, clip_on=False)
    ax.text(0.5, bracket_top + 0.06, "ns", ha="center", va="bottom", color=TEXT_COLOR)
    ax.set_ylim(0, bracket_top + 0.42)


def add_photo_panel(ax):
    ax.set_facecolor("black")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    image_axes = [
        ax.inset_axes([0.0, 0.070, 0.555, 0.86]),
        ax.inset_axes([0.595, 0.155, 0.35, 0.68]),
    ]
    for image_ax, group in zip(image_axes, ORDER):
        image = Image.open(PHOTOS[group]).convert("RGBA")
        if group == "Cas9 WT":
            # The supplied crop has a small screenshot-overlay icon at the bottom.
            image = image.crop((0, 0, image.width, round(image.height * 0.90)))
        alpha_bounds = image.getchannel("A").getbbox()
        if alpha_bounds:
            image = image.crop(alpha_bounds)
        image.thumbnail((1500, 2000), Image.Resampling.LANCZOS)
        image_ax.imshow(image)
        image_ax.set_facecolor("black")
        image_ax.set_anchor("S")
        image_ax.set_axis_off()

    label_style = {
        "transform": ax.transAxes, "ha": "center", "va": "center",
        "color": "white", "fontsize": 18, "fontfamily": "Helvetica Neue",
        "fontweight": 300,
    }
    ax.text(0.255, 0.925, "Cas9 WT", **label_style)
    ax.text(0.745, 0.925, "core-704-6", fontstyle="italic", **label_style)

    # The user-supplied calibration line spans approximately 13.3% of the
    # displayed photo-panel width and represents 10 cm. Draw the final scale
    # bar in the lower-right black margin so it does not obscure either plant.
    bar_right = 0.955
    bar_left = bar_right - 0.133
    bar_y = 0.060
    ax.plot(
        [bar_left, bar_right], [bar_y, bar_y], transform=ax.transAxes,
        color="white", linewidth=4.0, solid_capstyle="butt",
        clip_on=False, zorder=20,
    )
    ax.text(
        (bar_left + bar_right) / 2, bar_y + 0.020, "10 cm",
        transform=ax.transAxes, ha="center", va="bottom",
        color="white", fontsize=14, fontfamily="Helvetica Neue",
        fontweight=300, clip_on=False, zorder=20,
    )


def make_figure(data, summary, test):
    fig = plt.figure(figsize=(13.6, 6.9), facecolor="white")
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.06, 0.94], left=0.085, right=0.985,
        bottom=0.14, top=0.89, wspace=0.12,
    )
    photo_ax = fig.add_subplot(gs[0, 0])
    plot_ax = fig.add_subplot(gs[0, 1])
    add_plot(plot_ax, data, summary, test)
    add_photo_panel(photo_ax)
    for ax, letter in [(photo_ax, "A"), (plot_ax, "B")]:
        box = ax.get_position()
        fig.text(
            box.x0 - 0.012, 0.975, letter, ha="left", va="top", fontsize=20,
            fontweight="bold", fontfamily="Helvetica Neue", color=TEXT_COLOR,
        )
    stem = "202607_corevigs_gus_5dpi_cas9_wt_vs_core_704_6_with_photos"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    data = load_data()
    summary, test = analyse(data)
    make_figure(data, summary, test)
+ ". "
        "Normalised GFP was calculated as GFP minus background. Bars show mean ± SEM "
        "with all three measurements overlaid. Cas9 WT and core-704-6 were compared "
        "using the ordinary two-sided unpaired t-test specified in the source Prism "
        f"project (t = {test.t:.4f}, df = {int(test.df)}, P = {test.p:.4f}).\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(test.to_string())


if __name__ == "__main__":
    main()
