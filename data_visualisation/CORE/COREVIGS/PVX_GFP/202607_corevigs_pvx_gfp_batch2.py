"""Batch 2 PVX-GFP signal in WT and core-1-2 at 5 and 12 dpi."""

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

STYLE_DIR = HERE.parents[2] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/PVX_GFP/PVX-GFP_batch1and2/PVX_GFP.xlsx"
)
OUT = HERE / "202607_corevigs_pvx_gfp_batch2_figures"
LOG = HERE / "202607_corevigs_pvx_gfp_batch2_logs"
PHOTO_DIR = HERE / "source_photos"
for folder in (OUT, LOG, PHOTO_DIR):
    folder.mkdir(parents=True, exist_ok=True)

ORDER = ["WT", "core-1-2"]
COLORS = {"WT": "#A6B0BA", "core-1-2": "#287FB8"}
BAR_OUTLINE = "#6F7780"

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": 16, "axes.labelsize": 17, "axes.titlesize": 18,
    "xtick.labelsize": 15, "ytick.labelsize": 15,
    "figure.dpi": 120, "savefig.dpi": 300,
})


def read_sheet(sheet, dpi):
    data = pd.read_excel(SOURCE, sheet_name=sheet, header=1)
    data = data.rename(columns={
        "Sample": "source_genotype", "GFP type": "construct",
        "Background": "background", "GFP Signal": "gfp_signal",
    })
    data = data[data["source_genotype"].isin(["wt", "core-1-2"])].copy()
    data["genotype"] = data["source_genotype"].replace({"wt": "WT"})
    data["construct"] = data.get("construct", pd.Series(index=data.index, dtype=object)).fillna("PVX-GFP")
    data["construct"] = data["construct"].replace({"PVX": "PVX-GFP"})
    data["normalised_gfp"] = data["gfp_signal"] - data["background"]
    data["dpi"] = dpi
    data["replicate"] = data.groupby(["genotype", "construct"]).cumcount() + 1
    return data[[
        "dpi", "construct", "genotype", "replicate", "background",
        "gfp_signal", "normalised_gfp",
    ]]


def load_data():
    data = pd.concat([
        read_sheet("batch2_5dpi", 5), read_sheet("batch2_12dpi", 12)
    ], ignore_index=True)
    return data[data["construct"].eq("PVX-GFP")].copy()


def analyse(data):
    summary = (
        data.groupby(["dpi", "construct", "genotype"])["normalised_gfp"]
        .agg(n="size", mean="mean", sd="std", sem="sem").reset_index()
    )
    rows = []
    for (dpi, construct), subset in data.groupby(["dpi", "construct"], sort=True):
        wt = subset.loc[subset["genotype"].eq("WT"), "normalised_gfp"].to_numpy()
        core = subset.loc[subset["genotype"].eq("core-1-2"), "normalised_gfp"].to_numpy()
        result = stats.ttest_ind(wt, core, equal_var=False, alternative="two-sided")
        rows.append({
            "dpi": dpi, "construct": construct, "comparison": "WT vs core-1-2",
            "test": "two-sided Welch t-test", "n_wt": len(wt), "n_core": len(core),
            "t": result.statistic, "df": result.df, "p": result.pvalue,
        })
    tests = pd.DataFrame(rows)
    with pd.ExcelWriter(LOG / "corevigs_pvx_gfp_batch2_processed_data_and_statistics.xlsx") as writer:
        data.to_excel(writer, sheet_name="processed_data", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        tests.to_excel(writer, sheet_name="welch_tests", index=False)
    return summary, tests


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    apply_axis_style(ax)


def grouped_plot(ax, data, summary, tests):
    timepoints = [5, 12]
    centers = np.arange(len(timepoints), dtype=float)
    offsets = {"WT": -0.18, "core-1-2": 0.18}
    width = 0.34
    rng = np.random.default_rng(20260715)
    for genotype in ORDER:
        group_summary = summary[
            summary["genotype"].eq(genotype)
        ].set_index("dpi").loc[timepoints]
        positions = centers + offsets[genotype]
        ax.bar(
            positions, group_summary["mean"] / 1000, width=width,
            yerr=group_summary["sem"] / 1000, capsize=4,
            color=COLORS[genotype], alpha=0.60, edgecolor="#303030", linewidth=1.2,
            error_kw={"elinewidth": 1.3, "ecolor": TEXT_COLOR, "capthick": 1.3},
            label=genotype,
        )
        for position, dpi in zip(positions, timepoints):
            values = data[
                data["dpi"].eq(dpi) & data["genotype"].eq(genotype)
            ]["normalised_gfp"].to_numpy() / 1000
            point_x = np.full(len(values), position) + rng.uniform(-0.045, 0.045, len(values))
            ax.scatter(point_x, values, s=42, color=COLORS[genotype], alpha=0.82,
                       edgecolor="#202020", linewidth=0.65, zorder=4)
    ax.set_xticks(centers, ["5 dpi", "12 dpi"])
    ax.set_ylabel("Normalised GFP signal (×10³ a.u.)")
    ax.legend(frameon=False, loc="upper right")
    style_axis(ax)
    top = max(14.5, data["normalised_gfp"].max() / 1000 + 2.0)
    for index, dpi in enumerate(timepoints):
        subset = data[data["dpi"].eq(dpi)]
        y = subset["normalised_gfp"].max() / 1000 + 0.55
        ax.plot(
            [index - 0.18, index - 0.18, index + 0.18, index + 0.18],
            [y - 0.15, y, y, y - 0.15], color=TEXT_COLOR, lw=1.2,
        )
        ax.text(index, y + 0.12, "ns", ha="center", va="bottom", fontsize=14)
    ax.set_ylim(0, top)


def photo_panel(ax):
    ax.set_facecolor("black")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    paths = [PHOTO_DIR / "wt.png", PHOTO_DIR / "core-1-2.png"]
    image_axes = [
        ax.inset_axes([0.045, 0.13, 0.425, 0.62]),
        ax.inset_axes([0.53, 0.13, 0.425, 0.62]),
    ]
    for image_ax, path in zip(image_axes, paths):
        image = Image.open(path).convert("RGBA")
        bounds = image.getchannel("A").getbbox()
        if bounds:
            image = image.crop(bounds)
        image.thumbnail((1500, 1500), Image.Resampling.LANCZOS)
        image_ax.imshow(image)
        image_ax.set_facecolor("black")
        image_ax.set_anchor("C")
        image_ax.set_axis_off()
    label_style = {
        "transform": ax.transAxes, "ha": "center", "va": "center",
        "color": "white", "fontsize": 18, "fontfamily": "Helvetica Neue",
        "fontweight": 300,
    }
    ax.text(0.258, 0.79, "WT", **label_style)
    ax.text(0.742, 0.79, "core-1-2", fontstyle="italic", **label_style)


def make_figure(data, summary, tests):
    fig = plt.figure(figsize=(13.6, 6.9), facecolor="white")
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.12, 0.88], left=0.085, right=0.985,
        bottom=0.14, top=0.89, wspace=0.12,
    )
    axes = [fig.add_subplot(gs[0, index]) for index in range(2)]
    photo_panel(axes[0])
    grouped_plot(axes[1], data, summary, tests)
    for ax, letter in zip(axes, "AB"):
        box = ax.get_position()
        fig.text(
            box.x0 - 0.012, 0.965, letter, ha="left", va="top", fontsize=20,
            fontweight="bold", fontfamily="Helvetica Neue", color=TEXT_COLOR,
        )
    stem = "202607_corevigs_pvx_gfp_batch2_5dpi_12dpi_with_photos"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    data = load_data()
    summary, tests = analyse(data)
    make_figure(data, summary, tests)
print(summary.to_string(index=False))
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
