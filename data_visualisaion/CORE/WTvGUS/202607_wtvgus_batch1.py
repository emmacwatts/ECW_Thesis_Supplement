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
    "8. Miscellaneous/WTvGUS"
)
DATA_FILE = SOURCE_DIR / "Batch1" / "WTvGUS_6wk_20240909.xlsx"
PHOTO_DIR = HERE / "source_photos"
PHOTO = {
    ("4 wk", "WT"): PHOTO_DIR / "4wk_WT.png",
    ("4 wk", "GUS"): PHOTO_DIR / "4wk_GUS.png",
    ("6 wk", "WT"): PHOTO_DIR / "6wk_WT.png",
    ("6 wk", "GUS"): PHOTO_DIR / "6wk_GUS.png",
}
OUT = HERE / "202607_wtvgus_batch1_figures"
LOG = HERE / "202607_wtvgus_batch1_logs"
OUT.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

COLORS = {"WT": "#8E9AA7", "GUS": "#D68C68"}
ORDER = ["WT", "GUS"]
DISPLAY_LABELS = {"WT": "WT", "GUS": "TRV2::GUS"}
FORMATS = ("png", "pdf", "svg")

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 17,
    "axes.labelsize": 18,
    "axes.titlesize": 19,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})


def load_data():
    data = pd.read_excel(DATA_FILE, sheet_name="Sheet1")
    data = data.rename(columns={
        "VIGS Target": "genotype",
        "Replicate": "measurement",
        "Background": "background",
        "GFP Signal": "gfp_signal",
        "GFP Signal - Background": "normalized_gfp",
    })
    data = data.loc[data["genotype"].isin(ORDER)].copy()
    for col in ["measurement", "background", "gfp_signal", "normalized_gfp"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    return data.dropna(subset=["normalized_gfp"])


def analyse(data):
    wt = data.loc[data.genotype.eq("WT"), "normalized_gfp"].to_numpy()
    gus = data.loc[data.genotype.eq("GUS"), "normalized_gfp"].to_numpy()
    test = stats.ttest_ind(wt, gus, equal_var=False, alternative="two-sided")
    summary = data.groupby("genotype", sort=False).normalized_gfp.agg(
        n="size", mean="mean", sd="std", sem="sem"
    ).reindex(ORDER).reset_index()
    stats_out = pd.DataFrame([{
        "comparison": "WT vs GUS",
        "test": "two-sided Welch t-test",
        "t": test.statistic,
        "p": test.pvalue,
        "df": getattr(test, "df", np.nan),
    }])
    with pd.ExcelWriter(LOG / "batch1_6wk_processed_data_and_statistics.xlsx") as writer:
        data.to_excel(writer, sheet_name="source_data", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        stats_out.to_excel(writer, sheet_name="welch_test", index=False)
    return summary, stats_out.iloc[0]


def add_plot(ax, data, summary, test):
    x = np.arange(2)
    means = summary["mean"].to_numpy() / 1000
    sems = summary["sem"].to_numpy() / 1000
    ax.bar(x, means, width=0.58, color=[COLORS[g] for g in ORDER], alpha=0.82,
           edgecolor="#6F7780", linewidth=0.8, yerr=sems, capsize=5,
           error_kw={"elinewidth": 1.5, "ecolor": TEXT_COLOR, "capthick": 1.5})
    rng = np.random.default_rng(260715)
    for i, group in enumerate(ORDER):
        vals = data.loc[data.genotype.eq(group), "normalized_gfp"].to_numpy() / 1000
        jitter = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=31, facecolor=COLORS[group],
                   edgecolor="none", alpha=0.55, zorder=4)
    ax.set_xticks(x, [DISPLAY_LABELS[g] for g in ORDER])
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_ylabel("Normalised GFP signal (×10³ a.u.)")
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    ax.spines[["top", "right"]].set_visible(False)
    ymax = max((data.normalized_gfp / 1000).max(), (means + sems).max())
    y0, y1 = ymax + 0.25, ymax + 0.55
    ax.plot([0, 0, 1, 1], [y0, y1, y1, y0], color=TEXT_COLOR, lw=1.5, clip_on=False)
    label = "****" if test.p < 0.0001 else f"P = {test.p:.3g}"
    ax.text(0.5, y1 + 0.08, label, ha="center", va="bottom", color=TEXT_COLOR)
    ax.set_ylim(0, y1 + 0.85)


def photo_axis(ax, path, label, row_label=None, italic=False):
    image = Image.open(path).convert("RGB")
    image.thumbnail((1400, 2100), Image.Resampling.LANCZOS)
    ax.imshow(image)
    ax.set_facecolor("black")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(0.5, 0.965, label, transform=ax.transAxes, ha="center", va="top",
            color="white", fontsize=24, fontstyle="italic" if italic else "normal")
    if row_label:
        ax.text(-0.05, 0.5, row_label, transform=ax.transAxes, rotation=90,
                ha="right", va="center", color=TEXT_COLOR, fontsize=17)


def add_photo_panel(ax, left_path, right_path):
    """Draw both photographs inside one full-height black panel."""
    ax.set_facecolor("black")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Keep a compact label band above larger photographs.
    image_axes = [
        ax.inset_axes([0.015, 0.025, 0.477, 0.83]),
        ax.inset_axes([0.508, 0.025, 0.477, 0.83]),
    ]
    for image_ax, path in zip(image_axes, [left_path, right_path]):
        image = Image.open(path).convert("RGB")
        image.thumbnail((1400, 2100), Image.Resampling.LANCZOS)
        image_ax.imshow(image)
        image_ax.set_facecolor("black")
        image_ax.set_axis_off()

    label_style = {
        "transform": ax.transAxes,
        "ha": "center",
        "va": "center",
        "color": "white",
        "fontsize": 18,
        "fontfamily": "Helvetica Neue",
        "fontweight": 300,
    }
    ax.text(0.254, 0.925, DISPLAY_LABELS["WT"], **label_style)
    ax.text(0.746, 0.925, DISPLAY_LABELS["GUS"], fontstyle="italic", **label_style)


def save(fig, stem):
    for fmt in FORMATS:
        fig.savefig(OUT / f"{stem}.{fmt}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_6wk(data, summary, test):
    fig = plt.figure(figsize=(13.6, 6.9), constrained_layout=False, facecolor="white")
    gs = fig.add_gridspec(
        1, 2, width_ratios=[0.94, 1.06], left=0.085, right=0.985,
        bottom=0.14, top=0.89, wspace=0.12,
    )
    ax = fig.add_subplot(gs[0, 0])
    photo_panel = fig.add_subplot(gs[0, 1])
    add_plot(ax, data, summary, test)
    add_photo_panel(
        photo_panel, PHOTO[("6 wk", "WT")], PHOTO[("6 wk", "GUS")]
    )

    # Panel letters share the same baseline in the clear upper margin.
    for panel_ax, letter in [(ax, "A"), (photo_panel, "B")]:
        box = panel_ax.get_position()
        fig.text(
            box.x0 - 0.012, 0.975, letter, ha="left", va="top", fontsize=20,
            fontweight="bold", fontfamily="Helvetica Neue", color=TEXT_COLOR,
        )
    save(fig, "batch1_6wk_wt_vs_gus_with_photos")


def make_age_context(data, summary, test):
    fig = plt.figure(figsize=(9.2, 15.8), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[0.82, 1, 1], hspace=0.08, wspace=0.05)
    ax = fig.add_subplot(gs[0, :])
    add_plot(ax, data, summary, test)
    for row, age in enumerate(["4 wk", "6 wk"], start=1):
        for col, genotype in enumerate(ORDER):
            photo_axis(fig.add_subplot(gs[row, col]), PHOTO[(age, genotype)],
                       genotype if row == 1 else "", age if col == 0 else None)
    save(fig, "batch1_4wk_6wk_phenotype_context")


if __name__ == "__main__":
    data = load_data()
    summary, test = analyse(data)
    make_6wk(data, summary, test)
    (LOG / "methods_note.txt").write_text(
        "Batch 1, 6-week background-subtracted GFP values were read from " + str(DATA_FILE) +
        ". Bars show mean ± SEM with all source measurements overlaid. WT and GUS were "
        "compared using the same two-sided Welch t-test used in the source Prism analysis. "
        "The figure includes the corresponding 6-week Batch 1 WT and TRV2::GUS phenotype "
        "photographs in a shared black panel.\n"
    )
