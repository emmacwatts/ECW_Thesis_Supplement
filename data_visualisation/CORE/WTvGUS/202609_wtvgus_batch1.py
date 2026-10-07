from pathlib import Path
import os
import sys

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "8. Miscellaneous/WTvGUS/Batch1/WTvGUS_6wk_20240909.xlsx"
)
STEM = Path(__file__).stem
OUT = HERE / f"{STEM}_figures"
LOG = HERE / f"{STEM}_logs"
OUT.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

ORDER = ["WT", "GUS"]
BAR_COLOR = "#86A5C1"
FORMATS = ("png", "pdf", "svg")

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 18,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})


def load_data():
    data = pd.read_excel(SOURCE, sheet_name="Sheet1")
    data = data.rename(columns={
        "VIGS Target": "genotype",
        "Replicate": "measurement",
        "Background": "background",
        "GFP Signal": "gfp_signal",
        "GFP Signal - Background": "normalised_gfp",
    })
    data = data.loc[data.genotype.isin(ORDER)].copy()
    for column in ["measurement", "background", "gfp_signal", "normalised_gfp"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(subset=["normalised_gfp"])
    data["source_file"] = str(SOURCE)
    data["source_sheet"] = "Sheet1"
    return data.reset_index(drop=True)


def analyse(data):
    summary = (
        data.groupby("genotype")["normalised_gfp"]
        .agg(n="size", mean="mean", sd="std", sem="sem")
        .reindex(ORDER).reset_index()
    )
    wt = data.loc[data.genotype.eq("WT"), "normalised_gfp"]
    gus = data.loc[data.genotype.eq("GUS"), "normalised_gfp"]
    test = stats.ttest_ind(wt, gus, equal_var=False, alternative="two-sided")
    result = pd.DataFrame([{
        "comparison": "WT vs GUS",
        "test": "two-sided Welch t-test",
        "wt_n": len(wt), "gus_n": len(gus),
        "wt_mean": wt.mean(), "gus_mean": gus.mean(),
        "t": test.statistic, "df": getattr(test, "df", np.nan), "p": test.pvalue,
    }])
    return summary, result


def p_label(p):
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def make_figure(data, summary, result):
    fig, ax = plt.subplots(figsize=(7.6, 5.4), facecolor="white")
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.96)
    x = np.arange(2)
    means = summary["mean"].to_numpy() / 1000
    sems = summary["sem"].to_numpy() / 1000
    ax.bar(
        x, means, yerr=sems, width=0.60, color=BAR_COLOR, alpha=0.9,
        edgecolor="#66727D", linewidth=1.25, capsize=5,
        error_kw={"elinewidth": 1.5, "ecolor": TEXT_COLOR, "capthick": 1.5},
    )
    rng = np.random.default_rng(260928)
    for i, genotype in enumerate(ORDER):
        values = data.loc[data.genotype.eq(genotype), "normalised_gfp"].to_numpy() / 1000
        jitter = rng.uniform(-0.105, 0.105, len(values))
        ax.scatter(
            np.full(len(values), i) + jitter, values, s=34,
            facecolor="#587C9C", edgecolor="white", linewidth=0.55,
            alpha=0.82, zorder=4,
        )

    ax.set_xticks(x, ["WT", "TRV2::GUS"])
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_ylabel("Relative GFP intensity\n(a.u. × 10³)")
    ax.set_xlim(-0.62, 1.62)
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    ax.spines[["top", "right"]].set_visible(False)

    ymax = max((data.normalised_gfp / 1000).max(), (means + sems).max())
    bracket_y = ymax + 0.85
    tick = 0.34
    ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
            color=TEXT_COLOR, lw=1.3, clip_on=False)
    ax.text(0.5, bracket_y + 0.10, p_label(float(result.iloc[0].p)),
            ha="center", va="bottom", fontsize=18)
    ax.set_ylim(0, bracket_y + 1.05)

    for fmt in FORMATS:
        fig.savefig(OUT / f"{STEM}_wt_vs_trv2_gus.{fmt}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_logs(data, summary, result):
    with pd.ExcelWriter(LOG / "wtvgus_batch1_processed_data_and_statistics.xlsx") as writer:
        data.to_excel(writer, sheet_name="source_values", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        result.to_excel(writer, sheet_name="welch_test", index=False)
if __name__ == "__main__":
    data = load_data()
    summary, result = analyse(data)
    save_logs(data, summary, result)
    make_figure(data, summary, result)
    print(summary.to_string(index=False))
    print(result.to_string(index=False))
