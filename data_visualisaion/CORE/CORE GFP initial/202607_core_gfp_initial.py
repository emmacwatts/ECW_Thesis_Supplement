# %% [markdown]
# # CORE initial GFP endpoint
#
# Replot the 2026-03-30 CORE GFP endpoint data from the source workbook using
# the current shared thesis figure style.

# %%
from pathlib import Path
import os
import re
import sys

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
else:
    HERE = Path.cwd() / "CORE" / "CORE GFP initial"

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))

from figure_styles import (  # noqa: E402
    ANALYSIS_RCPARAMS,
    ERRORBAR_COLOR,
    FONT_STACK,
    INDIVIDUAL_POINT_ALPHA,
    TEXT_COLOR,
    apply_axis_style,
)

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_core_gfp_initial"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/20260330/20260330.xlsx"
)

SOURCE_PRISM = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/20260330/20260330.prism"
)

SAVE_FMTS = ["png", "pdf", "svg"]
SHEETS = ["batch1", "batch2"]
DISPLAY_COMPARISONS = [
    {"batch": "batch1", "p19": "no p19", "label": "No p19"},
    {"batch": "batch1", "p19": "p19", "label": "p19"},
]
SAMPLE_ORDER = ["cas9 WT", "core 1-2"]
SAMPLE_LABELS = {"cas9 WT": "Cas9 WT", "core 1-2": "core-1-2"}
BAR_COLOR = "#6F8FAF"
SIG_ALPHA = 0.05
BAR_ALPHA = 0.74
BAR_WIDTH = 0.34
GROUP_GAP = 0.52
Y_SCALE = 1000

STATS_METHODS = (
    "Stats/processing: source data were read from 20260330.xlsx, using the "
    "background-subtracted Normalised GFP column as the response. The plotted figure "
    "uses batch1 only for both no p19 and p19 conditions, comparing Cas9 WT with "
    "core 1-2. Non-CORE rows, including prp27 #1, "
    "were excluded before plotting and statistics. Bars show means +/- SEM with "
    "individual biological values overlaid and n annotated above each bar. Within "
    "each displayed batch x p19 condition, Cas9 WT and core 1-2 were compared using "
    "standard two-sided Welch two-sample t-tests; p-values were Holm-adjusted across "
    "the two displayed comparisons. Non-significant comparisons are annotated as ns "
    "on the plot. No outlier removal or transformation was applied. "
    f"The previous Prism file used as a processing reference was: {SOURCE_PRISM}."
)

CORE_RCPARAMS = {
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 18,
    "axes.titlesize": 19,
    "axes.labelsize": 18,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "legend.fontsize": 17,
    "legend.title_fontsize": 17,
    "figure.dpi": 120,
    "savefig.dpi": 300,
}
plt.rcParams.update(CORE_RCPARAMS)


# %%
def sem(values):
    values = pd.Series(values).dropna()
    if len(values) <= 1:
        return 0.0
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def p_to_label(p_value):
    if pd.isna(p_value):
        return ""
    if p_value >= SIG_ALPHA:
        return "ns"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3, direction="out")
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.88)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def load_raw_data(src):
    frames = []
    for sheet_name in SHEETS:
        df = pd.read_excel(src, sheet_name=sheet_name, header=1).iloc[:, 1:].copy()
        df = df.rename(
            columns={
                "Sample Type": "sample_type",
                "p19?": "p19",
                "Background": "background",
                "GFP": "gfp",
                "Normalised GFP": "normalized_gfp",
            }
        )
        for column in ["background", "gfp", "normalized_gfp"]:
            df[column] = pd.to_numeric(df[column], errors="coerce")
        df["sample_type"] = df["sample_type"].astype(str).str.strip()
        df["p19"] = df["p19"].astype(str).str.strip()
        df["batch"] = sheet_name
        df["source_file"] = str(src)
        df["source_sheet"] = sheet_name
        frames.append(df)
    raw = pd.concat(frames, ignore_index=True)
    return raw.dropna(subset=["sample_type", "p19", "normalized_gfp"]).reset_index(drop=True)


def analysis_data(raw):
    keep_conditions = {(item["batch"], item["p19"]) for item in DISPLAY_COMPARISONS}
    data = raw.loc[
        raw["sample_type"].isin(SAMPLE_ORDER)
        & raw[["batch", "p19"]].apply(tuple, axis=1).isin(keep_conditions)
    ].copy()
    data["condition"] = data["batch"] + " " + data["p19"]
    data["sample_type"] = pd.Categorical(data["sample_type"], categories=SAMPLE_ORDER, ordered=True)
    data["replicate"] = data.groupby(["batch", "p19", "sample_type"], observed=True).cumcount() + 1
    return data.sort_values(["batch", "p19", "sample_type", "replicate"]).reset_index(drop=True)


def summary_table(data):
    return (
        data.groupby(["batch", "p19", "sample_type"], observed=True)
        .agg(
            n=("normalized_gfp", "size"),
            mean_normalized_gfp=("normalized_gfp", "mean"),
            sem_normalized_gfp=("normalized_gfp", sem),
            sd_normalized_gfp=("normalized_gfp", lambda values: pd.Series(values).std(ddof=1)),
        )
        .reset_index()
    )


def run_stats(data):
    rows = []
    for item in DISPLAY_COMPARISONS:
        subset = data.loc[data["batch"].eq(item["batch"]) & data["p19"].eq(item["p19"])]
        wt = subset.loc[subset["sample_type"].eq("cas9 WT"), "normalized_gfp"].dropna()
        core = subset.loc[subset["sample_type"].eq("core 1-2"), "normalized_gfp"].dropna()
        if len(wt) >= 2 and len(core) >= 2:
            result = stats.ttest_ind(wt, core, equal_var=False, alternative="two-sided")
            t_statistic = float(result.statistic)
            p_raw = float(result.pvalue)
        else:
            t_statistic = np.nan
            p_raw = np.nan
        rows.append(
            {
                "batch": item["batch"],
                "p19": item["p19"],
                "comparison": "Cas9 WT vs core 1-2",
                "test": "two-sided Welch t-test",
                "wt_n": int(len(wt)),
                "core_1_2_n": int(len(core)),
                "wt_mean": float(wt.mean()) if len(wt) else np.nan,
                "core_1_2_mean": float(core.mean()) if len(core) else np.nan,
                "wt_sem": sem(wt),
                "core_1_2_sem": sem(core),
                "estimate_core_1_2_minus_wt": float(core.mean() - wt.mean()) if len(wt) and len(core) else np.nan,
                "t_statistic": t_statistic,
                "p_raw": p_raw,
            }
        )

    table = pd.DataFrame(rows)
    table["p_holm_displayed_family"] = np.nan
    valid = table["p_raw"].notna()
    if valid.any():
        _, adjusted, _, _ = multipletests(table.loc[valid, "p_raw"], method="holm")
        table.loc[valid, "p_holm_displayed_family"] = adjusted
    table["annotation"] = table["p_holm_displayed_family"].map(p_to_label)
    return table


def bar_layout():
    positions = {}
    xticks = []
    xlabels = []
    condition_centers = []
    x = 0.0
    for condition_index, condition in enumerate(DISPLAY_COMPARISONS):
        condition_positions = []
        for sample in SAMPLE_ORDER:
            positions[(condition["batch"], condition["p19"], sample)] = x
            xticks.append(x)
            xlabels.append(SAMPLE_LABELS[sample])
            condition_positions.append(x)
            x += BAR_WIDTH + 0.32
        condition_centers.append(float(np.mean(condition_positions)))
        if condition_index < len(DISPLAY_COMPARISONS) - 1:
            x += GROUP_GAP
    return positions, xticks, xlabels, condition_centers


def add_n_label(ax, x, y, label):
    ax.text(
        x,
        y,
        str(label),
        ha="center",
        va="bottom",
        fontsize=16,
        color="#8A8A8A",
        zorder=8,
    )


def add_stat_label(ax, x1, x2, y, label):
    if not label:
        return y
    ax.plot([x1, x2], [y, y], color=TEXT_COLOR, linewidth=0.9, clip_on=False)
    ax.text(
        (x1 + x2) / 2,
        y + 0.45,
        label,
        ha="center",
        va="bottom",
        fontsize=18,
        color=TEXT_COLOR,
        clip_on=False,
    )
    return y + 1.15


def plot_initial_gfp(data, summary, stats_table):
    fig, ax = plt.subplots(figsize=(8.2, 5.4), constrained_layout=False)
    positions, xticks, xlabels, condition_centers = bar_layout()
    ymax = 0.0
    n_label_tops = {}

    for condition in DISPLAY_COMPARISONS:
        for sample in SAMPLE_ORDER:
            subset = data.loc[
                data["batch"].eq(condition["batch"])
                & data["p19"].eq(condition["p19"])
                & data["sample_type"].eq(sample)
            ]
            values = subset["normalized_gfp"].dropna()
            if values.empty:
                continue
            scaled_values = values / Y_SCALE
            mean = values.mean() / Y_SCALE
            err = sem(values) / Y_SCALE
            x = positions[(condition["batch"], condition["p19"], sample)]
            ax.bar(
                x,
                mean,
                width=BAR_WIDTH,
                color=BAR_COLOR,
                alpha=BAR_ALPHA,
                edgecolor="black",
                linewidth=1.0,
                zorder=2,
            )
            ax.errorbar(
                x,
                mean,
                yerr=err,
                fmt="none",
                ecolor=ERRORBAR_COLOR,
                elinewidth=1.1,
                capsize=3.0,
                capthick=1.1,
                zorder=5,
            )
            jitter = np.linspace(-0.045, 0.045, len(values)) if len(values) > 1 else np.array([0.0])
            ax.scatter(
                x + jitter,
                scaled_values,
                s=24,
                facecolor=BAR_COLOR,
                edgecolor="none",
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=4,
            )
            label_y = max(mean + err, scaled_values.max()) + 0.70
            add_n_label(ax, x, label_y, len(values))
            n_label_tops[(condition["batch"], condition["p19"])] = max(
                n_label_tops.get((condition["batch"], condition["p19"]), 0.0),
                label_y,
            )
            ymax = max(ymax, label_y + 0.85)

    star_ymax = ymax
    for condition in DISPLAY_COMPARISONS:
        row = stats_table.loc[
            stats_table["batch"].eq(condition["batch"]) & stats_table["p19"].eq(condition["p19"])
        ]
        if row.empty or not row["annotation"].iloc[0]:
            continue
        x1 = positions[(condition["batch"], condition["p19"], "cas9 WT")]
        x2 = positions[(condition["batch"], condition["p19"], "core 1-2")]
        y = n_label_tops.get((condition["batch"], condition["p19"]), 0.0) + 1.95
        star_ymax = max(star_ymax, add_stat_label(ax, x1, x2, y, row["annotation"].iloc[0]))

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    for tick_label, sample_label in zip(ax.get_xticklabels(), xlabels):
        if sample_label.startswith("core"):
            tick_label.set_fontstyle("italic")

    for center, condition in zip(condition_centers, DISPLAY_COMPARISONS):
        ax.text(
            center,
            -0.165,
            condition["label"],
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=17,
            color=TEXT_COLOR,
            fontfamily=FONT_STACK,
        )

    ax.set_ylabel("Relative GFP intensity (AU x1000)")
    y_upper = max(30, np.ceil(star_ymax * 1.14 / 10) * 10)
    ax.set_ylim(0, y_upper)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    style_axis(ax)

    fig.subplots_adjust(left=0.17, right=0.98, top=0.93, bottom=0.25)
    save(fig, "core_gfp_initial")
    plt.close(fig)


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{output_stem}.{file_format}", dpi=300, bbox_inches="tight")


def write_outputs(raw, data, summary, stats_table):
    raw.to_excel(LOG_DIR / "core_gfp_initial_raw_long.xlsx", index=False)
    data.to_excel(LOG_DIR / "core_gfp_initial_analysis_long.xlsx", index=False)
    summary.to_excel(LOG_DIR / "core_gfp_initial_summary_means.xlsx", index=False)
    stats_table.to_excel(LOG_DIR / "core_gfp_initial_welch_tests.xlsx", index=False)
    (LOG_DIR / "core_gfp_initial_methods_note.txt").write_text(STATS_METHODS + "\n", encoding="utf-8")


# %%
def main():
    raw = load_raw_data(SRC)
    data = analysis_data(raw)
    summary = summary_table(data)
    stats_table = run_stats(data)
    write_outputs(raw, data, summary, stats_table)
    plot_initial_gfp(data, summary, stats_table)
    print(f"Analysis rows: {len(data)}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")
    print(stats_table[["batch", "p19", "wt_n", "core_1_2_n", "p_raw", "p_holm_displayed_family", "annotation"]].to_string(index=False))
    return raw, data, summary, stats_table


if __name__ == "__main__":
    main()
