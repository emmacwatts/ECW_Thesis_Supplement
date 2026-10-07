# %% [markdown]
# # Combined CORE GFP endpoint
#
# Plot the combined-batch endpoint workbook using the shared thesis figure style.

# %%
from pathlib import Path
import os
import sys

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from PIL import Image
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

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "combinedGFP"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/CombinedTest/GFP_20250227/Combined_Batch_GFP.xlsx"
)
SPOT_DIR = HERE / "combined_spots"
SHEET = "Sheet1"
SAMPLE_ORDER = ["WT", "2-8", "2-1", "704-6", "714-1"]
SAMPLE_LABELS = {
    "WT": "Cas9 WT",
    "2-8": "core-2-8",
    "2-1": "core-1-2",
    "704-6": "core-704-6",
    "714-1": "core-714-1",
}
EXCLUDED_SAMPLES = {"705-1"}
SAVE_FMTS = ["png", "pdf", "svg"]
BAR_COLOR = "#6F8FAF"
WT_COLOR = "#6D6D6D"
BAR_ALPHA = 0.74
BAR_WIDTH = 0.58
SIG_ALPHA = 0.05
Y_SCALE = 1000

SPOT_FILES = {
    "WT": "WT.png",
    "2-8": "2-8.png",
    "2-1": "1-2.png",
    "704-6": "704-6.png",
    "714-1": "714-1.png",
}

STATS_METHODS = (
    "Stats/processing: source data were read from Sheet1 of Combined_Batch_GFP.xlsx "
    "using the background-subtracted Normalised GFP column as the response. CORE "
    "705-1 was excluded before plotting and statistics; CORE 2-1 is displayed as "
    "core-1-2 for consistency with the other CORE analyses. Bars show means +/- SEM "
    "with individual biological values overlaid and n annotated above each bar. Each "
    "CORE line was compared with WT using a two-sided Welch two-sample t-test, with "
    "Holm adjustment across the four comparisons. The Cas9 WT control is shown in "
    "grey and CORE lines in muted blue. No outlier removal or transformation "
    "was applied. The plotted y-axis is divided by 1,000 and labelled as AU x1000."
)

plt.rcParams.update(
    {
        **ANALYSIS_RCPARAMS,
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 15,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    }
)


def sem(values):
    values = pd.Series(values).dropna()
    return 0.0 if len(values) <= 1 else float(values.std(ddof=1) / np.sqrt(len(values)))


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


def load_signal_spots(spot_dir):
    """Trim black margins and return equally framed, user-prepared spot images."""
    spots = {}
    for sample, filename in SPOT_FILES.items():
        image = np.asarray(Image.open(spot_dir / filename).convert("RGB"))
        foreground = np.max(image, axis=2) > 8
        rows, columns = np.where(foreground)
        if not len(rows):
            raise ValueError(f"No spot signal found in {spot_dir / filename}")
        crop = image[rows.min() : rows.max() + 1, columns.min() : columns.max() + 1]
        side = max(crop.shape[:2])
        square = np.zeros((side, side, 3), dtype=np.uint8)
        y0 = (side - crop.shape[0]) // 2
        x0 = (side - crop.shape[1]) // 2
        square[y0 : y0 + crop.shape[0], x0 : x0 + crop.shape[1]] = crop
        spots[sample] = square
    return spots


def load_raw_data(src):
    raw = pd.read_excel(src, sheet_name=SHEET, header=4).iloc[:, 1:5].copy()
    raw = raw.rename(
        columns={
            "Sample Type": "sample_type",
            "Background": "background",
            "GFP": "gfp",
            "Normalised GFP": "normalized_gfp",
        }
    )
    for column in ["background", "gfp", "normalized_gfp"]:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    raw["sample_type"] = raw["sample_type"].astype(str).str.strip()
    raw = raw.dropna(subset=["sample_type", "normalized_gfp"]).copy()
    raw["source_file"] = str(src)
    raw["source_sheet"] = SHEET
    return raw.reset_index(drop=True)


def analysis_data(raw):
    data = raw.loc[
        raw["sample_type"].isin(SAMPLE_ORDER)
        & ~raw["sample_type"].isin(EXCLUDED_SAMPLES)
    ].copy()
    data["sample_type"] = pd.Categorical(data["sample_type"], SAMPLE_ORDER, ordered=True)
    data["replicate"] = data.groupby("sample_type", observed=True).cumcount() + 1
    return data.sort_values(["sample_type", "replicate"]).reset_index(drop=True)


def summary_table(data):
    return (
        data.groupby("sample_type", observed=True)
        .agg(
            n=("normalized_gfp", "size"),
            mean_normalized_gfp=("normalized_gfp", "mean"),
            sem_normalized_gfp=("normalized_gfp", sem),
            sd_normalized_gfp=("normalized_gfp", lambda x: pd.Series(x).std(ddof=1)),
        )
        .reset_index()
    )


def run_stats(data):
    wt = data.loc[data["sample_type"].eq("WT"), "normalized_gfp"].dropna()
    rows = []
    for sample in SAMPLE_ORDER[1:]:
        mutant = data.loc[data["sample_type"].eq(sample), "normalized_gfp"].dropna()
        result = stats.ttest_ind(wt, mutant, equal_var=False, alternative="two-sided")
        rows.append(
            {
                "comparison": f"WT vs {SAMPLE_LABELS[sample]}",
                "test": "two-sided Welch t-test",
                "wt_n": len(wt),
                "core_n": len(mutant),
                "wt_mean": wt.mean(),
                "core_mean": mutant.mean(),
                "estimate_core_minus_wt": mutant.mean() - wt.mean(),
                "t_statistic": result.statistic,
                "p_raw": result.pvalue,
                "sample_type": sample,
            }
        )
    table = pd.DataFrame(rows)
    _, adjusted, _, _ = multipletests(table["p_raw"], method="holm")
    table["p_holm"] = adjusted
    table["annotation"] = table["p_holm"].map(p_to_label)
    return table


def draw_bars(ax, data, stats_table, show_xlabels=True, font_increase=0):
    x_positions = np.arange(len(SAMPLE_ORDER), dtype=float)
    tops = {}

    for x, sample in zip(x_positions, SAMPLE_ORDER):
        values = data.loc[data["sample_type"].eq(sample), "normalized_gfp"].dropna() / Y_SCALE
        mean = values.mean()
        error = sem(values)
        fill_color = WT_COLOR if sample == "WT" else BAR_COLOR
        ax.bar(x, mean, width=BAR_WIDTH, color=fill_color, alpha=BAR_ALPHA,
               edgecolor="black", linewidth=1.0, zorder=2)
        ax.errorbar(x, mean, yerr=error, fmt="none", ecolor=ERRORBAR_COLOR,
                    elinewidth=1.1, capsize=3, capthick=1.1, zorder=5)
        jitter = np.linspace(-0.07, 0.07, len(values)) if len(values) > 1 else np.array([0.0])
        ax.scatter(x + jitter, values, s=25, facecolor=fill_color, edgecolor="none",
                   alpha=INDIVIDUAL_POINT_ALPHA, zorder=4)
        label_y = max(mean + error, values.max()) + 0.08
        ax.text(x, label_y, str(len(values)), ha="center", va="bottom",
                fontsize=14 + font_increase, color="#8A8A8A")
        tops[sample] = label_y

    annotation_y = max(tops.values()) + 0.22
    for row_index, row in stats_table.iterrows():
        mutant_x = SAMPLE_ORDER.index(row["sample_type"])
        y = annotation_y + row_index * 0.16
        ax.plot([0, mutant_x], [y, y], color=TEXT_COLOR, linewidth=0.9, clip_on=False)
        ax.text(mutant_x / 2, y + 0.025, row["annotation"], ha="center", va="bottom",
                fontsize=15 + font_increase, color=TEXT_COLOR, clip_on=False)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([SAMPLE_LABELS[x] for x in SAMPLE_ORDER] if show_xlabels else [])
    for label, sample in zip(ax.get_xticklabels(), SAMPLE_ORDER):
        if sample != "WT":
            label.set_fontstyle("italic")
    ax.set_ylabel("Relative GFP intensity (AU x1000)")
    ax.set_ylim(0, annotation_y + len(stats_table) * 0.16 + 0.18)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
    style_axis(ax)


def plot_combined_gfp(data, stats_table):
    fig, ax = plt.subplots(figsize=(9.2, 5.6), constrained_layout=False)
    draw_bars(ax, data, stats_table)
    fig.subplots_adjust(left=0.14, right=0.98, top=0.94, bottom=0.20)
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{SCRIPT_STEM}_combined_gfp.{file_format}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_combined_gfp_with_dots(data, stats_table):
    dots = load_signal_spots(SPOT_DIR)
    fig = plt.figure(figsize=(10.4, 8.4), constrained_layout=False)
    grid = GridSpec(2, 1, figure=fig, height_ratios=[4.25, 1.55], hspace=0.10)
    ax = fig.add_subplot(grid[0])
    dot_ax = fig.add_subplot(grid[1], sharex=ax)
    draw_bars(ax, data, stats_table, show_xlabels=False, font_increase=4)
    ax.yaxis.label.set_fontsize(20)
    ax.tick_params(axis="y", labelsize=19)
    ax.tick_params(axis="x", bottom=False)

    dot_ax.set_facecolor("black")
    dot_ax.set_xlim(-0.5, len(SAMPLE_ORDER) - 0.5)
    dot_ax.set_ylim(0, 1)
    dot_ax.set_yticks([])
    dot_ax.set_xticks([])
    for spine in dot_ax.spines.values():
        spine.set_visible(False)

    half_width = 0.40
    for x, sample in enumerate(SAMPLE_ORDER):
        dot_ax.imshow(
            dots[sample],
            extent=(x - half_width, x + half_width, 0.28, 0.94),
            interpolation="lanczos",
            aspect="auto",
            zorder=2,
        )
        dot_ax.text(
            x,
            0.12,
            SAMPLE_LABELS[sample],
            ha="center",
            va="center",
            fontsize=18,
            color="white",
            fontstyle="italic" if sample != "WT" else "normal",
            zorder=3,
        )

    fig.subplots_adjust(left=0.15, right=0.98, top=0.965, bottom=0.07)
    for file_format in SAVE_FMTS:
        fig.savefig(
            SAVE_DIR / f"{SCRIPT_STEM}_combined_gfp_with_signal_dots.{file_format}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def write_outputs(raw, data, summary, stats_table):
    raw.to_excel(LOG_DIR / "combined_gfp_raw_long.xlsx", index=False)
    data.to_excel(LOG_DIR / "combined_gfp_analysis_long.xlsx", index=False)
    summary.to_excel(LOG_DIR / "combined_gfp_summary_means.xlsx", index=False)
    stats_table.to_excel(LOG_DIR / "combined_gfp_welch_tests.xlsx", index=False)
    (LOG_DIR / "combined_gfp_methods_note.txt").write_text(STATS_METHODS + "\n", encoding="utf-8")


# %%
def main():
    raw = load_raw_data(SRC)
    data = analysis_data(raw)
    summary = summary_table(data)
    stats_table = run_stats(data)
    write_outputs(raw, data, summary, stats_table)
    plot_combined_gfp(data, stats_table)
    plot_combined_gfp_with_dots(data, stats_table)
    print(f"Analysis rows: {len(data)}")
    print(f"Excluded rows (705-1): {(raw['sample_type'] == '705-1').sum()}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")
    print(stats_table[["comparison", "p_raw", "p_holm", "annotation"]].to_string(index=False))
    return raw, data, summary, stats_table


if __name__ == "__main__":
    main()
