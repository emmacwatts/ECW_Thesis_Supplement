# %% [markdown]
# # CORE OE batch 1 GFP analysis
#
# Recreate the CORE overexpression batch 1 plot from the raw workbook rows,
# excluding selected plant lines and using the shared thesis figure style.

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
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, TextArea
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
else:
    HERE = Path.cwd() / "CORE" / "CORE OE"

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

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_core_oe_batch1"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/CombinedTest/LMU_LargeBatch_20250426/"
    "core exp/20250426_CORE.xlsx"
)

SHEET = "CORE OE"
BATCH = 1
EXCLUDED_PLANT_TYPES = {"cas9 WT", "2-8", "705-1"}
PLANT_ORDER = ["714-1", "704-6", "1-2"]
PLANT_PLOT_LABELS = ["core-714-1", "core-704-6", "core-1-2"]
TREATMENT_ORDER = ["EV", "CORE OE"]
TREATMENT_COLORS = {
    "EV": "#4E79A7",  # muted blue
    "CORE OE": "#D77A73",  # dulled salmon
}
SAVE_FMTS = ["png", "pdf", "svg"]
SIG_ALPHA = 0.05
BAR_ALPHA = 0.74
BAR_WIDTH = 0.34
Y_SCALE = 1000

STATS_METHODS = (
    "Stats/processing: source data were read from the CORE OE sheet of "
    "20250426_CORE.xlsx. Raw rows were filtered to batch 1, middle leaf CORE OE/EV "
    "measurements, with Cas9 WT and plant types 2-8 and 705-1 excluded before "
    "plotting or statistics. The response variable was Normalised GFP, calculated in the source "
    "workbook as GFP minus background. Bars show plant-type treatment means +/- SEM, "
    "with individual replicate values overlaid and n annotated above each bar. "
    "The primary/plotted comparison tests the pre-specified directional expectation "
    "that CORE OE reduces GFP signal relative to EV within each plant type. This used "
    "one-sided Welch two-sample t-tests testing EV > CORE OE, with Holm correction "
    "across the batch 1 plant-type comparison family. Two-sided Welch p-values are "
    "retained in the stats table as sensitivity checks. No outlier removal or "
    "transformation was applied."
)

CORE_RCPARAMS = {
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 16,
    "axes.titlesize": 17,
    "axes.labelsize": 16,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 15,
    "legend.title_fontsize": 15,
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


def p_to_stars(p_value):
    if pd.isna(p_value) or p_value >= SIG_ALPHA:
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


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
    raw = pd.read_excel(src, sheet_name=SHEET).iloc[:, :8].copy()
    raw = raw.rename(
        columns={
            "Plant Type": "plant_type",
            "Leaf Age": "leaf_age",
            "Batch": "batch",
            "Plant": "plant",
            "pre-infil": "treatment",
            "Background": "background",
            "GFP": "gfp",
            "Normalised GFP": "normalized_gfp",
        }
    )
    for column in ["batch", "plant", "background", "gfp", "normalized_gfp"]:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    raw["plant_type"] = raw["plant_type"].astype(str).str.strip()
    raw["treatment"] = raw["treatment"].astype(str).str.strip()
    raw["leaf_age"] = raw["leaf_age"].astype(str).str.strip()
    raw = raw.dropna(subset=["plant_type", "batch", "plant", "treatment", "normalized_gfp"])
    raw = raw[raw["treatment"].isin(TREATMENT_ORDER)].copy()
    raw["source_file"] = str(src)
    raw["source_sheet"] = SHEET
    return raw.reset_index(drop=True)


def analysis_data(raw):
    data = raw.loc[
        raw["batch"].eq(BATCH)
        & raw["plant_type"].isin(PLANT_ORDER)
        & ~raw["plant_type"].isin(EXCLUDED_PLANT_TYPES)
    ].copy()
    data["plant_type"] = pd.Categorical(data["plant_type"], categories=PLANT_ORDER, ordered=True)
    data["treatment"] = pd.Categorical(data["treatment"], categories=TREATMENT_ORDER, ordered=True)
    data["replicate"] = data.groupby(["plant_type", "treatment"], observed=True).cumcount() + 1
    return data.sort_values(["plant_type", "plant", "treatment"]).reset_index(drop=True)


def summary_table(data):
    return (
        data.groupby(["plant_type", "treatment"], observed=True)
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
    for plant_type, subset in data.groupby("plant_type", observed=True):
        ev = subset.loc[subset["treatment"].eq("EV"), "normalized_gfp"].dropna()
        core = subset.loc[subset["treatment"].eq("CORE OE"), "normalized_gfp"].dropna()
        if len(ev) >= 2 and len(core) >= 2:
            primary = stats.ttest_ind(ev, core, equal_var=False, alternative="greater")
            two_sided = stats.ttest_ind(ev, core, equal_var=False, alternative="two-sided")
            t_statistic = float(primary.statistic)
            p_raw = float(primary.pvalue)
            p_two_sided = float(two_sided.pvalue)
        else:
            t_statistic = np.nan
            p_raw = np.nan
            p_two_sided = np.nan

        rows.append(
            {
                "batch": BATCH,
                "plant_type": str(plant_type),
                "comparison": "EV > CORE OE",
                "test": "one-sided Welch t-test",
                "sensitivity_test": "two-sided Welch t-test",
                "ev_n": int(len(ev)),
                "core_oe_n": int(len(core)),
                "ev_mean": float(ev.mean()) if len(ev) else np.nan,
                "core_oe_mean": float(core.mean()) if len(core) else np.nan,
                "ev_sem": sem(ev),
                "core_oe_sem": sem(core),
                "estimate_core_oe_minus_ev": float(core.mean() - ev.mean()) if len(ev) and len(core) else np.nan,
                "t_statistic": t_statistic,
                "p_raw": p_raw,
                "p_welch_two_sided": p_two_sided,
            }
        )

    table = pd.DataFrame(rows)
    table["p_holm_batch_family"] = np.nan
    valid = table["p_raw"].notna()
    if valid.any():
        _, adjusted, _, _ = multipletests(table.loc[valid, "p_raw"], method="holm")
        table.loc[valid, "p_holm_batch_family"] = adjusted
    table["stars"] = table["p_holm_batch_family"].map(p_to_stars)
    return table


def bar_positions():
    centers = np.arange(len(PLANT_ORDER), dtype=float)
    return {"EV": centers - BAR_WIDTH / 2, "CORE OE": centers + BAR_WIDTH / 2}, centers


def add_n_label(ax, x, y, label):
    ax.text(
        x,
        y,
        str(label),
        ha="center",
        va="bottom",
        fontsize=14,
        color="#8A8A8A",
        zorder=8,
    )


def add_star_label(ax, x1, x2, y, stars):
    if not stars:
        return y
    ax.plot([x1, x2], [y, y], color=TEXT_COLOR, linewidth=0.9, clip_on=False)
    ax.text(
        (x1 + x2) / 2,
        y + 0.45,
        stars,
        ha="center",
        va="bottom",
        fontsize=16,
        color=TEXT_COLOR,
        clip_on=False,
    )
    return y + 1.15


def legend_item(color, label_parts):
    swatch = DrawingArea(46, 18, 0, 0)
    swatch.add_artist(
        plt.Rectangle(
            (0, 2),
            38,
            14,
            facecolor=color,
            edgecolor="black",
            linewidth=1.0,
            alpha=BAR_ALPHA,
        )
    )
    text_areas = [
        TextArea(
            text,
            textprops={
                "color": TEXT_COLOR,
                "fontfamily": FONT_STACK,
                "fontsize": CORE_RCPARAMS["legend.fontsize"],
                "fontstyle": fontstyle,
            },
        )
        for text, fontstyle in label_parts
    ]
    label = HPacker(children=text_areas, align="center", pad=0, sep=0)
    return HPacker(children=[swatch, label], align="center", pad=0, sep=8)


def add_manual_legend(ax):
    ev_item = legend_item(TREATMENT_COLORS["EV"], [("EV", "normal")])
    core_item = legend_item(TREATMENT_COLORS["CORE OE"], [("CORE", "italic"), (" OE", "normal")])
    legend_box = HPacker(children=[ev_item, core_item], align="center", pad=0, sep=54)
    anchored = AnchoredOffsetbox(
        loc="lower center",
        child=legend_box,
        pad=0,
        frameon=False,
        bbox_to_anchor=(0.5, -0.29),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    ax.add_artist(anchored)


def plot_batch1(data, summary, stats_table):
    fig, ax = plt.subplots(figsize=(8.6, 5.2), constrained_layout=False)
    positions, centers = bar_positions()
    ymax = 0.0
    n_label_tops = {}

    for treatment in TREATMENT_ORDER:
        for idx, plant_type in enumerate(PLANT_ORDER):
            subset = data.loc[data["plant_type"].eq(plant_type) & data["treatment"].eq(treatment)]
            values = subset["normalized_gfp"].dropna()
            if values.empty:
                continue
            scaled_values = values / Y_SCALE
            mean = values.mean() / Y_SCALE
            err = sem(values) / Y_SCALE
            x = positions[treatment][idx]
            color = TREATMENT_COLORS[treatment]
            ymax = max(ymax, float(scaled_values.max()), float(mean + err))
            ax.bar(
                x,
                mean,
                width=BAR_WIDTH * 0.82,
                color=color,
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
            jitter = np.linspace(-0.035, 0.035, len(values)) if len(values) > 1 else np.array([0.0])
            ax.scatter(
                x + jitter,
                scaled_values,
                s=24,
                facecolor=color,
                edgecolor="none",
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=4,
            )
            label_y = max(mean + err, scaled_values.max()) + 0.70
            add_n_label(ax, x, label_y, len(values))
            n_label_tops[plant_type] = max(n_label_tops.get(plant_type, 0.0), label_y)
            ymax = max(ymax, label_y + 0.85)

    star_ymax = ymax
    for idx, plant_type in enumerate(PLANT_ORDER):
        row = stats_table.loc[stats_table["plant_type"].eq(plant_type)]
        if row.empty or not row["stars"].iloc[0]:
            continue
        group_data = data.loc[data["plant_type"].eq(plant_type)]
        y = max(group_data["normalized_gfp"].max() / Y_SCALE, n_label_tops.get(plant_type, 0.0)) + 1.95
        star_ymax = max(star_ymax, add_star_label(ax, positions["EV"][idx], positions["CORE OE"][idx], y, row["stars"].iloc[0]))

    ax.set_xticks(centers)
    ax.set_xticklabels(PLANT_PLOT_LABELS)
    ax.set_ylabel("Relative GFP intensity (AU x1000)")
    ax.set_ylim(0, max(20, star_ymax * 1.10))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax)
    for tick_label in ax.get_xticklabels():
        tick_label.set_fontstyle("italic")
    add_manual_legend(ax)

    fig.subplots_adjust(left=0.14, right=0.98, top=0.95, bottom=0.25)
    save(fig, "batch1_core_oe_relative_gfp")
    plt.close(fig)


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{output_stem}.{file_format}", dpi=300, bbox_inches="tight")


def write_outputs(raw, data, summary, stats_table):
    raw.to_excel(LOG_DIR / "core_oe_raw_long.xlsx", index=False)
    data.to_excel(LOG_DIR / "core_oe_batch1_analysis_long.xlsx", index=False)
    summary.to_excel(LOG_DIR / "core_oe_batch1_summary_means.xlsx", index=False)
    stats_table.to_excel(LOG_DIR / "core_oe_batch1_directional_welch_tests.xlsx", index=False)
    (LOG_DIR / "core_oe_methods_note.txt").write_text(STATS_METHODS + "\n", encoding="utf-8")


# %%
def main():
    raw = load_raw_data(SRC)
    data = analysis_data(raw)
    summary = summary_table(data)
    stats_table = run_stats(data)
    write_outputs(raw, data, summary, stats_table)
    plot_batch1(data, summary, stats_table)
    print(f"Analysis rows: {len(data)}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")
    print(stats_table[["plant_type", "ev_n", "core_oe_n", "p_raw", "p_holm_batch_family", "stars"]].to_string(index=False))
    return raw, data, summary, stats_table


if __name__ == "__main__":
    main()
