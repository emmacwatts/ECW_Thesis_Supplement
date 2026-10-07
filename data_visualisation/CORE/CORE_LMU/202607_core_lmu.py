# %% [markdown]
# # CORE leaf-position (LMU) GFP analysis
#
# Replot the 2025-04-26 CORE lower/middle/upper leaf experiment with leaf
# position as the primary x-axis grouping and plant genotype encoded by colour.

# %%
from pathlib import Path
from itertools import combinations
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
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd() / "CORE" / "CORE_LMU"
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

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_core_lmu"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/CombinedTest/LMU_LargeBatch_20250426/"
    "core exp/20250426_CORE.xlsx"
)
SOURCE_PRISM = SRC.with_name("lmu_prism_20250426.prism")
SHEET = "LMU"
LEAF_ORDER = ["Lower", "Middle", "Upper"]
LEAF_LABELS = ["Leaf 2", "Leaf 5", "Leaf 7"]
PLANT_ORDER = ["cas9 WT", "714-1", "704-6", "1-2"]
PLANT_LABELS = {
    "cas9 WT": "Cas9 WT",
    "714-1": "core-714-1",
    "704-6": "core-704-6",
    "1-2": "core-1-2",
}
EXCLUDED_PLANTS = {"705-1", "2-8"}
P19_ORDER = ["no p19", "p19"]
PANEL_ORDER = [(1, "no p19"), (1, "p19")]
COLORS = {
    "cas9 WT": "#79C6A3",
    "714-1": "#38A6A5",
    "704-6": "#287FB8",
    "1-2": "#31539A",
}
SAVE_FMTS = ["png", "pdf", "svg"]
Y_SCALE = 1000
BAR_ALPHA = 0.76
BAR_WIDTH = 0.16

METHODS = (
    "Stats/processing: lower/middle/upper (LMU) endpoint data were read from the LMU "
    "sheet of 20250426_CORE.xlsx. Normalised GFP (source GFP minus background) was the "
    "response. CORE 705-1 and CORE 2-8 were excluded before summaries, plots, and "
    "inference, and the analysis was restricted to batch 1. The no-p19 and p19 "
    "conditions were analysed separately. Bars show mean +/- SEM with "
    "individual plant values and n labels. For each condition, an ordinary two-factor "
    "ANOVA (Type III sums of squares) tested leaf position, plant type, and their "
    "interaction, matching the structure used in the prior Prism processing. Because "
    "the same plants contributed multiple leaf positions, paired two-sided t-tests were "
    "also retained as targeted within-plant position contrasts for each genotype; Holm "
    "adjustment was applied across the valid position contrasts within each batch x p19 "
    "condition. The primary plotted tests asked whether each CORE line was higher than "
    "Cas9 WT within the same leaf position, using one-sided Welch tests (CORE > WT) with "
    "Holm adjustment across the nine comparisons within each p19 condition. Significant "
    "increases are marked above the corresponding CORE bars; decreases were not tested. Groups "
    "with fewer than two observations were plotted but not inferentially tested. No "
    "outlier removal or transformation was applied. The figure is scaled by 1,000. "
    f"Prior processing reference: {SOURCE_PRISM}."
)

plt.rcParams.update(
    {
        **ANALYSIS_RCPARAMS,
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": 17,
        "axes.titlesize": 19,
        "axes.labelsize": 18,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "legend.fontsize": 16,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    }
)


def sem(values):
    values = pd.Series(values).dropna()
    return np.nan if len(values) < 2 else float(values.sem())


def load_raw_data(src):
    raw = pd.read_excel(src, sheet_name=SHEET).iloc[:, :8].copy()
    raw = raw.rename(
        columns={
            "Plant Type": "plant_type",
            "Leaf Age": "leaf_position",
            "Batch": "batch",
            "Plant": "plant",
            "p19": "p19",
            "Background": "background",
            "GFP": "gfp",
            "Normalised GFP": "normalized_gfp",
        }
    )
    for column in ["batch", "plant", "background", "gfp", "normalized_gfp"]:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    for column in ["plant_type", "leaf_position", "p19"]:
        raw[column] = raw[column].astype(str).str.strip()
    raw = raw.dropna(subset=["plant_type", "leaf_position", "batch", "plant", "p19", "normalized_gfp"])
    raw = raw.loc[raw["leaf_position"].isin(LEAF_ORDER) & raw["p19"].isin(P19_ORDER)].copy()
    raw["source_file"] = str(src)
    raw["source_sheet"] = SHEET
    return raw.reset_index(drop=True)


def analysis_data(raw):
    data = raw.loc[
        raw["batch"].eq(1)
        & raw["plant_type"].isin(PLANT_ORDER)
        & ~raw["plant_type"].isin(EXCLUDED_PLANTS)
    ].copy()
    data["leaf_position"] = pd.Categorical(data["leaf_position"], LEAF_ORDER, ordered=True)
    data["plant_type"] = pd.Categorical(data["plant_type"], PLANT_ORDER, ordered=True)
    data["p19"] = pd.Categorical(data["p19"], P19_ORDER, ordered=True)
    data["plant_id"] = (
        "B" + data["batch"].astype(int).astype(str)
        + "_" + data["plant_type"].astype(str)
        + "_P" + data["plant"].astype(int).astype(str)
    )
    return data.sort_values(["batch", "p19", "leaf_position", "plant_type", "plant"]).reset_index(drop=True)


def summary_table(data):
    return (
        data.groupby(["batch", "p19", "leaf_position", "plant_type"], observed=True)
        .agg(
            n=("normalized_gfp", "size"),
            mean_normalized_gfp=("normalized_gfp", "mean"),
            sem_normalized_gfp=("normalized_gfp", sem),
            sd_normalized_gfp=("normalized_gfp", "std"),
        )
        .reset_index()
    )


def two_way_anova(data):
    rows = []
    for batch, p19 in PANEL_ORDER:
        subset = data.loc[data["batch"].eq(batch) & data["p19"].eq(p19)].copy()
        subset["leaf_position"] = subset["leaf_position"].astype(str)
        subset["plant_type"] = subset["plant_type"].astype(str)
        model = ols("normalized_gfp ~ C(leaf_position) * C(plant_type)", data=subset).fit()
        table = sm.stats.anova_lm(model, typ=3).reset_index().rename(columns={"index": "term"})
        table.insert(0, "p19", p19)
        table.insert(0, "batch", batch)
        table["model"] = "ordinary two-factor ANOVA, Type III SS"
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def paired_position_tests(data):
    rows = []
    for batch, p19 in PANEL_ORDER:
        condition = data.loc[data["batch"].eq(batch) & data["p19"].eq(p19)]
        for plant_type in PLANT_ORDER:
            subset = condition.loc[condition["plant_type"].eq(plant_type)]
            wide = subset.pivot_table(index="plant_id", columns="leaf_position", values="normalized_gfp")
            for position_a, position_b in combinations(LEAF_ORDER, 2):
                paired = wide[[position_a, position_b]].dropna() if set([position_a, position_b]).issubset(wide.columns) else pd.DataFrame()
                if len(paired) >= 2:
                    result = stats.ttest_rel(paired[position_a], paired[position_b])
                    p_raw = float(result.pvalue)
                    statistic = float(result.statistic)
                    mean_difference = float((paired[position_b] - paired[position_a]).mean())
                else:
                    p_raw = statistic = mean_difference = np.nan
                rows.append(
                    {
                        "batch": batch,
                        "p19": p19,
                        "plant_type": plant_type,
                        "position_a": position_a,
                        "position_b": position_b,
                        "paired_n": len(paired),
                        "mean_b_minus_a": mean_difference,
                        "t_statistic": statistic,
                        "p_raw": p_raw,
                        "test": "paired two-sided t-test",
                    }
                )
    table = pd.DataFrame(rows)
    table["p_holm_within_condition"] = np.nan
    for _, indexes in table.groupby(["batch", "p19"]).groups.items():
        valid = table.loc[indexes, "p_raw"].notna()
        valid_indexes = table.loc[indexes].index[valid]
        if len(valid_indexes):
            table.loc[valid_indexes, "p_holm_within_condition"] = multipletests(
                table.loc[valid_indexes, "p_raw"], method="holm"
            )[1]
    return table


def genotype_vs_wt_tests(data):
    rows = []
    for batch, p19 in PANEL_ORDER:
        condition = data.loc[data["batch"].eq(batch) & data["p19"].eq(p19)]
        for position in LEAF_ORDER:
            at_position = condition.loc[condition["leaf_position"].eq(position)]
            wt = at_position.loc[at_position["plant_type"].eq("cas9 WT"), "normalized_gfp"].dropna()
            for plant_type in PLANT_ORDER[1:]:
                test_values = at_position.loc[at_position["plant_type"].eq(plant_type), "normalized_gfp"].dropna()
                if len(wt) >= 2 and len(test_values) >= 2:
                    result = stats.ttest_ind(test_values, wt, equal_var=False, alternative="greater")
                    p_raw = float(result.pvalue)
                    statistic = float(result.statistic)
                else:
                    p_raw = statistic = np.nan
                rows.append(
                    {
                        "batch": batch,
                        "p19": p19,
                        "leaf_position": position,
                        "comparison": f"{PLANT_LABELS[plant_type]} > Cas9 WT",
                        "wt_n": len(wt),
                        "test_n": len(test_values),
                        "wt_mean": wt.mean() if len(wt) else np.nan,
                        "test_mean": test_values.mean() if len(test_values) else np.nan,
                        "t_statistic": statistic,
                        "p_raw": p_raw,
                        "test": "one-sided Welch t-test, CORE > WT",
                    }
                )
    table = pd.DataFrame(rows)
    table["p_holm_within_condition"] = np.nan
    for _, indexes in table.groupby(["batch", "p19"]).groups.items():
        valid_indexes = table.loc[indexes].index[table.loc[indexes, "p_raw"].notna()]
        if len(valid_indexes):
            table.loc[valid_indexes, "p_holm_within_condition"] = multipletests(
                table.loc[valid_indexes, "p_raw"], method="holm"
            )[1]
    return table


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(width=1.0, length=4.5, pad=3, direction="out")
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def panel_positions(plant_types):
    centers = np.arange(len(LEAF_ORDER), dtype=float)
    offsets = (np.arange(len(plant_types)) - (len(plant_types) - 1) / 2) * (BAR_WIDTH + 0.025)
    return centers, offsets


def p_to_stars(p_value):
    if pd.isna(p_value) or p_value >= 0.05:
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def plot_lmu(data, genotype_tests):
    fig, axes = plt.subplots(1, 2, figsize=(15.8, 6.8), sharey=True)
    global_ymax = data["normalized_gfp"].max() / Y_SCALE

    for ax, (batch, p19) in zip(axes, PANEL_ORDER):
        subset = data.loc[data["batch"].eq(batch) & data["p19"].eq(p19)]
        present_types = [plant for plant in PLANT_ORDER if subset["plant_type"].eq(plant).any()]
        centers, offsets = panel_positions(present_types)

        for plant_index, plant_type in enumerate(present_types):
            for leaf_index, leaf_position in enumerate(LEAF_ORDER):
                values = subset.loc[
                    subset["plant_type"].eq(plant_type)
                    & subset["leaf_position"].eq(leaf_position),
                    "normalized_gfp",
                ].dropna() / Y_SCALE
                if values.empty:
                    continue
                x = centers[leaf_index] + offsets[plant_index]
                mean = values.mean()
                error = values.sem() if len(values) > 1 else 0
                color = COLORS[plant_type]
                ax.bar(x, mean, width=BAR_WIDTH, color=color, alpha=BAR_ALPHA,
                       edgecolor="black", linewidth=0.9, zorder=2)
                ax.errorbar(x, mean, yerr=error, fmt="none", ecolor=ERRORBAR_COLOR,
                            elinewidth=1.0, capsize=2.8, capthick=1.0, zorder=5)
                jitter = np.linspace(-0.025, 0.025, len(values)) if len(values) > 1 else np.array([0.0])
                ax.scatter(x + jitter, values, s=21, color=color, edgecolor="none",
                           alpha=INDIVIDUAL_POINT_ALPHA, zorder=4)
                label_y = max(float(values.max()), mean + error) + 0.55
                ax.text(x, label_y, str(len(values)), ha="center", va="bottom",
                        fontsize=12.5, color="#878787")
                if plant_type != "cas9 WT":
                    test_row = genotype_tests.loc[
                        genotype_tests["batch"].eq(batch)
                        & genotype_tests["p19"].eq(p19)
                        & genotype_tests["leaf_position"].eq(leaf_position)
                        & genotype_tests["comparison"].eq(f"{PLANT_LABELS[plant_type]} > Cas9 WT")
                    ]
                    if not test_row.empty:
                        stars = p_to_stars(test_row["p_holm_within_condition"].iloc[0])
                        if stars:
                            ax.text(x, label_y + 1.0, stars, ha="center", va="bottom",
                                    fontsize=17, color=TEXT_COLOR, fontweight="bold")

        ax.set_xticks(centers)
        ax.set_xticklabels(LEAF_LABELS)
        ax.set_title("No p19" if p19 == "no p19" else "p19", pad=10)
        ax.set_ylim(0, max(30, np.ceil((global_ymax + 2) / 5) * 5))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        style_axis(ax)

    axes[0].set_ylabel("Relative GFP intensity (AU x1000)")
    for ax in axes:
        ax.set_xlabel("Leaf position")

    handles = [Patch(facecolor=COLORS[p], edgecolor="black", alpha=BAR_ALPHA, label=PLANT_LABELS[p]) for p in PLANT_ORDER]
    legend = fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.52, 1.005),
                        ncol=len(handles), frameon=False, handlelength=1.1, columnspacing=1.35)
    for text, plant in zip(legend.get_texts(), PLANT_ORDER):
        if plant != "cas9 WT":
            text.set_fontstyle("italic")

    fig.subplots_adjust(left=0.09, right=0.985, top=0.84, bottom=0.15, wspace=0.10)
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{SCRIPT_STEM}_leaf_position_grouped.{file_format}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_outputs(raw, data, summary, anova, paired, genotype_tests):
    raw.to_excel(LOG_DIR / "core_lmu_raw_long.xlsx", index=False)
    data.to_excel(LOG_DIR / "core_lmu_analysis_long.xlsx", index=False)
    summary.to_excel(LOG_DIR / "core_lmu_summary_means.xlsx", index=False)
    anova.to_excel(LOG_DIR / "core_lmu_two_way_anova.xlsx", index=False)
    paired.to_excel(LOG_DIR / "core_lmu_paired_leaf_position_tests.xlsx", index=False)
    genotype_tests.to_excel(LOG_DIR / "core_lmu_genotype_vs_wt_tests.xlsx", index=False)
def main():
    raw = load_raw_data(SRC)
    data = analysis_data(raw)
    summary = summary_table(data)
    anova = two_way_anova(data)
    paired = paired_position_tests(data)
    genotype_tests = genotype_vs_wt_tests(data)
    write_outputs(raw, data, summary, anova, paired, genotype_tests)
    plot_lmu(data, genotype_tests)
    print(f"Analysis rows: {len(data)}")
    print(f"Excluded 705-1 rows: {raw['plant_type'].eq('705-1').sum()}")
    print(f"Excluded 2-8 rows: {raw['plant_type'].eq('2-8').sum()}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")
    return raw, data, summary, anova, paired, genotype_tests


if __name__ == "__main__":
    main()
