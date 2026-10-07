# %% [markdown]
# # CORE ROS response by leaf age / position

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
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
STYLE_DIR = HERE.parents[2] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style  # noqa: E402

SCRIPT_STEM = Path(__file__).stem
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/ROS_byleafAge/20260324/leaf1-7-ros.xlsx"
)
PRIOR_PROCESSING = SOURCE.with_name("leaf1-7ROS.html")
SHEET = "Result sheet (1)"
LEAVES = list(range(1, 8))
ELICITORS = ["water", "flg22", "csp22"]
ELICITOR_LABELS = {"water": "Water", "flg22": "flg22", "csp22": "csp22"}
TIME_GRID = np.linspace(0, 1, 31)
AUC_MAD_Z_THRESHOLD = 1.5
Y_SCALE = 1000
SAVE_FMTS = ["png", "pdf", "svg"]

# Explicit layout reconstructed from the prior processing output.
# Each tuple is (row block, first column); the three successive columns are
# water, flg22 and csp22. Rows are replicate leaf discs.
LEAF_LAYOUT = {
    1: (list("ABCD"), 1), 2: (list("ABCD"), 4),
    3: (list("ABCD"), 7), 4: (list("ABCD"), 10),
    5: (list("EFGH"), 1), 6: (list("EFGH"), 4),
    7: (list("EFGH"), 7),
}
ELICITOR_OFFSET = {"water": 0, "flg22": 1, "csp22": 2}
# A warm yellow-to-magenta sequence keeps leaf position visibly ordered while
# remaining separate from the grey, blue, green and muted salmon used for the
# genotypes in the companion LMU plot.
LEAF_COLORS = dict(zip(LEAVES, [
    "#F0D43A",  # yellow
    "#F5B335",  # golden yellow
    "#F18A3D",  # orange
    "#E45F52",  # coral-orange
    "#D13F6A",  # raspberry
    "#AF327F",  # magenta
    "#812B86",  # deep purple-magenta
]))

METHODS = (
    "ROS data were read from leaf1-7-ros.xlsx using the plate map documented by the "
    "previous leaf1-7ROS.html processing. Leaves 1-4 use replicate rows A-D in column "
    "blocks 1-3, 4-6, 7-9 and 10-12; leaves 5-7 use replicate rows E-H in blocks "
    "1-3, 4-6 and 7-9. Within each block the columns are water, flg22 and csp22. "
    "Replicate leaf-disc wells were screened within leaf x elicitor using AUC-based "
    "robust MAD-z scores (threshold 1.5), matching the established CORE ROS pipeline. "
    "Curves show cleaned-well means +/- SEM; points in AUC panels are individual cleaned "
    "leaf-disc wells. AUC was integrated over 0-1 h. Because the source layout does not "
    "identify matched plants across all seven leaf positions, inference is unpaired: an "
    "ordinary one-way ANOVA tests leaf position within each elicitor, followed by all "
    "pairwise two-sided Welch tests with Holm adjustment within elicitor. No baseline "
    "subtraction or transformation was applied. Luminescence is displayed as AU x1000. "
    f"Prior processing reference: {PRIOR_PROCESSING}."
)

plt.rcParams.update({
    **ANALYSIS_RCPARAMS, "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    # This six-axis panel is reduced more strongly than the companion LMU
    # panel in the composite; larger source typography equalises final size.
    "font.size": 24, "axes.titlesize": 27, "axes.labelsize": 26,
    "xtick.labelsize": 24, "ytick.labelsize": 24, "legend.fontsize": 22,
    "figure.dpi": 120, "savefig.dpi": 300,
})


def read_plate():
    data = pd.read_excel(SOURCE, sheet_name=SHEET, skiprows=52)
    data = data.loc[:, ~data.columns.astype(str).str.contains("^Unnamed")]
    data = data.dropna(how="all").dropna(axis=1, how="all").copy()
    data["time_h"] = pd.to_numeric(data["Time [s]"], errors="coerce") / 3600
    return data.dropna(subset=["time_h"])


def curve_auc(times, values):
    mask = np.isfinite(times) & np.isfinite(values)
    return float(np.trapezoid(values[mask], x=times[mask])) if mask.sum() >= 2 else np.nan


def mad_z(values):
    values = np.asarray(values, dtype=float)
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    return np.zeros(len(values)) if not np.isfinite(mad) or mad == 0 else 0.6745 * (values - median) / mad


def process_plate():
    plate = read_plate()
    times = plate["time_h"].to_numpy(float)
    use = times <= 1.0001
    times = times[use]
    well_rows, outlier_rows, clean_rows = [], [], []
    for leaf in LEAVES:
        rows, start_column = LEAF_LAYOUT[leaf]
        for elicitor in ELICITORS:
            column = start_column + ELICITOR_OFFSET[elicitor]
            wells, curves, aucs = [], [], []
            for row in rows:
                well = f"{row}{column}"
                values = pd.to_numeric(plate[well], errors="coerce").to_numpy(float)[use]
                auc = curve_auc(times, values)
                wells.append(well); curves.append(values); aucs.append(auc)
                for time_h, value in zip(times, values):
                    well_rows.append({"leaf": leaf, "elicitor": elicitor, "well": well,
                                      "time_h": time_h, "luminescence": value})
            scores = mad_z(aucs)
            flagged = np.abs(scores) > AUC_MAD_Z_THRESHOLD
            if flagged.all():
                flagged[:] = False
            for well, auc, score, flag, curve in zip(wells, aucs, scores, flagged, curves):
                outlier_rows.append({"leaf": leaf, "elicitor": elicitor, "well": well,
                                     "auc_0_1h": auc, "mad_z": score, "flagged": bool(flag),
                                     "threshold": AUC_MAD_Z_THRESHOLD})
                if not flag:
                    interpolated = np.interp(TIME_GRID, times, curve)
                    for time_h, value in zip(TIME_GRID, interpolated):
                        clean_rows.append({"leaf": leaf, "elicitor": elicitor, "well": well,
                                           "time_h": time_h, "luminescence": value,
                                           "auc_0_1h": curve_auc(TIME_GRID, interpolated)})
    return pd.DataFrame(well_rows), pd.DataFrame(outlier_rows), pd.DataFrame(clean_rows)


def summaries(clean):
    curves = (clean.groupby(["leaf", "elicitor", "time_h"])
              .agg(n=("well", "nunique"), mean=("luminescence", "mean"), sem=("luminescence", "sem"))
              .reset_index())
    aucs = clean[["leaf", "elicitor", "well", "auc_0_1h"]].drop_duplicates()
    auc_summary = (aucs.groupby(["leaf", "elicitor"])
                   .agg(n=("well", "nunique"), mean=("auc_0_1h", "mean"),
                        sem=("auc_0_1h", "sem"), sd=("auc_0_1h", "std"))
                   .reset_index())
    return curves, aucs, auc_summary


def statistics(aucs):
    anova_rows, pair_rows = [], []
    for elicitor in ELICITORS:
        subset = aucs[aucs["elicitor"].eq(elicitor)]
        groups = [subset.loc[subset["leaf"].eq(leaf), "auc_0_1h"].to_numpy() for leaf in LEAVES]
        result = stats.f_oneway(*groups)
        anova_rows.append({"elicitor": elicitor, "test": "ordinary one-way ANOVA",
                           "factor": "leaf", "df_between": len(groups) - 1,
                           "df_within": sum(map(len, groups)) - len(groups),
                           "F": result.statistic, "p_value": result.pvalue})
        for leaf_a, leaf_b in combinations(LEAVES, 2):
            a, b = groups[leaf_a - 1], groups[leaf_b - 1]
            test = stats.ttest_ind(a, b, equal_var=False)
            pair_rows.append({"elicitor": elicitor, "leaf_a": leaf_a, "leaf_b": leaf_b,
                              "n_a": len(a), "n_b": len(b), "mean_b_minus_a": b.mean() - a.mean(),
                              "t_statistic": test.statistic, "p_raw": test.pvalue,
                              "test": "two-sided Welch t-test"})
    pairs = pd.DataFrame(pair_rows)
    pairs["p_holm_within_elicitor"] = np.nan
    for _, indexes in pairs.groupby("elicitor").groups.items():
        pairs.loc[indexes, "p_holm_within_elicitor"] = multipletests(pairs.loc[indexes, "p_raw"], method="holm")[1]
    return pd.DataFrame(anova_rows), pairs


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def draw_curves(ax, curve_summary, elicitor, legend=False, y_max=None):
    subset = curve_summary[curve_summary["elicitor"].eq(elicitor)]
    for leaf in LEAVES:
        data = subset[subset["leaf"].eq(leaf)]
        x = data["time_h"].to_numpy(); mean = data["mean"].to_numpy() / Y_SCALE
        sem = data["sem"].fillna(0).to_numpy() / Y_SCALE
        ax.plot(x, mean, color=LEAF_COLORS[leaf], linewidth=1.8, label=f"Leaf {leaf}")
        ax.fill_between(x, mean - sem, mean + sem, color=LEAF_COLORS[leaf], alpha=0.12, linewidth=0)
    ax.set(xlabel="Time after elicitation (h)", ylabel="Luminescence (AU x1000)",
           title=ELICITOR_LABELS[elicitor], xlim=(0, 1), ylim=(0, None))
    if y_max is not None:
        ax.set_ylim(0, y_max)
    ax.xaxis.set_major_locator(MaxNLocator(6)); style_axis(ax)
    if legend:
        ax.legend(frameon=False, ncol=2, title="Leaf position", loc="upper right")


def draw_auc(ax, aucs, summary, elicitor, y_max=None):
    raw = aucs[aucs["elicitor"].eq(elicitor)]
    stats_df = summary[summary["elicitor"].eq(elicitor)].set_index("leaf")
    x = np.arange(7); means = stats_df.loc[LEAVES, "mean"].to_numpy() / Y_SCALE
    errors = stats_df.loc[LEAVES, "sem"].to_numpy() / Y_SCALE
    ax.plot(x, means, color="#425D78", linewidth=1.5)
    ax.errorbar(x, means, yerr=errors, fmt="o", color="#425D78", ecolor="#333333",
                markersize=5, elinewidth=1, capsize=3, zorder=4)
    for index, leaf in enumerate(LEAVES):
        values = raw.loc[raw["leaf"].eq(leaf), "auc_0_1h"].to_numpy() / Y_SCALE
        ax.scatter(index + np.linspace(-0.09, 0.09, len(values)), values, s=28,
                   color=LEAF_COLORS[leaf], edgecolor="none", alpha=0.8, zorder=3)
    ax.set_xticks(x, LEAVES)
    ax.set(xlabel="Leaf position", ylabel="Integrated ROS response\n(AU·h x1000)",
           title=f"{ELICITOR_LABELS[elicitor]} AUC (0–1 h)", ylim=(0, None))
    if y_max is not None:
        ax.set_ylim(0, y_max)
    style_axis(ax)


def save(fig, suffix):
    for fmt in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{SCRIPT_STEM}_{suffix}.{fmt}", dpi=300, bbox_inches="tight")


def make_figures(curves, aucs, auc_summary):
    fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.8), gridspec_kw={"width_ratios": [1.42, 1]})
    draw_curves(axes[0], curves, "csp22", legend=True)
    draw_auc(axes[1], aucs, auc_summary, "csp22")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.88, bottom=0.16, wspace=0.20)
    save(fig, "csp22_leaf_position"); plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(17.2, 9.4))
    for column, elicitor in enumerate(ELICITORS):
        matched_curve_max = 140 if elicitor in {"water", "flg22"} else None
        draw_curves(
            axes[0, column], curves, elicitor, legend=False,
            y_max=matched_curve_max,
        )
        draw_auc(axes[1, column], aucs, auc_summary, elicitor)
        if column:
            axes[0, column].set_ylabel(""); axes[1, column].set_ylabel("")
    # Use the naturally autoscaled flg22 AUC range for water as well. This
    # preserves the complete flg22 data/error bars while keeping those two
    # panels directly comparable; csp22 remains independently autoscaled.
    flg_auc_upper = axes[1, 1].get_ylim()[1]
    axes[1, 0].set_ylim(0, flg_auc_upper)
    axes[1, 1].set_ylim(0, flg_auc_upper)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=7,
               loc="upper center", bbox_to_anchor=(0.53, 1.015),
               columnspacing=1.15, handlelength=1.7)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.82, bottom=0.09, wspace=0.15, hspace=0.58)
    save(fig, "all_elicitors_leaf_position"); plt.close(fig)


def main():
    wells, outliers, clean = process_plate()
    curves, aucs, auc_summary = summaries(clean)
    anova, pairs = statistics(aucs)
    wells.to_excel(LOG_DIR / "leaf_age_ros_raw_well_timecourses.xlsx", index=False)
    outliers.to_excel(LOG_DIR / "leaf_age_ros_auc_outlier_log.xlsx", index=False)
    clean.to_excel(LOG_DIR / "leaf_age_ros_clean_well_timecourses.xlsx", index=False)
    curves.to_excel(LOG_DIR / "leaf_age_ros_curve_summary.xlsx", index=False)
    aucs.to_excel(LOG_DIR / "leaf_age_ros_clean_well_auc_values.xlsx", index=False)
    auc_summary.to_excel(LOG_DIR / "leaf_age_ros_auc_summary.xlsx", index=False)
    anova.to_excel(LOG_DIR / "leaf_age_ros_one_way_anova.xlsx", index=False)
    pairs.to_excel(LOG_DIR / "leaf_age_ros_pairwise_leaf_tests.xlsx", index=False)
make_figures(curves, aucs, auc_summary)
    print(f"Flagged wells: {int(outliers['flagged'].sum())} / {len(outliers)}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(anova.to_string(index=False))
    return wells, outliers, clean, curves, aucs, auc_summary, anova, pairs


if __name__ == "__main__":
    main()
