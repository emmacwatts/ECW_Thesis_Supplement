# %% [markdown]
# # CORE ROS response across plant age
#
# Combine ROS plate-reader runs collected every two days from 28 to 40 days.
# The plate layout contains three biological plants, each with water, flg22,
# and csp22 columns and eight technical wells per plant-condition.

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
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd() / "CORE" / "ROS_timecourse"
STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))

from figure_styles import (  # noqa: E402
    ANALYSIS_RCPARAMS,
    FONT_STACK,
    TEXT_COLOR,
    apply_axis_style,
)

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_core_ros_age_timecourse"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SOURCE_ROOT = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/ROS_Timecourse/Att3_Worked"
)
EXPERIMENTS = [
    (28, "20250624_28d", "ROS_20250624_28d.xlsx", "Result sheet"),
    (30, "20250626_30d", "ROS_30d.xlsx", "Result sheet (1)"),
    (32, "20250628_32d", "ROS_timecourse.xlsx", "Result sheet (1)"),
    (34, "20250630_34d", "ROS_timecourse.xlsx", "Result sheet (1)"),
    (36, "20250702_36d", "ROSTimeCourse_andoagroLiveBoiledLast3Lanes.xlsx", "Result sheet (1)"),
    (38, "20250704_38d", "time_andAgroaLaB_Long.xlsx", "Result sheet"),
    (40, "20250706_40d", "40d.xlsx", "Result sheet"),
]
AGES = [item[0] for item in EXPERIMENTS]
ELICITOR_ORDER = ["water", "flg22", "csp22"]
ELICITOR_LABELS = {"water": "Water", "flg22": "flg22", "csp22": "csp22"}
ROWS = list("ABCDEFGH")
PLANTS = [1, 2, 3]
WELL_COLUMNS = {
    1: {"water": 1, "flg22": 2, "csp22": 3},
    2: {"water": 4, "flg22": 5, "csp22": 6},
    3: {"water": 7, "flg22": 8, "csp22": 9},
}
TIME_GRID = np.linspace(0, 1, 31)
AUC_MAD_Z_THRESHOLD = 1.5
Y_SCALE = 1000
SAVE_FMTS = ["png", "pdf", "svg"]

# Match the high-contrast ordered palette used in the companion leaf-position
# ROS timecourse. Seven distinct hues remain after excluding the failed 42 d run.
AGE_COLORS = dict(zip(AGES, [
    "#F0D43A",  # yellow
    "#F5B335",  # golden yellow
    "#F18A3D",  # orange
    "#E45F52",  # coral-orange
    "#D13F6A",  # raspberry
    "#AF327F",  # magenta
    "#812B86",  # deep purple-magenta
]))

METHODS = (
    "Stats/processing: ROS plate-reader workbooks collected from 28 to 40 days were "
    "read using the plate layout recorded in the original analysis notebooks. Columns "
    "1-3, 4-6, and 7-9 represent biological plants 1, 2, and 3, respectively; within "
    "each plant the columns are water, flg22, and csp22. Rows A-H are technical leaf-disc "
    "wells. Technical wells were screened within each age x plant x elicitor group using "
    "AUC-based robust MAD-z scores (threshold 1.5), matching the established CORE ROS "
    "pipeline. Clean technical wells were averaged to one curve per biological plant. "
    "Plant curves were interpolated to a common 0-1 h grid, and plotted means +/- SEM "
    "therefore use plants as the independent unit (n=3). Plant-level AUC was calculated "
    "over 0-1 h. For each elicitor, a one-factor repeated-measures ANOVA tested age with "
    "plant as subject; paired two-sided age contrasts were Holm-adjusted within elicitor. "
    "The 42-day experiment was excluded because the assay did not work reliably. "
    "No baseline subtraction or transformation was applied. The main figure focuses on "
    "csp22, consistent with the prior daily exports; all elicitors are retained in logs "
    "and a supplementary overview. Luminescence is displayed as AU x1000."
)

plt.rcParams.update(
    {
        **ANALYSIS_RCPARAMS,
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": 16,
        "axes.titlesize": 18,
        "axes.labelsize": 17,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "legend.fontsize": 13,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    }
)


def read_plate(path, sheet_name):
    data = pd.read_excel(path, sheet_name=sheet_name, skiprows=52)
    data = data.loc[:, ~data.columns.astype(str).str.contains("^Unnamed")]
    data = data.dropna(how="all").dropna(axis=1, how="all").copy()
    if "Time [s]" not in data.columns:
        raise ValueError(f"No Time [s] column after skiprows=52: {path} / {sheet_name}")
    data["time_h"] = pd.to_numeric(data["Time [s]"], errors="coerce") / 3600
    data = data.dropna(subset=["time_h"])
    return data


def curve_auc(times, values):
    mask = np.isfinite(times) & np.isfinite(values)
    return float(np.trapezoid(values[mask], x=times[mask])) if mask.sum() >= 2 else np.nan


def mad_z(aucs):
    values = np.asarray(aucs, dtype=float)
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    if not np.isfinite(mad) or mad == 0:
        return np.zeros(len(values))
    return 0.6745 * (values - median) / mad


def process_experiments():
    well_rows = []
    outlier_rows = []
    plant_curve_rows = []
    plant_auc_rows = []

    for age, folder, filename, sheet_name in EXPERIMENTS:
        path = SOURCE_ROOT / folder / filename
        plate = read_plate(path, sheet_name)
        times = plate["time_h"].to_numpy(dtype=float)
        use_time = times <= 1.0001

        for plant in PLANTS:
            for elicitor in ELICITOR_ORDER:
                column_number = WELL_COLUMNS[plant][elicitor]
                wells = [f"{row}{column_number}" for row in ROWS]
                curves = []
                aucs = []
                valid_wells = []
                for well in wells:
                    if well not in plate.columns:
                        continue
                    values = pd.to_numeric(plate[well], errors="coerce").to_numpy(dtype=float)[use_time]
                    well_times = times[use_time]
                    auc = curve_auc(well_times, values)
                    curves.append(values)
                    aucs.append(auc)
                    valid_wells.append(well)
                    for time_h, value in zip(well_times, values):
                        well_rows.append(
                            {"age_days": age, "plant": plant, "elicitor": elicitor,
                             "well": well, "time_h": time_h, "luminescence": value,
                             "source_file": str(path), "source_sheet": sheet_name}
                        )

                z_scores = mad_z(aucs)
                flagged = np.abs(z_scores) > AUC_MAD_Z_THRESHOLD
                for well, auc, score, is_flagged in zip(valid_wells, aucs, z_scores, flagged):
                    outlier_rows.append(
                        {"age_days": age, "plant": plant, "elicitor": elicitor,
                         "well": well, "auc_0_1h": auc, "mad_z": score,
                         "flagged": bool(is_flagged), "threshold": AUC_MAD_Z_THRESHOLD}
                    )

                clean_curves = [curve for curve, is_flagged in zip(curves, flagged) if not is_flagged]
                if not clean_curves:
                    clean_curves = curves
                mean_curve = np.nanmean(np.vstack(clean_curves), axis=0)
                interpolated = np.interp(TIME_GRID, times[use_time], mean_curve)
                plant_auc = curve_auc(TIME_GRID, interpolated)
                plant_auc_rows.append(
                    {"age_days": age, "plant": plant, "elicitor": elicitor,
                     "auc_0_1h": plant_auc, "technical_wells_total": len(curves),
                     "technical_wells_clean": len(clean_curves)}
                )
                for time_h, value in zip(TIME_GRID, interpolated):
                    plant_curve_rows.append(
                        {"age_days": age, "plant": plant, "elicitor": elicitor,
                         "time_h": time_h, "luminescence": value,
                         "technical_wells_clean": len(clean_curves)}
                    )

    return (
        pd.DataFrame(well_rows),
        pd.DataFrame(outlier_rows),
        pd.DataFrame(plant_curve_rows),
        pd.DataFrame(plant_auc_rows),
    )


def curve_summary(plant_curves):
    return (
        plant_curves.groupby(["age_days", "elicitor", "time_h"])
        .agg(n_plants=("plant", "nunique"), mean_luminescence=("luminescence", "mean"),
             sem_luminescence=("luminescence", "sem"))
        .reset_index()
    )


def auc_summary(plant_aucs):
    return (
        plant_aucs.groupby(["age_days", "elicitor"])
        .agg(n_plants=("plant", "nunique"), mean_auc=("auc_0_1h", "mean"),
             sem_auc=("auc_0_1h", "sem"), sd_auc=("auc_0_1h", "std"))
        .reset_index()
    )


def repeated_measures_anova(plant_aucs):
    rows = []
    for elicitor in ELICITOR_ORDER:
        subset = plant_aucs.loc[plant_aucs["elicitor"].eq(elicitor)].copy()
        result = AnovaRM(subset, depvar="auc_0_1h", subject="plant", within=["age_days"]).fit()
        table = result.anova_table.reset_index().rename(columns={"index": "term"})
        table.insert(0, "elicitor", elicitor)
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def paired_age_tests(plant_aucs):
    rows = []
    for elicitor in ELICITOR_ORDER:
        subset = plant_aucs.loc[plant_aucs["elicitor"].eq(elicitor)]
        wide = subset.pivot(index="plant", columns="age_days", values="auc_0_1h")
        for age_a, age_b in combinations(AGES, 2):
            paired = wide[[age_a, age_b]].dropna()
            result = stats.ttest_rel(paired[age_a], paired[age_b]) if len(paired) >= 2 else None
            rows.append(
                {"elicitor": elicitor, "age_a": age_a, "age_b": age_b, "paired_n": len(paired),
                 "mean_b_minus_a": (paired[age_b] - paired[age_a]).mean() if len(paired) else np.nan,
                 "t_statistic": float(result.statistic) if result else np.nan,
                 "p_raw": float(result.pvalue) if result else np.nan,
                 "test": "paired two-sided t-test"}
            )
    table = pd.DataFrame(rows)
    table["p_holm_within_elicitor"] = np.nan
    for elicitor, indexes in table.groupby("elicitor").groups.items():
        valid = table.loc[indexes].index[table.loc[indexes, "p_raw"].notna()]
        table.loc[valid, "p_holm_within_elicitor"] = multipletests(table.loc[valid, "p_raw"], method="holm")[1]
    return table


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def draw_curves(ax, curve_stats, elicitor, show_legend=True):
    subset = curve_stats.loc[curve_stats["elicitor"].eq(elicitor)]
    for age in AGES:
        age_data = subset.loc[subset["age_days"].eq(age)]
        x = age_data["time_h"].to_numpy()
        mean = age_data["mean_luminescence"].to_numpy() / Y_SCALE
        error = age_data["sem_luminescence"].fillna(0).to_numpy() / Y_SCALE
        color = AGE_COLORS[age]
        ax.plot(x, mean, color=color, linewidth=1.8, label=f"{age} d")
        ax.fill_between(x, mean - error, mean + error, color=color, alpha=0.12, linewidth=0)
    ax.set_xlabel("Time after elicitation (h)")
    ax.set_ylabel("Luminescence (AU x1000)")
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.set_title(ELICITOR_LABELS[elicitor])
    ax.xaxis.set_major_locator(MaxNLocator(6))
    style_axis(ax)
    if show_legend:
        ax.legend(frameon=False, ncol=2, title="Plant age", loc="upper right")


def draw_auc(ax, plant_aucs, auc_stats, elicitor):
    raw = plant_aucs.loc[plant_aucs["elicitor"].eq(elicitor)]
    summary = auc_stats.loc[auc_stats["elicitor"].eq(elicitor)].set_index("age_days")
    x = np.arange(len(AGES))
    means = summary.loc[AGES, "mean_auc"].to_numpy() / Y_SCALE
    errors = summary.loc[AGES, "sem_auc"].to_numpy() / Y_SCALE
    ax.plot(x, means, color="#425D78", linewidth=1.5, zorder=2)
    ax.errorbar(x, means, yerr=errors, fmt="o", color="#425D78", ecolor="#333333",
                markersize=5, elinewidth=1.0, capsize=3, zorder=4)
    for index, age in enumerate(AGES):
        values = raw.loc[raw["age_days"].eq(age), "auc_0_1h"].to_numpy() / Y_SCALE
        ax.scatter(index + np.linspace(-0.08, 0.08, len(values)), values, s=27,
                   color=AGE_COLORS[age], edgecolor="none", alpha=0.78, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(AGES)
    ax.set_xlabel("Plant age (days)")
    ax.set_ylabel("Integrated ROS response\n(AU·h x1000)")
    ax.set_ylim(bottom=0)
    ax.set_title(f"{ELICITOR_LABELS[elicitor]} AUC (0–1 h)")
    style_axis(ax)


def save(fig, name):
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{SCRIPT_STEM}_{name}.{file_format}", dpi=300, bbox_inches="tight")


def plot_main(curve_stats, plant_aucs, auc_stats):
    fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.8), gridspec_kw={"width_ratios": [1.42, 1.0]})
    draw_curves(axes[0], curve_stats, "csp22")
    draw_auc(axes[1], plant_aucs, auc_stats, "csp22")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.88, bottom=0.16, wspace=0.20)
    save(fig, "csp22_age_timecourse")
    plt.close(fig)


def plot_all_elicitors(curve_stats, plant_aucs, auc_stats):
    fig, axes = plt.subplots(2, 3, figsize=(17.2, 9.4))
    for column, elicitor in enumerate(ELICITOR_ORDER):
        draw_curves(axes[0, column], curve_stats, elicitor, show_legend=(column == 2))
        draw_auc(axes[1, column], plant_aucs, auc_stats, elicitor)
        if column > 0:
            axes[0, column].set_ylabel("")
            axes[1, column].set_ylabel("")
    # Match Water directly to flg22 for both the timecourse and AUC rows;
    # csp22 retains its independent autoscaling.
    flg_curve_limits = axes[0, 1].get_ylim()
    flg_auc_limits = axes[1, 1].get_ylim()
    axes[0, 0].set_ylim(flg_curve_limits)
    axes[1, 0].set_ylim(flg_auc_limits)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.94, bottom=0.09, wspace=0.15, hspace=0.30)
    save(fig, "all_elicitors_age_timecourse")
    plt.close(fig)


def write_outputs(wells, outliers, plant_curves, plant_aucs, curve_stats, auc_stats, anova, pairs):
    wells.to_excel(LOG_DIR / "ros_age_raw_well_timecourses.xlsx", index=False)
    outliers.to_excel(LOG_DIR / "ros_age_auc_outlier_log.xlsx", index=False)
    plant_curves.to_excel(LOG_DIR / "ros_age_plant_mean_curves.xlsx", index=False)
    plant_aucs.to_excel(LOG_DIR / "ros_age_plant_auc_values.xlsx", index=False)
    curve_stats.to_excel(LOG_DIR / "ros_age_curve_summary.xlsx", index=False)
    auc_stats.to_excel(LOG_DIR / "ros_age_auc_summary.xlsx", index=False)
    anova.to_excel(LOG_DIR / "ros_age_repeated_measures_anova.xlsx", index=False)
    pairs.to_excel(LOG_DIR / "ros_age_paired_age_contrasts.xlsx", index=False)
def main():
    wells, outliers, plant_curves, plant_aucs = process_experiments()
    curve_stats = curve_summary(plant_curves)
    auc_stats = auc_summary(plant_aucs)
    anova = repeated_measures_anova(plant_aucs)
    pairs = paired_age_tests(plant_aucs)
    write_outputs(wells, outliers, plant_curves, plant_aucs, curve_stats, auc_stats, anova, pairs)
    plot_main(curve_stats, plant_aucs, auc_stats)
    plot_all_elicitors(curve_stats, plant_aucs, auc_stats)
    print(f"Ages: {AGES}")
    print(f"Plant-level AUC rows: {len(plant_aucs)}")
    print(f"Flagged technical wells: {int(outliers['flagged'].sum())} / {len(outliers)}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")
    print(anova.to_string(index=False))
    return wells, outliers, plant_curves, plant_aucs, curve_stats, auc_stats, anova, pairs


if __name__ == "__main__":
    main()
