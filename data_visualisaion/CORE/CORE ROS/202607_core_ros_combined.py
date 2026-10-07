# %% [markdown]
# # CORE ROS combined timecourse
#
# Combine ROS plate-reader timecourses across selected experiments using the
# legacy well mappings, AUC-based replicate outlier screening, and shared figure
# styling.

# %%
from pathlib import Path
import math
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
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
else:
    HERE = Path.cwd() / "CORE" / "CORE ROS"

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))

from figure_styles import (  # noqa: E402
    ANALYSIS_RCPARAMS,
    FONT_STACK,
    PALETTE_SEQUENCE,
    TEXT_COLOR,
    apply_axis_style,
)

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_core_ros_combined"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SOURCE_ROOT = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/CSPR/ROS"
)
FOREVER_YOUNG_SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/10. Ethylene/ForeverYoung/"
    "ROSx3_core46_28_FY copy.xlsx"
)

SAVE_FMTS = ["png", "pdf", "svg"]
AUC_MAD_Z_THRESHOLD = 1.5
TIME_GRID = np.linspace(0, 1, 31)
MAIN_GENOTYPE_ORDER = ["WT", "CORE 1-2", "CORE 2-8", "CORE 704-6"]
COHORT_CONTROL = "WT (714-1 set)"
COHORT_MUTANT = "CORE 714-1"
GENOTYPE_ORDER = MAIN_GENOTYPE_ORDER + [COHORT_CONTROL, COHORT_MUTANT]
ELICITOR_ORDER = ["water", "flg22", "csp22"]
EXCLUDED_GENOTYPES = {"CORE 705-1"}
Y_LIMITS = {"water": (0, 15000), "flg22": (0, 100000), "csp22": (0, 15000)}
Y_SCALE = 1000
Y_SCALE_LABEL = "x1000"
LINE_COLORS = {
    "WT": PALETTE_SEQUENCE[2],
    "CORE 1-2": PALETTE_SEQUENCE[3],
    "CORE 2-8": PALETTE_SEQUENCE[4],
    "CORE 704-6": "#287FB8",
    COHORT_CONTROL: PALETTE_SEQUENCE[2],
    "CORE 714-1": PALETTE_SEQUENCE[6],
}

EXPERIMENTS = [
    {
        "experiment_id": "20240801_core_714_1_repeat",
        "path": Path(
            "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/"
            "CORE/GFP-COVA-VirB/EvanMutant/6wk/CORE/ROS/20240801/CORE_6wk_repeat.xlsx"
        ),
        "sheet_name": "Result sheet",
        "rows": ["A", "B", "C", "D", "E", "F", "G", "H"],
        "num_replicates": 8,
        "include_genotypes": [COHORT_CONTROL, COHORT_MUTANT],
        "sample_names": [
            "WT1 water",
            "WT1 flg22",
            "WT1 csp22",
            "WT2 water",
            "WT2 flg22",
            "WT2 csp22",
            "CORE1 water",
            "CORE1 flg22",
            "CORE1 csp22",
            "CORE2 water",
            "CORE2 flg22",
            "CORE2 csp22",
        ],
    },
    {
        "experiment_id": "20260320_2_ros2",
        "path": SOURCE_ROOT / "20260320_2" / "ros2_wt_21_7051_7064.xlsx",
        "sheet_name": 0,
        "rows": ["E", "F", "G", "H"],
        "num_replicates": 4,
        "include_genotypes": ["CORE 1-2"],
        "sample_names": [
            "WT1 water",
            "WT1 flg22",
            "WT1 csp22",
            "CORE 2-1 water",
            "CORE 2-1 flg22",
            "CORE 2-1 csp22",
            "CORE 705-1 water",
            "CORE 705-1 flg22",
            "CORE 705-1 csp22",
            "CSPR3-1-1 water",
            "CSPR3-1-1 flg22",
            "CSPR3-1-1 csp22",
        ],
    },
    {
        "experiment_id": "20260506_fy_sheet1",
        "path": FOREVER_YOUNG_SRC,
        "sheet_name": "Result sheet (1)",
        "rows": ["A", "B", "C", "D"],
        "num_replicates": 4,
        "include_genotypes": ["WT", "CORE 704-6", "CORE 2-8"],
        "sample_names": [
            "WT water",
            "WT flg22",
            "WT csp22",
            "CORE 704-6 water",
            "CORE 704-6 flg22",
            "CORE 704-6 csp22",
            "CORE 2-8 water",
            "CORE 2-8 flg22",
            "CORE 2-8 csp22",
            "FY water",
            "FY flg22",
            "FY csp22",
        ],
    },
    {
        "experiment_id": "20260506_fy_sheet2",
        "path": FOREVER_YOUNG_SRC,
        "sheet_name": "Result sheet (2)",
        "rows": ["A", "B", "C", "D"],
        "num_replicates": 4,
        "include_genotypes": ["WT", "CORE 704-6", "CORE 2-8"],
        "sample_names": [
            "WT water",
            "WT flg22",
            "WT csp22",
            "CORE 704-6 water",
            "CORE 704-6 flg22",
            "CORE 704-6 csp22",
            "CORE 2-8 water",
            "CORE 2-8 flg22",
            "CORE 2-8 csp22",
            "FY water",
            "FY flg22",
            "FY csp22",
        ],
    },
    {
        "experiment_id": "20260506_fy_sheet3",
        "path": FOREVER_YOUNG_SRC,
        "sheet_name": "Result sheet (3)",
        "rows": ["E", "F", "G", "H"],
        "num_replicates": 4,
        "include_genotypes": ["WT", "CORE 704-6", "CORE 2-8"],
        "sample_names": [
            "WT water",
            "WT flg22",
            "WT csp22",
            "CORE 704-6 water",
            "CORE 704-6 flg22",
            "CORE 704-6 csp22",
            "CORE 2-8 water",
            "CORE 2-8 flg22",
            "CORE 2-8 csp22",
            "FY water",
            "FY flg22",
            "FY csp22",
        ],
    },
]

STATS_METHODS = (
    "Stats/processing: ROS timecourse data were read from selected plate-reader "
    "workbooks using explicit sample-to-well mappings in this script. CORE 1-2 was "
    "taken only from the same-day 20260320_2 ROS set; WT, CORE 704-6, and CORE 2-8 "
    "were taken from the ForeverYoung ROSx3_core46_28_FY workbook. CORE 714-1 was taken "
    "from the validated 2024-08-01 CORE_6wk_repeat workbook (the CORE2 plate-map group) "
    "and is displayed with its same-day WT2 control as a separate paired set; "
    "the earlier attempt-2 workbook was not included. Sample names were "
    "canonicalised before plotting: WT1 is plotted as WT, and CORE 2-1 is treated as "
    "CORE 1-2. CORE 705-1, CSPR lines, and FY columns were excluded from the combined "
    "CORE figure. Within each experiment and sample, "
    "replicate wells were screened using the previous AUC-based MAD-z approach "
    f"(threshold = {AUC_MAD_Z_THRESHOLD}); flagged wells were removed before calculating "
    "experiment-level mean curves. Each experiment-level mean curve was interpolated to a "
    "common 0-1 h grid. The final figure shows the mean of experiment-level curves with "
    "a shaded SEM band across experiments, so experiments rather than wells are treated "
    "as the independent unit when multiple experiments are available. Where only one "
    "source experiment is available for a genotype/elicitor, the shaded band shows SEM "
    "across cleaned wells from that experiment for visual reference only. No inferential p-value annotations are drawn on this "
    "timecourse figure, and individual time-point markers are omitted as a figure-specific "
    "exception; cleaned well-level data, outlier logs, experiment-level curves, "
    "and AUC summaries are saved in the logs folder. Y-axis ranges are shared within "
    "each elicitor column: water 0-15,000; flg22 0-100,000; csp22 0-15,000 luminescence units. "
    "The plotted y-axis is divided by 1,000 and labelled as AU x1000."
)

ROS_RCPARAMS = {
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 15,
    "axes.titlesize": 15.5,
    "axes.labelsize": 15,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14.5,
    "legend.title_fontsize": 14.5,
    "figure.dpi": 120,
    "savefig.dpi": 300,
}
plt.rcParams.update(ROS_RCPARAMS)


# %%
def generate_sample_mapping(rows, sample_names, num_replicates):
    sample_mapping = {}
    for i, sample in enumerate(sample_names):
        sample_mapping[sample] = [f"{rows[rep]}{i + 1}" for rep in range(num_replicates)]
    return sample_mapping


def canonical_sample(sample_name):
    text = str(sample_name).strip()
    lowered = text.lower()
    elicitor = None
    for candidate in ELICITOR_ORDER:
        if lowered.endswith(candidate):
            elicitor = candidate
            genotype = text[: -len(candidate)].strip()
            break
    if elicitor is None:
        return None, None

    genotype = re.sub(r"\s+", " ", genotype)
    genotype = genotype.replace("WT2", COHORT_CONTROL)
    genotype = genotype.replace("WT1", "WT")
    genotype = genotype.replace("CORE 2-1", "CORE 1-2")
    genotype = genotype.replace("CORE2", "CORE 714-1")
    genotype = genotype.replace("706-4", "CORE 704-6")
    genotype = genotype.replace("CORE 706-4", "CORE 704-6")
    genotype = genotype.replace("CORE 704-6", "CORE 704-6")
    if genotype.startswith("CSPR"):
        return None, None
    return genotype, elicitor


def read_plate_data(path, sheet_name=0):
    data = pd.read_excel(path, sheet_name=sheet_name, skiprows=52)
    data = data.loc[:, ~data.columns.astype(str).str.contains("^Unnamed")]
    data = data.dropna(how="all").dropna(axis=1, how="all").copy()
    if "Time [s]" not in data.columns:
        raise ValueError(f"No Time [s] column after skiprows=52: {path}")
    data["time_h"] = pd.to_numeric(data["Time [s]"], errors="coerce") / 3600.0
    data = data.dropna(subset=["time_h"])
    for column in ["Cycle Nr.", "Temp. [°C]"]:
        if column in data.columns:
            data = data.drop(columns=column)
    return data


def replicate_auc(values, times):
    mask = np.isfinite(values) & np.isfinite(times)
    if mask.sum() < 2:
        return np.nan
    return float(np.trapezoid(values[mask], x=times[mask]))


def mad_z_flags(aucs, threshold=AUC_MAD_Z_THRESHOLD):
    aucs = np.asarray(aucs, dtype=float)
    median_auc = np.nanmedian(aucs)
    mad = np.nanmedian(np.abs(aucs - median_auc))
    if not np.isfinite(mad) or mad == 0:
        return np.zeros(len(aucs), dtype=bool), np.zeros(len(aucs), dtype=float), median_auc, mad
    mad_scaled = mad * 1.4826
    z_scores = np.abs((aucs - median_auc) / mad_scaled)
    return z_scores > threshold, z_scores, median_auc, mad


def parse_experiment(experiment):
    data = read_plate_data(experiment["path"], experiment.get("sheet_name", 0))
    times = data["time_h"].to_numpy(dtype=float)
    mapping = generate_sample_mapping(experiment["rows"], experiment["sample_names"], experiment["num_replicates"])
    include_genotypes = set(experiment.get("include_genotypes", GENOTYPE_ORDER))
    well_rows = []
    outlier_rows = []
    experiment_curve_rows = []

    for sample, wells in mapping.items():
        genotype, elicitor = canonical_sample(sample)
        if genotype is None or genotype in EXCLUDED_GENOTYPES or genotype not in GENOTYPE_ORDER:
            continue
        if genotype not in include_genotypes:
            continue
        present_wells = [well for well in wells if well in data.columns]
        if not present_wells:
            continue

        well_values = []
        aucs = []
        for well in present_wells:
            values = pd.to_numeric(data[well], errors="coerce").to_numpy(dtype=float)
            well_values.append(values)
            aucs.append(replicate_auc(values, times))

        flags, z_scores, median_auc, mad = mad_z_flags(aucs)
        clean_values = []
        for well, values, auc, is_outlier, z_score in zip(present_wells, well_values, aucs, flags, z_scores):
            outlier_rows.append(
                {
                    "experiment_id": experiment["experiment_id"],
                    "source_file": str(experiment["path"]),
                    "source_sheet": experiment.get("sheet_name", 0),
                    "sample_raw": sample,
                    "genotype": genotype,
                    "elicitor": elicitor,
                    "well": well,
                    "auc": auc,
                    "median_auc": median_auc,
                    "mad_auc": mad,
                    "mad_z": z_score,
                    "outlier": bool(is_outlier),
                    "threshold": AUC_MAD_Z_THRESHOLD,
                }
            )
            for time_h, value in zip(times, values):
                well_rows.append(
                    {
                        "experiment_id": experiment["experiment_id"],
                        "source_file": str(experiment["path"]),
                        "source_sheet": experiment.get("sheet_name", 0),
                        "sample_raw": sample,
                        "genotype": genotype,
                        "elicitor": elicitor,
                        "well": well,
                        "time_h": time_h,
                        "luminescence": value,
                        "auc": auc,
                        "outlier": bool(is_outlier),
                    }
                )
            if not is_outlier:
                clean_values.append(values)

        if not clean_values:
            continue
        stack = np.vstack(clean_values)
        mean_curve = np.nanmean(stack, axis=0)
        sem_curve = np.nanstd(stack, axis=0, ddof=1) / math.sqrt(stack.shape[0]) if stack.shape[0] > 1 else np.zeros(stack.shape[1])
        interp_curve = np.interp(TIME_GRID, times, mean_curve)
        interp_sem_curve = np.interp(TIME_GRID, times, sem_curve)
        for time_h, value, sem_value in zip(TIME_GRID, interp_curve, interp_sem_curve):
            experiment_curve_rows.append(
                {
                    "experiment_id": experiment["experiment_id"],
                    "source_file": str(experiment["path"]),
                    "source_sheet": experiment.get("sheet_name", 0),
                    "genotype": genotype,
                    "elicitor": elicitor,
                    "time_h": time_h,
                    "experiment_mean_luminescence": value,
                    "experiment_sem_luminescence": sem_value,
                    "n_wells_clean": len(clean_values),
                    "n_wells_original": len(present_wells),
                }
            )

    return pd.DataFrame(well_rows), pd.DataFrame(outlier_rows), pd.DataFrame(experiment_curve_rows)


def combine_experiments():
    all_wells = []
    all_outliers = []
    all_curves = []
    for experiment in EXPERIMENTS:
        wells, outliers, curves = parse_experiment(experiment)
        all_wells.append(wells)
        all_outliers.append(outliers)
        all_curves.append(curves)
    well_data = pd.concat(all_wells, ignore_index=True) if all_wells else pd.DataFrame()
    outlier_log = pd.concat(all_outliers, ignore_index=True) if all_outliers else pd.DataFrame()
    experiment_curves = pd.concat(all_curves, ignore_index=True) if all_curves else pd.DataFrame()
    return well_data, outlier_log, experiment_curves


def sem(values):
    values = pd.Series(values).dropna()
    if len(values) <= 1:
        return 0.0
    return float(values.std(ddof=1) / math.sqrt(len(values)))


def combined_curve_summary(experiment_curves):
    rows = []
    for (genotype, elicitor, time_h), subset in experiment_curves.groupby(["genotype", "elicitor", "time_h"], sort=False):
        values = subset["experiment_mean_luminescence"].dropna()
        n_experiments = values.size
        sem_luminescence = sem(values) if n_experiments > 1 else float(subset["experiment_sem_luminescence"].iloc[0])
        rows.append(
            {
                "genotype": genotype,
                "elicitor": elicitor,
                "time_h": time_h,
                "n_experiments": n_experiments,
                "mean_luminescence": values.mean(),
                "sem_luminescence": sem_luminescence,
                "sd_luminescence": values.std(ddof=1) if n_experiments > 1 else 0.0,
                "sem_basis": "experiment_means" if n_experiments > 1 else "cleaned_wells_single_experiment",
            }
        )
    return pd.DataFrame(rows)


def auc_summary(experiment_curves):
    rows = []
    for (experiment_id, genotype, elicitor), subset in experiment_curves.groupby(["experiment_id", "genotype", "elicitor"]):
        values = subset.sort_values("time_h")
        auc = replicate_auc(
            values["experiment_mean_luminescence"].to_numpy(dtype=float),
            values["time_h"].to_numpy(dtype=float),
        )
        rows.append(
            {
                "experiment_id": experiment_id,
                "genotype": genotype,
                "elicitor": elicitor,
                "auc_0_1h": auc,
                "max_luminescence_0_1h": values["experiment_mean_luminescence"].max(),
                "n_wells_clean": values["n_wells_clean"].iloc[0],
                "n_wells_original": values["n_wells_original"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.0, pad=3, direction="out")
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.88)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def rounded_upper_limit(value):
    if not np.isfinite(value) or value <= 0:
        return 1
    magnitude = 10 ** math.floor(math.log10(value))
    return math.ceil(value / magnitude) * magnitude


def autoscale_column_limits(summary):
    limits = {}
    for elicitor, subset in summary.groupby("elicitor"):
        upper = (subset["mean_luminescence"] + subset["sem_luminescence"]).max()
        limits[elicitor] = (0, rounded_upper_limit(upper * 1.08))
    return limits


def scaled_limits(y_limits):
    return {key: (low / Y_SCALE, high / Y_SCALE) for key, (low, high) in y_limits.items()}


def display_genotype(genotype):
    """Use lowercase italic labels for CORE mutants while keeping WT upright."""
    if genotype.startswith("CORE "):
        suffix = genotype.removeprefix("CORE ").replace("-", r"\!-\!")
        return rf"$\mathit{{core\!-\!{suffix}}}$"
    return "WT"


def plot_combined(summary, y_limits, output_name):
    cohort_limits = {}
    cohort_data = summary[summary["genotype"].isin([COHORT_CONTROL, COHORT_MUTANT])]
    for elicitor in ELICITOR_ORDER:
        subset = cohort_data[cohort_data["elicitor"].eq(elicitor)]
        upper = (subset["mean_luminescence"] + subset["sem_luminescence"]).max()
        cohort_limits[elicitor] = (0, rounded_upper_limit(upper * 1.08))
    # Keep water and csp22 directly comparable within the second experiment,
    # and use the requested 0--80 display range for flg22 (AU x10).
    water_csp_upper = max(cohort_limits["water"][1], cohort_limits["csp22"][1])
    cohort_limits["water"] = (0, water_csp_upper)
    cohort_limits["csp22"] = (0, water_csp_upper)
    cohort_limits["flg22"] = (0, 800)

    fig, axes = plt.subplots(
        len(MAIN_GENOTYPE_ORDER) + 2,
        len(ELICITOR_ORDER),
        figsize=(11.5, 12.6),
        sharex=True,
        constrained_layout=False,
    )

    for row_idx, genotype in enumerate(MAIN_GENOTYPE_ORDER):
        color = LINE_COLORS[genotype]
        for col_idx, elicitor in enumerate(ELICITOR_ORDER):
            ax = axes[row_idx, col_idx]
            subset = summary.loc[summary["genotype"].eq(genotype) & summary["elicitor"].eq(elicitor)].sort_values("time_h")
            if subset.empty:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes, color=TEXT_COLOR)
            else:
                x = subset["time_h"].to_numpy(dtype=float)
                y = subset["mean_luminescence"].to_numpy(dtype=float) / Y_SCALE
                err = subset["sem_luminescence"].to_numpy(dtype=float) / Y_SCALE
                ax.fill_between(x, y - err, y + err, color=color, alpha=0.20, linewidth=0, zorder=1)
                ax.plot(x, y, color=color, linewidth=2.2, zorder=3)

            title = f"{display_genotype(genotype)} {elicitor}"
            ax.set_title(title, pad=4)
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(*scaled_limits(y_limits)[elicitor])
            ax.set_xticks([0, 0.5, 1.0])
            ax.set_xticklabels(["0", "0.5", "1"])
            style_axis(ax)

    cohort_scale = 10
    for cohort_offset, genotype in enumerate([COHORT_CONTROL, COHORT_MUTANT]):
        cohort_row = len(MAIN_GENOTYPE_ORDER) + cohort_offset
        for col_idx, elicitor in enumerate(ELICITOR_ORDER):
            ax = axes[cohort_row, col_idx]
            subset = summary.loc[
                summary["genotype"].eq(genotype) & summary["elicitor"].eq(elicitor)
            ].sort_values("time_h")
            if not subset.empty:
                x = subset["time_h"].to_numpy(dtype=float)
                y = subset["mean_luminescence"].to_numpy(dtype=float) / cohort_scale
                err = subset["sem_luminescence"].to_numpy(dtype=float) / cohort_scale
                color = LINE_COLORS[genotype]
                ax.fill_between(x, y - err, y + err, color=color, alpha=0.18, linewidth=0, zorder=1)
                ax.plot(x, y, color=color, linewidth=2.0, zorder=3)
            display_name = display_genotype(genotype)
            ax.set_title(f"{display_name} {elicitor}", pad=4)
            ax.set_xlim(-0.02, 1.02)
            low, high = cohort_limits[elicitor]
            ax.set_ylim(low / cohort_scale, high / cohort_scale)
            ax.set_xticks([0, 0.5, 1.0])
            ax.set_xticklabels(["0", "0.5", "1"])
            style_axis(ax)

    fig.supxlabel("Time (hours)", x=0.54, y=0.025, color=TEXT_COLOR, fontfamily=FONT_STACK)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.96, bottom=0.09, wspace=0.27, hspace=0.62)
    fig.text(0.018, 0.655, f"Luminescence (AU {Y_SCALE_LABEL})", rotation=90,
             ha="center", va="center", color=TEXT_COLOR, fontfamily=FONT_STACK)
    fig.text(0.018, 0.235, "Luminescence (AU x10)", rotation=90,
             ha="center", va="center", color=TEXT_COLOR, fontfamily=FONT_STACK)
    divider_y = (axes[len(MAIN_GENOTYPE_ORDER) - 1, 0].get_position().y0
                 + axes[len(MAIN_GENOTYPE_ORDER), 0].get_position().y1) / 2
    fig.add_artist(Line2D([0.07, 0.985], [divider_y, divider_y], transform=fig.transFigure,
                          color="#A7ADB4", linewidth=1.0, alpha=0.9))
    save(fig, output_name)
    plt.close(fig)


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{output_stem}.{file_format}", dpi=300, bbox_inches="tight")


def write_outputs(well_data, outlier_log, experiment_curves, summary, aucs, column_limits):
    well_data.to_excel(LOG_DIR / "core_ros_well_timecourses.xlsx", index=False)
    outlier_log.to_excel(LOG_DIR / "core_ros_auc_outlier_log.xlsx", index=False)
    experiment_curves.to_excel(LOG_DIR / "core_ros_experiment_mean_curves.xlsx", index=False)
    summary.to_excel(LOG_DIR / "core_ros_combined_curve_summary.xlsx", index=False)
    aucs.to_excel(LOG_DIR / "core_ros_experiment_auc_summary.xlsx", index=False)
    manifest = pd.DataFrame(
        {
            "experiment_id": experiment["experiment_id"],
            "source_file": str(experiment["path"]),
            "source_sheet": experiment.get("sheet_name", 0),
            "rows": "".join(experiment["rows"]),
            "num_replicates": experiment["num_replicates"],
            "included_genotypes": "; ".join(experiment.get("include_genotypes", GENOTYPE_ORDER)),
            "sample_names": "; ".join(experiment["sample_names"]),
        }
        for experiment in EXPERIMENTS
    )
    manifest.to_excel(LOG_DIR / "core_ros_experiment_manifest.xlsx", index=False)
    pd.DataFrame(
        {"elicitor": elicitor, "y_min": limits[0], "y_max": limits[1]}
        for elicitor, limits in column_limits.items()
    ).to_excel(LOG_DIR / "core_ros_column_axis_limits.xlsx", index=False)
    (LOG_DIR / "core_ros_methods_note.txt").write_text(STATS_METHODS + "\n", encoding="utf-8")


# %%
def main():
    well_data, outlier_log, experiment_curves = combine_experiments()
    summary = combined_curve_summary(experiment_curves)
    aucs = auc_summary(experiment_curves)
    column_limits = Y_LIMITS
    write_outputs(well_data, outlier_log, experiment_curves, summary, aucs, column_limits)
    plot_combined(summary, column_limits, f"{SCRIPT_STEM}_timecourse")
    print(f"Well timecourse rows: {len(well_data)}")
    print(f"Experiment mean curve rows: {len(experiment_curves)}")
    print(f"Flagged outlier wells: {int(outlier_log['outlier'].sum()) if not outlier_log.empty else 0}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")
    print(
        aucs.groupby(["genotype", "elicitor"])
        .agg(n_experiments=("experiment_id", "nunique"), mean_auc=("auc_0_1h", "mean"))
        .reset_index()
        .to_string(index=False)
    )
    return well_data, outlier_log, experiment_curves, summary, aucs


if __name__ == "__main__":
    main()
