"""Regenerate the 2026-03-25 CSPR ROS assay in the current thesis style."""

from pathlib import Path
import math
import os
import sys

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STYLE_DIR = HERE.parents[2] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CSPR/ROS/Clearest CSPR ROS 20260325/"
    "ros2_wt_1-2_64_311_wfcca.xlsx"
)
OUT = HERE / "202607_cspr_ros_timecourse_figures"
LOG = HERE / "202607_cspr_ros_timecourse_logs"
for folder in (OUT, LOG):
    folder.mkdir(parents=True, exist_ok=True)

AUC_MAD_Z_THRESHOLD = 1.5
TIME_LIMIT_H = 1.0
Y_SCALE = 1000
GROUP_ORDER = ["WT", "core-1-2", "cspr-706-4", "cspr-3-1-1"]
CONDITION_ORDER = ["water", "flg22", "csp22"]
COLORS = {
    "WT": "#F0B429",
    "core-1-2": "#E76F51",
    "cspr-706-4": "#2A9D8F",
    "cspr-3-1-1": "#31539A",
}
MAPPING = {
    "WT": {
        "water": ["A1", "B1", "C1"], "flg22": ["A2", "B2", "C2"],
        "csp22": ["A3", "B3", "C3"],
        "csp22 + Agrobacterium": ["A4", "B4", "C4"],
    },
    "core-1-2": {
        "water": ["A5", "B5", "C5"], "flg22": ["A6", "B6", "C6"],
        "csp22": ["A7", "B7", "C7"],
        "csp22 + Agrobacterium": ["A8", "B8", "C8"],
    },
    "cspr-706-4": {
        "water": ["A9", "B9", "C9"], "flg22": ["A10", "B10", "C10"],
        "csp22": ["A11", "B11", "C11"],
        "csp22 + Agrobacterium": ["A12", "B12", "C12"],
    },
    "cspr-3-1-1": {
        "water": ["D1", "E1", "F1"], "flg22": ["D2", "E2", "F2"],
        "csp22": ["D3", "E3", "F3"],
        "csp22 + Agrobacterium": ["D4", "E4", "F4"],
    },
}

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": 15, "axes.titlesize": 15.5, "axes.labelsize": 15,
    "xtick.labelsize": 13.5, "ytick.labelsize": 13.5,
    "legend.fontsize": 13.5, "figure.dpi": 120, "savefig.dpi": 300,
})


def auc(values, times):
    valid = np.isfinite(values) & np.isfinite(times)
    return float(np.trapezoid(values[valid], x=times[valid])) if valid.sum() >= 2 else np.nan


def mad_outliers(aucs):
    values = np.asarray(aucs, dtype=float)
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    if not np.isfinite(mad) or mad == 0:
        return np.zeros(len(values), dtype=bool), np.zeros(len(values)), median, mad
    z = np.abs((values - median) / (1.4826 * mad))
    return z > AUC_MAD_Z_THRESHOLD, z, median, mad


def process():
    raw = pd.read_excel(SOURCE, sheet_name="Result sheet (1)", skiprows=52)
    raw = raw.loc[:, ~raw.columns.astype(str).str.contains("^Unnamed")]
    raw = raw.dropna(how="all").dropna(axis=1, how="all")
    raw["time_h"] = pd.to_numeric(raw["Time [s]"], errors="coerce") / 3600
    raw = raw.dropna(subset=["time_h"])
    raw = raw.loc[raw["time_h"].le(TIME_LIMIT_H + 1e-8)].copy()
    times = raw["time_h"].to_numpy(float)

    well_rows, outlier_rows, curve_rows = [], [], []
    for genotype in GROUP_ORDER:
        for condition in CONDITION_ORDER:
            wells = MAPPING[genotype][condition]
            stacks = [
                pd.to_numeric(raw[well], errors="coerce").to_numpy(float)
                for well in wells
            ]
            aucs = [auc(values, times) for values in stacks]
            flags, scores, median_auc, mad = mad_outliers(aucs)
            kept = []
            for well, values, well_auc, flag, score in zip(wells, stacks, aucs, flags, scores):
                outlier_rows.append({
                    "genotype": genotype, "condition": condition, "well": well,
                    "auc_0_1h": well_auc, "median_auc": median_auc,
                    "mad_auc": mad, "absolute_mad_z": score,
                    "outlier": bool(flag), "threshold": AUC_MAD_Z_THRESHOLD,
                })
                for time_h, value in zip(times, values):
                    well_rows.append({
                        "genotype": genotype, "condition": condition, "well": well,
                        "time_h": time_h, "luminescence": value, "outlier": bool(flag),
                    })
                if not flag:
                    kept.append(values)
            stack = np.vstack(kept)
            means = np.nanmean(stack, axis=0)
            sems = (
                np.nanstd(stack, axis=0, ddof=1) / math.sqrt(len(stack))
                if len(stack) > 1 else np.zeros(len(times))
            )
            for time_h, mean, sem in zip(times, means, sems):
                curve_rows.append({
                    "genotype": genotype, "condition": condition, "time_h": time_h,
                    "mean_luminescence": mean, "sem_luminescence": sem,
                    "n_clean_wells": len(kept), "n_original_wells": len(wells),
                })
    return pd.DataFrame(well_rows), pd.DataFrame(outlier_rows), pd.DataFrame(curve_rows)


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    apply_axis_style(ax)


def round_up(value):
    magnitude = 10 ** math.floor(math.log10(value))
    return math.ceil(value / magnitude) * magnitude


def plot(curves):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.5), sharex=True)
    maxima = {}
    for condition in CONDITION_ORDER:
        subset = curves[curves["condition"].eq(condition)]
        maxima[condition] = (
            subset["mean_luminescence"] + subset["sem_luminescence"]
        ).max()
    # Water and flg22 use the same scale, matching the recent ROS figure convention.
    shared_control_max = round_up(max(maxima["water"], maxima["flg22"]) * 1.08)
    limits = {
        "water": shared_control_max, "flg22": shared_control_max,
        "csp22": round_up(maxima["csp22"] * 1.08),
    }

    handles = []
    for ax, condition in zip(axes, CONDITION_ORDER):
        for genotype in GROUP_ORDER:
            subset = curves[
                curves["condition"].eq(condition) & curves["genotype"].eq(genotype)
            ].sort_values("time_h")
            x = subset["time_h"].to_numpy(float)
            y = subset["mean_luminescence"].to_numpy(float) / Y_SCALE
            error = subset["sem_luminescence"].to_numpy(float) / Y_SCALE
            ax.fill_between(
                x, y - error, y + error, color=COLORS[genotype],
                alpha=0.17, linewidth=0,
            )
            line, = ax.plot(x, y, color=COLORS[genotype], linewidth=2.0)
            if ax is axes[0]:
                handles.append(line)
        ax.set_title(condition)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, limits[condition] / Y_SCALE)
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        style_axis(ax)

    axes[0].set_ylabel("Luminescence (AU ×1000)")
    fig.supxlabel("Time after elicitation (h)", y=0.06)
    legend_labels = ["WT", r"$\it{core}$-1-2", r"$\it{cspr}$-706-4", r"$\it{cspr}$-3-1-1"]
    fig.legend(
        handles, legend_labels, loc="upper center", bbox_to_anchor=(0.53, 1.01),
        frameon=False, ncol=4, handlelength=2.0, columnspacing=1.4,
    )
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.18, top=0.80, wspace=0.25)
    stem = "202607_cspr_ros_wt_core_cspr_timecourse"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return limits


def main():
    wells, outliers, curves = process()
    limits = plot(curves)
    with pd.ExcelWriter(LOG / "cspr_ros_processed_data_and_outliers.xlsx") as writer:
        wells.to_excel(writer, sheet_name="well_timecourses", index=False)
        outliers.to_excel(writer, sheet_name="auc_outlier_log", index=False)
        curves.to_excel(writer, sheet_name="cleaned_curve_summary", index=False)
        pd.DataFrame([
            {"condition": condition, "y_min": 0, "y_max": upper}
            for condition, upper in limits.items()
        ]).to_excel(writer, sheet_name="axis_limits", index=False)
worksheet of "
        f"{SOURCE}. The explicit well mapping was reproduced from "
        "revised_indiv-plateReader.ipynb in the source folder. Curves are restricted "
        "to 0-1 h. Replicate wells were screened within each genotype-condition group "
        f"using AUC-based absolute MAD z-scores (threshold {AUC_MAD_Z_THRESHOLD}); "
        "flagged wells were excluded. Lines show the mean and shaded regions show SEM "
        "across retained wells. Water and flg22 panels share a y-axis range. This is a "
        "single biological experiment, so the uncertainty bands describe within-plate "
        "well variation and no inferential significance test is shown.\n",
        encoding="utf-8",
    )
    print(outliers.groupby(["genotype", "condition"])["outlier"].agg(["sum", "count"]))
    print(f"Figures: {OUT}")


if __name__ == "__main__":
    main()
