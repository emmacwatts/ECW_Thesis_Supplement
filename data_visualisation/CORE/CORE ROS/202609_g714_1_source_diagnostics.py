"""Diagnostic plots for the two candidate G714-1 CORE ROS source files.

The plate map was recovered and verified against CORE_repeat_6wk.html:
columns 1--12 are WT1, WT2, CORE1, and CORE2, with water, flg22,
and csp22 in that order; rows A--H are replicate wells.

This script is deliberately separate from 202607_core_ros_combined.py and does
not alter the combined CORE ROS analysis.
"""

from pathlib import Path
import os
import sys

HERE = Path(__file__).resolve().parent
MPLCONFIGDIR = HERE.parents[1] / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

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

SOURCE_DIR = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/"
    "CORE/GFP-COVA-VirB/EvanMutant/6wk/CORE/ROS"
)
SOURCES = {
    "2024-08-01 repeat": SOURCE_DIR / "20240801" / "CORE_6wk_repeat.xlsx",
    "attempt 2 (2024-07-30)": SOURCE_DIR / "CORE_6wk_ROS_att2.xlsx",
}

OUT_DIR = HERE / "202609_g714_1_source_diagnostics"
OUT_DIR.mkdir(parents=True, exist_ok=True)
ELICITORS = ["water", "flg22", "csp22"]
GROUPS = ["WT1", "WT2", "CORE1", "CORE2"]
ROWS = list("ABCDEFGH")
COLORS = dict(zip(GROUPS, PALETTE_SEQUENCE[2:6]))

plt.rcParams.update(
    {
        **ANALYSIS_RCPARAMS,
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": 11,
        "axes.titlesize": 11.5,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    }
)


def read_source(label, path):
    data = pd.read_excel(path, sheet_name="Result sheet", skiprows=52)
    data = data.loc[:, ~data.columns.astype(str).str.contains("^Unnamed")]
    data = data.dropna(how="all").copy()
    data["time_h"] = pd.to_numeric(data["Time [s]"], errors="coerce") / 3600
    data = data.dropna(subset=["time_h"])

    long_rows = []
    for group_idx, group in enumerate(GROUPS):
        for elicitor_idx, elicitor in enumerate(ELICITORS):
            column = group_idx * 3 + elicitor_idx + 1
            for row in ROWS:
                well = f"{row}{column}"
                values = pd.to_numeric(data[well], errors="coerce")
                for time_h, value in zip(data["time_h"], values):
                    long_rows.append(
                        {
                            "source": label,
                            "source_file": str(path),
                            "group": group,
                            "elicitor": elicitor,
                            "well": well,
                            "time_h": time_h,
                            "luminescence": value,
                        }
                    )
    return pd.DataFrame(long_rows)


def summarise(long_data):
    return (
        long_data.groupby(["source", "group", "elicitor", "time_h"], sort=False)
        .agg(
            mean_luminescence=("luminescence", "mean"),
            sd_luminescence=("luminescence", "std"),
            sem_luminescence=("luminescence", "sem"),
            n_wells=("luminescence", "count"),
        )
        .reset_index()
    )


def axis_style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.88)
    ax.set_axisbelow(True)
    ax.tick_params(width=1, length=4, direction="out")
    apply_axis_style(ax)


def plot_source(source, long_data, summary, y_limits):
    fig, axes = plt.subplots(4, 3, figsize=(11.7, 10.3), sharex=True)
    for row_idx, group in enumerate(GROUPS):
        color = COLORS[group]
        for col_idx, elicitor in enumerate(ELICITORS):
            ax = axes[row_idx, col_idx]
            wells = long_data.loc[
                long_data["group"].eq(group) & long_data["elicitor"].eq(elicitor)
            ]
            curve = summary.loc[
                summary["group"].eq(group) & summary["elicitor"].eq(elicitor)
            ].sort_values("time_h")
            for _, well_data in wells.groupby("well", sort=False):
                well_data = well_data.sort_values("time_h")
                ax.plot(
                    well_data["time_h"],
                    well_data["luminescence"] / 1000,
                    color=color,
                    alpha=0.16,
                    linewidth=0.7,
                )
            x = curve["time_h"].to_numpy(float)
            mean = curve["mean_luminescence"].to_numpy(float) / 1000
            sem = curve["sem_luminescence"].to_numpy(float) / 1000
            ax.fill_between(x, mean - sem, mean + sem, color=color, alpha=0.22, linewidth=0)
            ax.plot(x, mean, color=color, linewidth=2.1)
            ax.set_title(f"{group} {elicitor}", pad=4)
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(0, y_limits[elicitor] / 1000)
            ax.set_xticks([0, 0.5, 1])
            axis_style(ax)

    fig.suptitle(f"Candidate G714-1 source diagnostic: {source}\nmean ± SEM; faint lines are individual wells (n=8)", y=0.995)
    fig.supylabel("Luminescence (AU x1000)", x=0.018, color=TEXT_COLOR)
    fig.supxlabel("Time (hours)", x=0.53, y=0.025, color=TEXT_COLOR)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.92, bottom=0.08, wspace=0.25, hspace=0.48)
    return fig


def main():
    frames = [read_source(label, path) for label, path in SOURCES.items()]
    long_data = pd.concat(frames, ignore_index=True)
    summary = summarise(long_data)

    figures = []
    axis_rows = []
    for source in SOURCES:
        source_data = long_data.loc[long_data["source"].eq(source)]
        y_limits = {}
        for elicitor in ELICITORS:
            upper = source_data.loc[source_data["elicitor"].eq(elicitor), "luminescence"].max()
            step = 10000 if upper >= 10000 else 1000
            y_limits[elicitor] = max(step, int(np.ceil(upper * 1.05 / step) * step))
            axis_rows.append(
                {"source": source, "elicitor": elicitor, "y_min": 0, "y_max": y_limits[elicitor]}
            )
        fig = plot_source(
            source,
            source_data,
            summary.loc[summary["source"].eq(source)],
            y_limits,
        )
        safe_name = source.lower().replace(" ", "_").replace("(", "").replace(")", "")
        fig.savefig(OUT_DIR / f"{safe_name}.png", bbox_inches="tight")
        figures.append(fig)

    with PdfPages(OUT_DIR / "g714_1_candidate_source_diagnostics.pdf") as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    for fig in figures:
        plt.close(fig)

    long_data.to_excel(OUT_DIR / "g714_1_candidate_raw_long.xlsx", index=False)
    summary.to_excel(OUT_DIR / "g714_1_candidate_curve_summary.xlsx", index=False)
    pd.DataFrame(axis_rows).to_excel(OUT_DIR / "g714_1_candidate_axis_limits.xlsx", index=False)

    peaks = (
        summary.groupby(["source", "group", "elicitor"], sort=False)["mean_luminescence"]
        .max()
        .rename("peak_mean_luminescence")
        .reset_index()
    )
    peaks.to_excel(OUT_DIR / "g714_1_candidate_peak_summary.xlsx", index=False)
    print(f"Diagnostic report: {OUT_DIR / 'g714_1_candidate_source_diagnostics.pdf'}")
    print(f"Raw rows: {len(long_data):,}; summary rows: {len(summary):,}")
    print(peaks.to_string(index=False))


if __name__ == "__main__":
    main()
