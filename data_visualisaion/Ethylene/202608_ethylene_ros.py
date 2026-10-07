"""Plot three independent ROS time-course repeats in the thesis house style.

The workbook contains one repeat per worksheet. Four plate wells represent each
genotype-by-elicitor condition. As in the source analysis notebook, wells are
screened within condition using AUC MAD-z (cut-off 1.5), and CORE-704-6 is not
included in the final figures.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator


ROOT = Path(__file__).resolve().parent
SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/Ethylene/ROS_FY/ROSx3_core46_28_FY.xlsx"
)
OUT = ROOT / "202608_ethylene_ros_figures"
LOGS = ROOT / "202608_ethylene_ros_logs"

PANEL_BG = "#E9EFF6"
GRID = "#FFFFFF"
TEXT = "#4A4A4A"
EDGE = "#000000"
COLORS = {"WT": "#4E79A7", "core-2-8": "#F28E2B", "FY": "#59A14F"}
FONT = "Helvetica Neue"
ELICITORS = ("water", "flg22", "csp22")
GENOTYPES = ("WT", "core-2-8", "FY")
SOURCE_NAMES = {"WT": "WT1", "core-2-8": "CORE-2-8", "FY": "FY"}
AUC_Z_THRESHOLD = 1.5

plt.rcParams.update(
    {
        "font.family": FONT,
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.labelweight": "bold",
        "axes.edgecolor": EDGE,
        "axes.linewidth": 0.8,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "axes.labelcolor": TEXT,
        "text.color": TEXT,
        "savefig.dpi": 300,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def read_repeat(sheet_index: int) -> pd.DataFrame:
    data = pd.read_excel(SOURCE, sheet_name=sheet_index, skiprows=52)
    data = data.loc[:, ~data.columns.astype(str).str.startswith("Unnamed")]
    data = data.dropna(how="all").dropna(axis=1, how="all")
    data["Time [h]"] = pd.to_numeric(data["Time [s]"], errors="coerce") / 3600
    return data


def mapping_for(data: pd.DataFrame) -> dict[tuple[str, str], list[str]]:
    first_well = next(col for col in data.columns if str(col)[0] in "ABCDEFGH" and str(col)[1:].isdigit())
    row_letters = "ABCD" if str(first_well)[0] == "A" else "EFGH"
    mapping = {}
    for genotype in GENOTYPES:
        source_genotype = SOURCE_NAMES[genotype]
        start = {"WT1": 1, "CORE-2-8": 7, "FY": 10}[source_genotype]
        for offset, elicitor in enumerate(ELICITORS):
            mapping[(genotype, elicitor)] = [f"{row}{start + offset}" for row in row_letters]
    return mapping


def auc(values: np.ndarray, time: np.ndarray) -> float:
    valid = np.isfinite(values) & np.isfinite(time)
    return float(np.trapezoid(values[valid], x=time[valid]))


def clean_wells(
    data: pd.DataFrame, mapping: dict[tuple[str, str], list[str]], repeat: int
) -> tuple[dict[tuple[str, str], list[str]], list[dict[str, object]]]:
    time = data["Time [h]"].to_numpy(float)
    cleaned = {}
    audit = []
    for condition, wells in mapping.items():
        aucs = np.array([auc(data[well].to_numpy(float), time) for well in wells])
        median = float(np.median(aucs))
        mad = float(np.median(np.abs(aucs - median)))
        scaled_mad = mad * 1.4826 if mad > 0 else 1e-9
        scores = (aucs - median) / scaled_mad
        cleaned[condition] = []
        for well, well_auc, score in zip(wells, aucs, scores):
            excluded = abs(score) > AUC_Z_THRESHOLD
            if not excluded:
                cleaned[condition].append(well)
            audit.append(
                {
                    "repeat": repeat,
                    "genotype": condition[0],
                    "elicitor": condition[1],
                    "well": well,
                    "auc_rlu_h": well_auc,
                    "median_auc_rlu_h": median,
                    "auc_mad_z": score,
                    "excluded": excluded,
                }
            )
    return cleaned, audit


def nice_upper(maximum: float) -> float:
    """Return a clean upper bound with about five major tick intervals."""
    if not np.isfinite(maximum) or maximum <= 0:
        return 1
    raw_step = maximum / 4.6
    magnitude = 10 ** np.floor(np.log10(raw_step))
    step = next(value for value in (1, 2, 2.5, 5, 10) if value * magnitude >= raw_step)
    return float(np.ceil(maximum / (step * magnitude)) * step * magnitude)


def axis_style(ax) -> None:
    ax.set_facecolor(PANEL_BG)
    ax.set_axisbelow(True)
    ax.grid(axis="both", color=GRID, linewidth=1.1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=0.8)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


def plot_repeat(
    data: pd.DataFrame,
    cleaned: dict[tuple[str, str], list[str]],
    repeat: int,
) -> list[dict[str, object]]:
    time = data["Time [h]"].to_numpy(float)
    summaries = []
    curves = {}
    for condition, wells in cleaned.items():
        values = data[wells].to_numpy(float).T / 1000
        curves[condition] = (values.mean(axis=0), values.std(axis=0, ddof=0))
        for well, curve in zip(wells, values):
            for time_h, value in zip(time, curve):
                summaries.append(
                    {
                        "repeat": repeat,
                        "genotype": condition[0],
                        "elicitor": condition[1],
                        "well": well,
                        "time_h": time_h,
                        "rlu_x1000": value,
                    }
                )

    # Each elicitor is a panel. Water deliberately inherits flg22's scale;
    # csp22 retains its own scale, as requested.
    panel_upper = {}
    for elicitor in ELICITORS:
        extent = max(
            np.nanmax(curves[(genotype, elicitor)][0] + curves[(genotype, elicitor)][1])
            for genotype in GENOTYPES
        )
        panel_upper[elicitor] = nice_upper(extent * 1.04)
    panel_upper["water"] = panel_upper["flg22"]

    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.65), sharex=True, facecolor="white")
    for ax, elicitor in zip(axes, ELICITORS):
        axis_style(ax)
        for genotype in GENOTYPES:
            mean, sd = curves[(genotype, elicitor)]
            color = COLORS[genotype]
            ax.fill_between(time, mean - sd, mean + sd, color=color, alpha=0.16, linewidth=0)
            ax.plot(time, mean, color=color, linewidth=2.1, label=genotype)
        ax.set_title(elicitor, fontweight="bold", pad=8)
        ax.set_ylim(0, panel_upper[elicitor])
        ax.set_xlim(time.min(), time.max())
        ax.set_xlabel("Time (h)")
    axes[0].set_ylabel("Luminescence (RLU x1000)")
    legend = axes[1].legend(
        frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.13),
        ncol=3, handlelength=1.8, columnspacing=1.3,
    )
    legend.get_texts()[1].set_fontstyle("italic")
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.17, top=0.76, wspace=0.27)

    stem = f"202608_ethylene_ROS_repeat{repeat}"
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return summaries


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(f"Source workbook not found: {SOURCE}")
    LOGS.mkdir(parents=True, exist_ok=True)
    all_audit = []
    all_clean = []
    sheet_names = pd.ExcelFile(SOURCE).sheet_names
    for sheet_index, sheet_name in enumerate(sheet_names):
        repeat = sheet_index + 1
        data = read_repeat(sheet_index)
        mapping = mapping_for(data)
        cleaned, audit = clean_wells(data, mapping, repeat)
        all_audit.extend(audit)
        all_clean.extend(plot_repeat(data, cleaned, repeat))

    pd.DataFrame(all_audit).to_csv(LOGS / "ros_auc_outlier_audit.csv", index=False)
    pd.DataFrame(all_clean).to_csv(LOGS / "ros_cleaned_long.csv", index=False)
    excluded = [row for row in all_audit if row["excluded"]]
    note = [
        "Ethylene ROS time-course processing note",
        "========================================",
        f"Source: {SOURCE}",
        f"Worksheets/repeats: {', '.join(sheet_names)}",
        "Plate mapping: four wells per condition; columns encode genotype and elicitor.",
        "Final plot includes WT, core-2-8 and FY; CORE-704-6 was omitted to match the prior subset figures.",
        f"Outlier rule: within-condition trapezoidal AUC MAD-z > {AUC_Z_THRESHOLD}; mean +/- population SD after exclusion.",
        "Scaling: RLU divided by 1000; within each repeat water uses the flg22 y-axis limits, while csp22 is independently scaled.",
        "No inferential statistical tests were applied to these time courses.",
        "",
        "Excluded wells:",
    ]
    note.extend(
        f"Repeat {row['repeat']}: {row['genotype']} {row['elicitor']} {row['well']} (MAD-z {row['auc_mad_z']:.3f})"
        for row in excluded
    )
    (LOGS / "ros_methods_note.txt").write_text("\n".join(note) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
