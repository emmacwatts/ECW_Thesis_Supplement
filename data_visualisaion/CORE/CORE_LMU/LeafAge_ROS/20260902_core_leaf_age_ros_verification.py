"""Plot two verification repeats of the CORE ROS leaf-position experiment.

The plots intentionally mirror the established 202607 analysis but remain
separate from the manuscript composite.
"""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/ROS_byleafAge/20260902/ros_x2_leafAge.xlsx"
)
BASE_SCRIPT = HERE / "202607_core_leaf_age_ros.py"
SHEETS = ["Sheet6", "Sheet5"]
OUTPUT_ROOT = HERE / "20260902_core_leaf_age_ros_verification"

spec = spec_from_file_location("core_leaf_age_ros_base", BASE_SCRIPT)
base = module_from_spec(spec)
spec.loader.exec_module(base)

base.LEAVES = list(range(1, 9))
base.LEAF_LAYOUT = {
    1: (list("ABCD"), 1), 2: (list("ABCD"), 4),
    3: (list("ABCD"), 7), 4: (list("ABCD"), 10),
    5: (list("EFGH"), 1), 6: (list("EFGH"), 4),
    7: (list("EFGH"), 7), 8: (list("EFGH"), 10),
}
base.LEAF_COLORS = dict(zip(base.LEAVES, [
    "#F0D43A", "#F5B335", "#F18A3D", "#E45F52",
    "#D13F6A", "#AF327F", "#812B86", "#54246F",
]))


def read_transposed_plate(sheet):
    """Return one Tecan sheet as time rows and well columns."""
    raw = pd.read_excel(SOURCE, sheet_name=sheet, header=None)
    time_matches = raw.index[raw.iloc[:, 0].astype(str).eq("Time [s]")]
    if len(time_matches) != 1:
        raise ValueError(f"Expected one Time [s] row in {sheet}; found {len(time_matches)}")
    time_row = int(time_matches[0])
    times = pd.to_numeric(raw.iloc[time_row, 1:], errors="coerce")
    values = {"Time [s]": times.to_numpy(float)}
    for index in range(time_row + 2, len(raw)):
        well = str(raw.iloc[index, 0]).strip()
        if not well or well == "nan":
            continue
        if well[0] in "ABCDEFGH" and well[1:].isdigit():
            values[well] = pd.to_numeric(raw.iloc[index, 1:], errors="coerce").to_numpy(float)
    plate = pd.DataFrame(values).dropna(subset=["Time [s]"])
    expected = {f"{row}{column}" for row in "ABCDEFGH" for column in range(1, 13)}
    missing = sorted(expected.difference(plate.columns))
    if missing:
        raise ValueError(f"Missing wells in {sheet}: {', '.join(missing)}")
    plate["time_h"] = plate["Time [s]"] / 3600
    return plate


def plate_column(leaf, elicitor, sheet):
    column = base.LEAF_LAYOUT[leaf][1] + base.ELICITOR_OFFSET[elicitor]
    if sheet == "Sheet5" and column in (11, 12):
        column = 23 - column
    return column


def process_plate(sheet):
    plate = read_transposed_plate(sheet)
    source_times = plate["time_h"].to_numpy(float)
    if source_times.min() > 0 or source_times.max() < 1:
        raise ValueError(f"{sheet} does not span the full 0-1 h interval")
    well_rows, outlier_rows, clean_rows = [], [], []
    for leaf in base.LEAVES:
        rows, _ = base.LEAF_LAYOUT[leaf]
        for elicitor in base.ELICITORS:
            column = plate_column(leaf, elicitor, sheet)
            wells, curves, aucs = [], [], []
            for row in rows:
                well = f"{row}{column}"
                source_values = pd.to_numeric(plate[well], errors="coerce").to_numpy(float)
                curve = np.interp(base.TIME_GRID, source_times, source_values)
                auc = base.curve_auc(base.TIME_GRID, curve)
                wells.append(well)
                curves.append(curve)
                aucs.append(auc)
                for time_h, value in zip(base.TIME_GRID, curve):
                    well_rows.append({"sheet": sheet, "leaf": leaf, "elicitor": elicitor,
                                      "well": well, "time_h": time_h, "luminescence": value})
            scores = base.mad_z(aucs)
            flagged = np.abs(scores) > base.AUC_MAD_Z_THRESHOLD
            if flagged.all():
                flagged[:] = False
            for well, auc, score, flag, curve in zip(wells, aucs, scores, flagged, curves):
                outlier_rows.append({"sheet": sheet, "leaf": leaf, "elicitor": elicitor,
                                     "well": well, "auc_0_1h": auc, "mad_z": score,
                                     "flagged": bool(flag), "threshold": base.AUC_MAD_Z_THRESHOLD})
                if not flag:
                    for time_h, value in zip(base.TIME_GRID, curve):
                        clean_rows.append({"sheet": sheet, "leaf": leaf, "elicitor": elicitor,
                                           "well": well, "time_h": time_h, "luminescence": value,
                                           "auc_0_1h": auc})
    return pd.DataFrame(well_rows), pd.DataFrame(outlier_rows), pd.DataFrame(clean_rows)


def draw_auc(ax, aucs, summary, elicitor, y_max=None):
    raw = aucs[aucs["elicitor"].eq(elicitor)]
    stats_df = summary[summary["elicitor"].eq(elicitor)].set_index("leaf")
    x = np.arange(len(base.LEAVES))
    means = stats_df.loc[base.LEAVES, "mean"].to_numpy() / base.Y_SCALE
    errors = stats_df.loc[base.LEAVES, "sem"].to_numpy() / base.Y_SCALE
    ax.plot(x, means, color="#425D78", linewidth=1.5)
    ax.errorbar(x, means, yerr=errors, fmt="o", color="#425D78", ecolor="#333333",
                markersize=5, elinewidth=1, capsize=3, zorder=4)
    for index, leaf in enumerate(base.LEAVES):
        values = raw.loc[raw["leaf"].eq(leaf), "auc_0_1h"].to_numpy() / base.Y_SCALE
        ax.scatter(index + np.linspace(-0.09, 0.09, len(values)), values, s=28,
                   color=base.LEAF_COLORS[leaf], edgecolor="none", alpha=0.8, zorder=3)
    ax.set_xticks(x, base.LEAVES)
    ax.set(xlabel="Leaf position", ylabel="Integrated ROS response\n(AU·h x1000)",
           title=f"{base.ELICITOR_LABELS[elicitor]} AUC (0–1 h)", ylim=(0, None))
    if y_max is not None:
        ax.set_ylim(0, y_max)
    base.style_axis(ax)


def make_figures(curves, aucs, auc_summary):
    """Use the established layouts, expanding the legend for eight leaves."""
    fig, axes = base.plt.subplots(1, 2, figsize=(14.8, 5.8),
                                 gridspec_kw={"width_ratios": [1.42, 1]})
    base.draw_curves(axes[0], curves, "csp22", legend=True)
    draw_auc(axes[1], aucs, auc_summary, "csp22")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.88, bottom=0.16, wspace=0.20)
    base.save(fig, "csp22_leaf_position")
    base.plt.close(fig)

    fig, axes = base.plt.subplots(2, 3, figsize=(17.2, 9.4))
    for column, elicitor in enumerate(base.ELICITORS):
        base.draw_curves(axes[0, column], curves, elicitor, legend=False)
        draw_auc(axes[1, column], aucs, auc_summary, elicitor)
        if column:
            axes[0, column].set_ylabel("")
            axes[1, column].set_ylabel("")
    # Water and flg22 retain a common scale, but let that scale follow each
    # repeat's observed curves (including SEM) so the flg22 peak stays visible.
    water_flg_upper = max(axes[0, 0].get_ylim()[1], axes[0, 1].get_ylim()[1])
    axes[0, 0].set_ylim(0, water_flg_upper)
    axes[0, 1].set_ylim(0, water_flg_upper)
    flg_auc_upper = axes[1, 1].get_ylim()[1]
    axes[1, 0].set_ylim(0, flg_auc_upper)
    axes[1, 1].set_ylim(0, flg_auc_upper)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=len(base.LEAVES),
               loc="upper center", bbox_to_anchor=(0.53, 1.015),
               columnspacing=1.05, handlelength=1.55)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.82, bottom=0.09,
                        wspace=0.15, hspace=0.58)
    base.save(fig, "all_elicitors_leaf_position")
    base.plt.close(fig)


def run_sheet(sheet):
    output = OUTPUT_ROOT / sheet
    base.SCRIPT_STEM = f"20260902_core_leaf_age_ros_{sheet.lower()}"
    base.SAVE_DIR = output / "figures"
    base.LOG_DIR = output / "logs"
    base.SAVE_DIR.mkdir(parents=True, exist_ok=True)
    base.LOG_DIR.mkdir(parents=True, exist_ok=True)

    wells, outliers, clean = process_plate(sheet)
    curves, aucs, auc_summary = base.summaries(clean)
    anova, pairs = base.statistics(aucs)
    tables = {
        "raw_well_timecourses": wells, "auc_outlier_log": outliers,
        "flagged_outliers_only": outliers.loc[outliers["flagged"]].copy(),
        "clean_well_timecourses": clean, "curve_summary": curves,
        "clean_well_auc_values": aucs, "auc_summary": auc_summary,
        "one_way_anova": anova, "pairwise_leaf_tests": pairs,
    }
    for name, table in tables.items():
        table.to_excel(base.LOG_DIR / f"{base.SCRIPT_STEM}_{name}.xlsx", index=False)
    note = (
        f"Verification repeat from {SOURCE}, sheet {sheet}. Leaves 1-4 occupy rows A-D and "
        "leaves 5-8 rows E-H; leaf blocks run left-to-right in columns 1-3, 4-6, 7-9, "
        "and 10-12. Each leaf x treatment has four technical replicate wells. Within each "
        "block the order is water, flg22, csp22. For Sheet5 only, physical columns 11 and "
        "12 were swapped during plate setup and are corrected during import. Curves were "
        "linearly interpolated to include an endpoint at exactly 1 h and AUC was integrated "
        "over exactly 0-1 h. Outlier handling, summaries, statistics, and styling match the "
        "202607 first-replicate analysis. These plots are verification-only.\n"
    )
    (base.LOG_DIR / f"{base.SCRIPT_STEM}_methods_note.txt").write_text(note, encoding="utf-8")
    make_figures(curves, aucs, auc_summary)
    print(f"{sheet}: flagged {int(outliers['flagged'].sum())}/{len(outliers)} wells; {output}")


def main():
    for sheet in SHEETS:
        run_sheet(sheet)


if __name__ == "__main__":
    main()
