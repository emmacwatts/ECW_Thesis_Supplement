"""Plot labelled CORE ELISA datasets, using each experiment's actual date."""

from pathlib import Path
import os
import sys

HERE = Path(__file__).resolve().parent
ROOT = Path("/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/CORE")
OUT = HERE / "ELISA_diagnostics"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(HERE.parents[1] / ".matplotlib-cache"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, apply_axis_style  # noqa: E402

plt.rcParams.update({
    **ANALYSIS_RCPARAMS, "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": 12, "axes.titlesize": 14, "axes.labelsize": 13,
    "figure.dpi": 120, "savefig.dpi": 300,
})

COLORS = ["#79C6A3", "#38A6A5", "#287FB8", "#536FA8", "#8E6FA8", "#D06C75"]
DILUTIONS = [("1:3000", "D", "H"), ("1:300", "C", "G"), ("1:30", "B", "F"), ("1:3", "A", "E")]


def style(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    apply_axis_style(ax)


def save(fig, filename):
    fig.savefig(OUT / filename, bbox_inches="tight")
    plt.close(fig)


def read_plate(path, sheet):
    data = pd.read_excel(path, sheet_name=sheet, header=None)
    result = {}
    for _, row in data.iterrows():
        well = str(row.iloc[0]).strip()
        if len(well) >= 2 and well[0] in "ABCDEFGH" and well[1:].isdigit():
            result[well] = pd.to_numeric(row.iloc[1], errors="coerce")
    return result


def plot_bradford_workbook_elisa(path):
    # The second sheet is an ELISA and already includes sample labels,
    # protein concentrations and A450 values in columns B-D.
    data = pd.read_excel(path, sheet_name="Result sheet (2)", header=None).iloc[41:134, :4].copy()
    data.columns = ["well", "sample", "protein_concentration", "a450"]
    data["sample"] = data["sample"].astype(str).str.strip()
    data["protein_concentration"] = pd.to_numeric(data["protein_concentration"], errors="coerce")
    data["a450"] = pd.to_numeric(data["a450"], errors="coerce")
    data = data.dropna(subset=["protein_concentration", "a450"])

    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    for color, (sample, group) in zip(COLORS, data.groupby("sample", sort=False)):
        summary = group.groupby("protein_concentration")["a450"].agg(["mean", "std"]).sort_index()
        x = np.arange(len(summary))
        ax.errorbar(x, summary["mean"], yerr=summary["std"], marker="o", capsize=3,
                    linewidth=2, color=color, label=sample)
    concentrations = sorted(data["protein_concentration"].unique())
    ax.set_xticks(range(len(concentrations)), [f"{x:g}" for x in concentrations])
    ax.set_xlabel("Protein concentration (µg/mL)")
    ax.set_ylabel(r"Absorbance ($A_{450}$)")
    ax.set_title("29 March 2026 ELISA — mean ± SD")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    style(ax)
    fig.tight_layout()
    save(fig, "20260329_CORE_EV_COVA_ELISA.png")


def plot_gfpx2(path):
    groups = {"WT": [2, 3, 4], "CORE": [5, 6, 7], "PRp27": [8, 9, 10]}
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), sharey=True)
    for ax, sheet in zip(axes, ["Result sheet (1)", "Result sheet (2)"]):
        plate = read_plate(path, sheet)
        for color, (group, columns) in zip(COLORS, groups.items()):
            corrected = []
            for _, ev_row, sample_row in DILUTIONS:
                ev = np.array([plate[f"{ev_row}{col}"] for col in columns], float)
                sample = np.array([plate[f"{sample_row}{col}"] for col in columns], float)
                corrected.append(sample - ev)
            means = np.array([values.mean() for values in corrected])
            sds = np.array([values.std(ddof=1) for values in corrected])
            x = np.arange(4)
            ax.errorbar(x, means, yerr=sds, marker="o", capsize=3, linewidth=2,
                        color=color, label=group)
        ax.set_xticks(range(4), [item[0] for item in DILUTIONS])
        ax.set_xlabel("Leaf extract dilution")
        ax.set_title(sheet.replace("Result sheet", "Read"))
        style(ax)
    axes[0].set_ylabel(r"EV-corrected absorbance ($A_{450}$)")
    axes[1].legend(frameon=False)
    fig.suptitle("6 May 2026 GFP ELISA — mean ± SD", y=1.02)
    fig.tight_layout()
    save(fig, "20260506_WT_CORE_PRp27_GFP_ELISA.png")


def plot_cova(path):
    plate = read_plate(path, "Result sheet (3)")
    genotypes = {
        "Cas9 WT": [1, 2], "CORE 1-2": [3, 4], "CORE 2-8": [5, 6],
        "CORE 704-6": [7, 8], "CORE 705-1": [9, 10], "CORE 714-1": [11, 12],
    }
    rows = [("1:3000", "G", "H"), ("1:300", "E", "F"), ("1:30", "C", "D"), ("1:3", "A", "B")]
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    for color, (genotype, columns) in zip(COLORS, genotypes.items()):
        replicate_means = []
        for _, row1, row2 in rows:
            values = []
            for col in columns:
                wells = [f"{row1}{col}", f"{row2}{col}"]
                vals = [plate[w] for w in wells if not (w in {"D1", "D10"})]
                values.append(np.nanmean(vals))
            replicate_means.append(np.array(values))
        means = np.array([x.mean() for x in replicate_means])
        sds = np.array([x.std(ddof=1) for x in replicate_means])
        x = np.arange(4)
        ax.errorbar(x, means, yerr=sds, marker="o", capsize=3, linewidth=2,
                    color=color, label=genotype)
    ax.set_xticks(range(4), [x[0] for x in rows])
    ax.set_xlabel("Leaf extract dilution")
    ax.set_ylabel(r"Absorbance ($A_{450}$)")
    ax.set_title("12 March 2025 COVA ELISA — mean ± SD")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    style(ax)
    fig.tight_layout()
    save(fig, "20250312_CORE_COVA_ELISA.png")


def main():
    plot_bradford_workbook_elisa(
        ROOT / "GFP-COVA-VirB/202603_batch1_Full/ELISA/Bradford_andE_ELISA_1.xlsx"
    )
    plot_gfpx2(ROOT / "GFP-COVA-VirB/ELISAS/20260506_eve_elisa/gfpx2-elisa.xlsx")
    # The two COVA workbook paths are byte-identical, so plot the experiment once.
    plot_cova(ROOT / "GFP-COVA-VirB/CombinedTest/COVA_ELISA_6wk_CORECombined.xlsx")
    print(f"Saved 3 labelled ELISA plots to {OUT}")


if __name__ == "__main__":
    main()
