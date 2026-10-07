"""Replot the COVA intact-antibody ELISA dilution verification."""

from pathlib import Path
import os
import sys

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style  # noqa: E402

SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "3. COVID Antibodies/ELISA/ELISATecanReads/20231105_COVATE.xlsx"
)
SOURCE_SHEET = "Processed"
PRIOR_PLOT = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "3. COVID Antibodies/ELISA/ELISATecanReads/Plots/"
    "COVA TE Dilutions EVHCLCHC trunc.pdf"
)
SCRIPT_STEM = Path(__file__).stem
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["EV", "HC", "LCHC"]
SAMPLE_LABELS = {"EV": "EV", "HC": "HC", "LCHC": "LCHC"}
COLORS = {"EV": "#6D6D6D", "HC": "#38A6A5", "LCHC": "#31539A"}
MARKERS = {"EV": "o", "HC": "s", "LCHC": "^"}
DILUTION_ORDER = ["1:1000", "1:200", "1:100", "1:20", "1:10"]
DILUTION_VALUE = {label: 1 / float(label.split(":")[1]) for label in DILUTION_ORDER}

METHODS = (
    "COVA intact-antibody ELISA absorbance values were read from the Processed sheet "
    "of 20231105_COVATE.xlsx. The displayed dilution range (1:1000 to 1:10) matches the "
    "previous truncated plot; the 1:2 measurement was deliberately omitted from this "
    "view. Points show individual plate replicates and symbols/lines show mean +/- SD. "
    "EV and HC have n=4 per dilution and LCHC has n=2 per dilution. No background "
    "subtraction, outlier removal or transformation was applied. The LCHC linear fit is "
    "provided in the logs as a descriptive calibration check, not an inferential model. "
    f"Prior plot reference: {PRIOR_PLOT}."
)

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 17,
    "axes.titlesize": 19,
    "axes.labelsize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 15,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})


def read_data():
    data = pd.read_excel(SOURCE, sheet_name=SOURCE_SHEET, header=3, usecols="A:D")
    data.columns = ["sample", "replicate", "dilution", "a450"]
    data = data.dropna(subset=["sample", "replicate", "dilution", "a450"]).copy()
    data["sample"] = data["sample"].astype(str).str.strip()
    data["dilution"] = data["dilution"].astype(str).str.strip()
    data["a450"] = pd.to_numeric(data["a450"], errors="coerce")
    data = data[
        data["sample"].isin(SAMPLE_ORDER) & data["dilution"].isin(DILUTION_ORDER)
    ].copy()
    data["dilution_fraction"] = data["dilution"].map(DILUTION_VALUE)
    return data.sort_values(["sample", "dilution_fraction", "replicate"]).reset_index(drop=True)


def summarise(data):
    return (
        data.groupby(["sample", "dilution", "dilution_fraction"], sort=False)
        .agg(n=("a450", "count"), mean_a450=("a450", "mean"),
             sd_a450=("a450", "std"), sem_a450=("a450", "sem"))
        .reset_index()
        .sort_values(["sample", "dilution_fraction"])
    )


def lchc_fit(summary):
    data = summary[summary["sample"].eq("LCHC")].sort_values("dilution_fraction")
    fit = stats.linregress(data["dilution_fraction"], data["mean_a450"])
    return pd.DataFrame([{
        "sample": "LCHC", "model": "mean A450 ~ dilution fraction",
        "slope": fit.slope, "intercept": fit.intercept,
        "r_squared": fit.rvalue ** 2, "p_value_slope": fit.pvalue,
        "n_dilutions": len(data),
    }])


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(True, color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def plot(data, summary):
    fig, ax = plt.subplots(figsize=(8.8, 6.6))
    rng = np.random.default_rng(202607)
    for sample in SAMPLE_ORDER:
        raw = data[data["sample"].eq(sample)]
        stats_df = summary[summary["sample"].eq(sample)].sort_values("dilution_fraction")
        x = stats_df["dilution_fraction"].to_numpy()
        mean = stats_df["mean_a450"].to_numpy()
        sd = stats_df["sd_a450"].to_numpy()
        color = COLORS[sample]

        for x_value in x:
            values = raw.loc[np.isclose(raw["dilution_fraction"], x_value), "a450"].to_numpy()
            jitter = rng.uniform(-0.0014, 0.0014, len(values))
            ax.scatter(np.full(len(values), x_value) + jitter, values, s=30,
                       color=color, alpha=0.48, edgecolor="none", zorder=3)

        ax.errorbar(
            x, mean, yerr=sd, color=color, marker=MARKERS[sample], markersize=7,
            markeredgecolor="white", markeredgewidth=0.7, linewidth=2.0,
            elinewidth=1.25, capsize=4,
            label=SAMPLE_LABELS[sample], zorder=4,
        )

    ax.set_xlim(-0.004, 0.108)
    ax.set_ylim(0, 1.0)
    ax.set_xticks([0, 0.01, 0.05, 0.10])
    ax.set_xticklabels(["0", "0.01", "0.05", "0.10"])
    ax.set_xlabel("Leaf extract dilution")
    ax.set_ylabel(r"Absorbance ($A_{450}$)")
    style_axis(ax)
    ax.legend(frameon=False, title="Sample", title_fontsize=18,
              loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.subplots_adjust(left=0.15, right=0.78, top=0.96, bottom=0.15)
    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(SAVE_DIR / f"{SCRIPT_STEM}_dilution_response.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    data = read_data()
    summary = summarise(data)
    fit = lchc_fit(summary)
    data.to_excel(LOG_DIR / "cova_elisa_raw_values.xlsx", index=False)
    summary.to_excel(LOG_DIR / "cova_elisa_summary_mean_sd.xlsx", index=False)
    fit.to_excel(LOG_DIR / "cova_elisa_lchc_linear_fit.xlsx", index=False)
plot(data, summary)
    print(summary.to_string(index=False))
    print(fit.to_string(index=False))
    print(f"Figures saved to: {SAVE_DIR}")
    return data, summary, fit


if __name__ == "__main__":
    main()
