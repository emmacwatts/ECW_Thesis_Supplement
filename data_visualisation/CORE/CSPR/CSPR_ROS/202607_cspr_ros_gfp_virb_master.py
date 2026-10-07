"""CSPR ROS, GFP and virB master figure."""

from pathlib import Path
import importlib.util
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

ROS_SCRIPT = HERE / "202607_cspr_ros_timecourse.py"
spec = importlib.util.spec_from_file_location("cspr_ros", ROS_SCRIPT)
cspr_ros = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cspr_ros)

SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CSPR/CSPR_6wk/Batch1.prism"
)
OUT = HERE / "202607_cspr_ros_gfp_virb_master_figures"
LOG = HERE / "202607_cspr_ros_gfp_virb_master_logs"
for folder in (OUT, LOG):
    folder.mkdir(parents=True, exist_ok=True)

GENOTYPES = ["WT", "cspr-706-4", "cspr-3-1-1"]
COLORS = {"WT": "#A6B0BA", "cspr-706-4": "#2A9D8F", "cspr-3-1-1": "#31539A"}
OUTLINE = "#6F7780"
GFP = {
    "No p19": {
        "WT": [5399.600, 11398.353, 4584.751, 6137.991],
        "cspr-706-4": [4937.675, 11018.915, 4997.989, 5647.688],
        "cspr-3-1-1": [10357.261, 6013.016, 6111.253, 5861.935],
    },
    "p19": {
        "WT": [14975.851, 18285.164, 16706.006, 17414.055],
        "cspr-706-4": [13249.938, 11209.940, 13138.512, 16936.692],
        "cspr-3-1-1": [13885.596, 20382.128, 12957.490, 18034.022],
    },
}
VIRB = {
    "WT": [237.689, 227.238, 245.452],
    "cspr-706-4": [200.649, 223.538, 204.346, 251.523],
    "cspr-3-1-1": [222.490, 310.035, 236.461, 178.716],
}

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": 15, "axes.titlesize": 15.5, "axes.labelsize": 15,
    "xtick.labelsize": 13.5, "ytick.labelsize": 13.5,
    "legend.fontsize": 13.5, "figure.dpi": 120, "savefig.dpi": 300,
})


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    apply_axis_style(ax)


def italicise_cspr_ticks(ax):
    for label in ax.get_xticklabels():
        if "cspr" in label.get_text():
            label.set_fontstyle("italic")


def bracket(ax, x1, x2, y, height, text="ns"):
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y],
            color=TEXT_COLOR, linewidth=1.0, clip_on=False)
    ax.text((x1 + x2) / 2, y + height * 1.18, text,
            ha="center", va="bottom", fontsize=12.5, color=TEXT_COLOR)


def draw_ros(fig, gs, curves):
    sub = gs.subgridspec(1, 3, wspace=0.25)
    axes = [fig.add_subplot(sub[0, i]) for i in range(3)]
    conditions = cspr_ros.CONDITION_ORDER
    groups = cspr_ros.GROUP_ORDER
    colors = cspr_ros.COLORS
    maxima = {
        condition: (
            curves.loc[curves["condition"].eq(condition), "mean_luminescence"]
            + curves.loc[curves["condition"].eq(condition), "sem_luminescence"]
        ).max()
        for condition in conditions
    }
    shared = cspr_ros.round_up(max(maxima["water"], maxima["flg22"]) * 1.08)
    limits = {
        "water": shared, "flg22": shared,
        "csp22": cspr_ros.round_up(maxima["csp22"] * 1.08),
    }
    handles = []
    for ax, condition in zip(axes, conditions):
        for genotype in groups:
            subset = curves[
                curves["condition"].eq(condition) & curves["genotype"].eq(genotype)
            ].sort_values("time_h")
            x = subset["time_h"].to_numpy(float)
            y = subset["mean_luminescence"].to_numpy(float) / 1000
            error = subset["sem_luminescence"].to_numpy(float) / 1000
            ax.fill_between(x, y - error, y + error, color=colors[genotype],
                            alpha=0.17, linewidth=0)
            line, = ax.plot(x, y, color=colors[genotype], linewidth=1.8)
            if ax is axes[0]:
                handles.append(line)
        ax.set_title(condition, pad=5)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, limits[condition] / 1000)
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        style_axis(ax)
    axes[0].set_ylabel("Luminescence (AU ×1000)")
    fig.text(0.52, axes[0].get_position().y0 - 0.045,
             "Time after elicitation (h)", ha="center", va="top", fontsize=15)
    labels = ["WT", r"$\it{core}$-1-2", r"$\it{cspr}$-706-4", r"$\it{cspr}$-3-1-1"]
    axes[1].legend(
        handles, labels, frameon=False, ncol=4, loc="lower center",
        bbox_to_anchor=(0.5, 1.18), columnspacing=1.2, handlelength=1.8,
    )
    return axes


def draw_gfp(ax):
    condition_centres = np.array([0.0, 1.25])
    offsets = np.array([-0.28, 0.0, 0.28])
    width = 0.24
    rng = np.random.default_rng(202607)
    for condition_index, condition in enumerate(GFP):
        for genotype_index, genotype in enumerate(GENOTYPES):
            values = np.asarray(GFP[condition][genotype]) / 1000
            x = condition_centres[condition_index] + offsets[genotype_index]
            mean, sem = values.mean(), values.std(ddof=1) / np.sqrt(len(values))
            ax.bar(x, mean, width=width, yerr=sem, capsize=3,
                   color=COLORS[genotype], alpha=0.84,
                   edgecolor=OUTLINE, linewidth=0.8,
                   error_kw={"ecolor": TEXT_COLOR, "elinewidth": 1.1, "capthick": 1.1})
            ax.scatter(
                np.full(len(values), x) + rng.uniform(-0.045, 0.045, len(values)),
                values, s=22, color=COLORS[genotype], alpha=0.58,
                edgecolor="none", zorder=4,
            )
    ax.set_xticks(condition_centres, ["No p19", "p19"])
    ax.set_ylabel("Normalised GFP signal (×10³ a.u.)")
    ax.set_ylim(0, 28)
    # Prism Dunn/Dunnett comparisons: WT vs each CSPR line, all adjusted P > 0.05.
    for centre in condition_centres:
        bracket(ax, centre - 0.28, centre, 22.2, 0.45)
        bracket(ax, centre - 0.28, centre + 0.28, 24.3, 0.45)
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=COLORS[g], edgecolor=OUTLINE, linewidth=0.8)
        for g in GENOTYPES
    ]
    labels = ["WT", r"$\it{cspr}$-706-4", r"$\it{cspr}$-3-1-1"]
    style_axis(ax)
    return handles, labels


def draw_virb(ax):
    x = np.arange(3)
    rng = np.random.default_rng(202608)
    for index, genotype in enumerate(GENOTYPES):
        values = np.asarray(VIRB[genotype])
        mean, sem = values.mean(), values.std(ddof=1) / np.sqrt(len(values))
        ax.bar(index, mean, width=0.62, yerr=sem, capsize=3,
               color=COLORS[genotype], alpha=0.84,
               edgecolor=OUTLINE, linewidth=0.8,
               error_kw={"ecolor": TEXT_COLOR, "elinewidth": 1.1, "capthick": 1.1})
        ax.scatter(
            np.full(len(values), index) + rng.uniform(-0.07, 0.07, len(values)),
            values, s=25, color=COLORS[genotype], alpha=0.60,
            edgecolor="none", zorder=4,
        )
    ax.set_xticks(x, ["WT", "cspr-706-4", "cspr-3-1-1"])
    italicise_cspr_ticks(ax)
    ax.set_ylabel("virB signal (a.u.)")
    ax.set_ylim(0, 410)
    bracket(ax, 0, 1, 330, 7)
    bracket(ax, 0, 2, 365, 7)
    style_axis(ax)


def make_figure(curves):
    fig = plt.figure(figsize=(14.0, 8.8), facecolor="white")
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.02, 1.0], width_ratios=[1.55, 1.0],
        left=0.07, right=0.985, bottom=0.145, top=0.88,
        hspace=0.48, wspace=0.25,
    )
    ros_axes = draw_ros(fig, gs[0, :], curves)
    gfp_ax = fig.add_subplot(gs[1, 0])
    virb_ax = fig.add_subplot(gs[1, 1])
    lower_handles, lower_labels = draw_gfp(gfp_ax)
    draw_virb(virb_ax)
    lower_left = gfp_ax.get_position().x0
    lower_right = virb_ax.get_position().x1
    fig.legend(
        lower_handles, lower_labels, frameon=False, ncol=3,
        loc="lower center",
        bbox_to_anchor=((lower_left + lower_right) / 2, 0.035),
        handlelength=1.3, columnspacing=1.4,
    )
    for ax, letter in [(ros_axes[0], "A"), (gfp_ax, "B"), (virb_ax, "C")]:
        box = ax.get_position()
        fig.text(box.x0 - 0.035, box.y1 + 0.035, letter,
                 fontsize=21, fontweight="bold", fontfamily="Helvetica Neue",
                 color=TEXT_COLOR, ha="left", va="bottom")
    stem = "202607_cspr_ros_gfp_virb_master"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_logs(outliers):
    rows = []
    for assay, groups in [("GFP no p19", GFP["No p19"]), ("GFP p19", GFP["p19"]), ("virB", VIRB)]:
        for genotype, values in groups.items():
            for replicate, value in enumerate(values, 1):
                rows.append({
                    "assay": assay, "genotype": genotype,
                    "replicate": replicate, "value": value, "source_file": str(SOURCE),
                })
    with pd.ExcelWriter(LOG / "cspr_batch1_gfp_virb_data_and_statistics.xlsx") as writer:
        pd.DataFrame(rows).to_excel(writer, sheet_name="extracted_values", index=False)
        pd.DataFrame([
            {"assay": "GFP no p19", "omnibus_test": "Kruskal-Wallis",
             "omnibus_p": 0.6667, "posthoc": "Dunn vs WT",
             "WT_vs_cspr_706_4_adjusted_p": ">0.9999",
             "WT_vs_cspr_3_1_1_adjusted_p": ">0.9999"},
            {"assay": "GFP p19", "omnibus_test": "ordinary one-way ANOVA",
             "omnibus_p": 0.2223, "posthoc": "Dunnett vs WT",
             "WT_vs_cspr_706_4_adjusted_p": "0.1894",
             "WT_vs_cspr_3_1_1_adjusted_p": "0.9404"},
            {"assay": "virB", "omnibus_test": "ordinary one-way ANOVA",
             "omnibus_p": 0.7705, "posthoc": "Dunnett vs WT",
             "WT_vs_cspr_706_4_adjusted_p": "0.7716",
             "WT_vs_cspr_3_1_1_adjusted_p": ">0.9999"},
        ]).to_excel(writer, sheet_name="Prism_statistics", index=False)
        outliers.to_excel(writer, sheet_name="ROS_auc_outliers", index=False)
def main():
    _, outliers, curves = cspr_ros.process()
    make_figure(curves)
    write_logs(outliers)
    print(f"Figure: {OUT / '202607_cspr_ros_gfp_virb_master.pdf'}")


if __name__ == "__main__":
    main()
