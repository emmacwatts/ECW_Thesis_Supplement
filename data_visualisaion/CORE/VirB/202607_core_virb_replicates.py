"""CORE virB expression across three independent experiments."""

from pathlib import Path
import os
import sys

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

OUT = HERE / "202607_core_virb_replicates_figures"
LOG = HERE / "202607_core_virb_replicates_logs"
for folder in (OUT, LOG):
    folder.mkdir(parents=True, exist_ok=True)

SOURCE_1_2 = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/20260330/20260330_virB/20260330_virB.prism"
)
SOURCE_3 = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/20260318core/20260318_core1-2_gfp_and_virB.prism"
)

# Values transcribed from the named virB Prism data tables only.
VIRB_VALUES = {
    "Replicate 1": {
        "Cas9 WT": [398.869, 693.737, 524.171, 643.001],
        "core-1-2": [817.422, 853.039, 867.536, 1225.208],
        "source": SOURCE_1_2, "source_table": "virB_batch1",
    },
    "Replicate 2": {
        "Cas9 WT": [994.223, 945.077, 820.156, 681.895, 1087.818],
        "core-1-2": [900.491, 895.549, 999.027, 873.758, 1214.636],
        "source": SOURCE_1_2, "source_table": "virB_batch2",
    },
    "Replicate 3": {
        "Cas9 WT": [388.822, 420.298, 334.089],
        "core-1-2": [732.132, 959.671, 963.894],
        "source": SOURCE_3, "source_table": "virB",
    },
}
EXPERIMENT_ORDER = ["Replicate 3", "Replicate 1", "Replicate 2"]
GROUP_ORDER = ["Cas9 WT", "core-1-2"]
COLORS = {"Cas9 WT": "#A6B0BA", "core-1-2": "#31539A"}
BAR_OUTLINE = "#6F7780"

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": 17, "axes.titlesize": 18, "axes.labelsize": 18,
    "xtick.labelsize": 16, "ytick.labelsize": 16,
    "figure.dpi": 120, "savefig.dpi": 300,
})


def build_data():
    rows = []
    for experiment, details in VIRB_VALUES.items():
        for genotype in GROUP_ORDER:
            for replicate, value in enumerate(details[genotype], start=1):
                rows.append({
                    "experiment": experiment, "genotype": genotype,
                    "biological_replicate": replicate, "virb": value,
                    "source_file": str(details["source"]),
                    "source_table": details["source_table"],
                })
    return pd.DataFrame(rows)


def analyse(data):
    summary = (
        data.groupby(["experiment", "genotype"])["virb"]
        .agg(n="size", mean="mean", sd="std", sem="sem").reset_index()
    )
    tests = []
    experiment_effects = []
    for experiment in EXPERIMENT_ORDER:
        subset = data[data["experiment"].eq(experiment)]
        wt = subset.loc[subset["genotype"].eq("Cas9 WT"), "virb"].to_numpy()
        core = subset.loc[subset["genotype"].eq("core-1-2"), "virb"].to_numpy()
        result = stats.ttest_ind(wt, core, equal_var=False, alternative="two-sided")
        tests.append({
            "experiment": experiment, "comparison": "Cas9 WT vs core-1-2",
            "test": "two-sided Welch t-test", "n_wt": len(wt), "n_core": len(core),
            "mean_wt": wt.mean(), "mean_core": core.mean(),
            "core_over_wt_fold_change": core.mean() / wt.mean(),
            "t": result.statistic, "df": result.df, "p_raw": result.pvalue,
        })
        experiment_effects.append({
            "experiment": experiment,
            "mean_wt": wt.mean(), "mean_core": core.mean(),
            "core_over_wt_fold_change": core.mean() / wt.mean(),
            "log2_fold_change": np.log2(core.mean() / wt.mean()),
        })
    tests = pd.DataFrame(tests)
    tests["p_holm_three_experiments"] = multipletests(tests["p_raw"], method="holm")[1]
    effects = pd.DataFrame(experiment_effects)
    one_sample = stats.ttest_1samp(
        effects["log2_fold_change"], popmean=0, alternative="two-sided"
    )
    combined = pd.DataFrame([{
        "analysis_unit": "experiment mean fold change",
        "n_experiments": len(effects),
        "geometric_mean_core_over_wt": 2 ** effects["log2_fold_change"].mean(),
        "mean_log2_fold_change": effects["log2_fold_change"].mean(),
        "t": one_sample.statistic, "df": one_sample.df, "p": one_sample.pvalue,
        "test": "two-sided one-sample t-test of experiment log2 fold changes vs 0",
    }])
    return summary, tests, effects, combined


def p_label(p):
    if p >= 0.05:
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    return "*"


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    apply_axis_style(ax)


def draw_experiment(
    ax, data, summary, tests, experiment, display_number,
    common_ylim, show_ylabel=False,
):
    subset = data[data["experiment"].eq(experiment)]
    stats_df = summary[summary["experiment"].eq(experiment)].set_index("genotype").loc[GROUP_ORDER]
    x = np.arange(2)
    means = stats_df["mean"].to_numpy()
    sems = stats_df["sem"].to_numpy()
    ax.bar(
        x, means, width=0.58, yerr=sems, capsize=4,
        color=[COLORS[group] for group in GROUP_ORDER], alpha=0.84,
        edgecolor=BAR_OUTLINE, linewidth=0.8,
        error_kw={"elinewidth": 1.3, "ecolor": TEXT_COLOR, "capthick": 1.3},
    )
    rng = np.random.default_rng(202607 + EXPERIMENT_ORDER.index(experiment))
    for index, genotype in enumerate(GROUP_ORDER):
        values = subset.loc[subset["genotype"].eq(genotype), "virb"].to_numpy()
        ax.scatter(
            np.full(len(values), index) + rng.uniform(-0.08, 0.08, len(values)),
            values, s=32, color=COLORS[genotype], edgecolor="none", alpha=0.60, zorder=4,
        )
        ax.text(
            index, values.max() + common_ylim * 0.055,
            f"{len(values)}", ha="center", va="bottom",
            fontsize=12, color="#929BA5",
        )
    top = subset["virb"].max()
    y = top + common_ylim * 0.14
    h = common_ylim * 0.025
    ax.plot([0, 0, 1, 1], [y, y + h, y + h, y], color=TEXT_COLOR, lw=1.2)
    p = tests.loc[tests["experiment"].eq(experiment), "p_holm_three_experiments"].iloc[0]
    ax.text(0.5, y + h + top * 0.025, p_label(p), ha="center", va="bottom")
    ax.set_ylim(0, common_ylim)
    ax.set_xticks(x, ["Cas9 WT", "core-1-2"])
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_title(str(display_number))
    if show_ylabel:
        ax.set_ylabel("virB signal (a.u.)")
    style_axis(ax)


def draw_summary(ax, effects, combined):
    x = np.arange(len(effects))
    values = effects.set_index("experiment").loc[EXPERIMENT_ORDER, "core_over_wt_fold_change"]
    ax.axhline(1, color=TEXT_COLOR, linestyle=":", linewidth=1.1)
    ax.scatter(x, values, s=65, color="#31539A", edgecolor="white", linewidth=0.7, zorder=4)
    ax.plot(x, values, color="#8794A1", linewidth=1.2, zorder=2)
    repeat_labels = ["Repeat 1", "Repeat 2", "Repeat 3"]
    ax.set_xticks(x, repeat_labels)
    ax.set_ylabel(r"$\it{core}$-1-2 / Cas9 WT")
    ax.set_title("Across repeats")
    ax.set_ylim(0, max(2.8, values.max() * 1.25))
    style_axis(ax)


def make_figure(data, summary, tests, effects, combined):
    fig = plt.figure(figsize=(12.0, 8.2), facecolor="white")
    gs = fig.add_gridspec(
        2, 3, height_ratios=[1.0, 0.72], left=0.085, right=0.985,
        bottom=0.10, top=0.92, wspace=0.18, hspace=0.36,
    )
    axes = [fig.add_subplot(gs[0, index]) for index in range(3)]
    common_ylim = max(
        data["virb"].max(),
        max(
            summary["mean"] + summary["sem"]
        ),
    ) * 1.34
    for index, (ax, experiment) in enumerate(zip(axes, EXPERIMENT_ORDER)):
        draw_experiment(
            ax, data, summary, tests, experiment, index + 1, common_ylim,
            show_ylabel=(index == 0),
        )
    summary_ax = fig.add_subplot(gs[1, :])
    draw_summary(summary_ax, effects, combined)
    for ax, letter in [(axes[0], "A"), (summary_ax, "B")]:
        box = ax.get_position()
        fig.text(
            box.x0 - 0.035, box.y1 + 0.035, letter, fontsize=21,
            fontweight="bold", fontfamily="Helvetica Neue", color=TEXT_COLOR,
        )
    stem = "202607_core_virb_three_replicates"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    data = build_data()
    summary, tests, effects, combined = analyse(data)
    with pd.ExcelWriter(LOG / "core_virb_three_replicates_data_and_statistics.xlsx") as writer:
        data.to_excel(writer, sheet_name="virB_values", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        tests.to_excel(writer, sheet_name="within_experiment_tests", index=False)
        effects.to_excel(writer, sheet_name="experiment_fold_changes", index=False)
        combined.to_excel(writer, sheet_name="cross_experiment_test", index=False)
    (LOG / "core_virb_three_replicates_methods_note.txt").write_text(
        "virB values were extracted from the virB_batch1 and virB_batch2 tables in "
        f"{SOURCE_1_2}, and the virB table in {SOURCE_3}. PRp27 and all GFP tables were "
        "excluded. Bars show mean ± SEM with all biological measurements overlaid. "
        "Within each experiment, Cas9 WT and core-1-2 were compared using two-sided "
        "Welch t-tests with Holm adjustment across the three experiments. Cross-experiment "
        "consistency was summarised as the core-1-2/Cas9 WT mean fold change per experiment; "
        "a two-sided one-sample t-test of experiment-level log2 fold changes against zero "
        "uses experiments as the independent unit. With only three experiments, the "
        "cross-experiment test is low-powered and should be interpreted alongside the "
        "individual experiment panels.\n",
        encoding="utf-8",
    )
    make_figure(data, summary, tests, effects, combined)
    print(tests.to_string(index=False))
    print(combined.to_string(index=False))


if __name__ == "__main__":
    main()
