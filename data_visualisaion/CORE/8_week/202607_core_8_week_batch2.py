"""Eight-week CORE GFP endpoint and matched-age variability analysis."""

from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location
import os
import sys

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

SOURCE_8WK = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/8wk_CORE/20260212_8wkCORE_gfp/8wk_CORE_batch2_gfp.xlsx"
)
SOURCE_EARLIER = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/GFP-COVA-VirB/CombinedTest/GFP_20250227/Combined_Batch_GFP.xlsx"
)
OUT = HERE / "202607_core_8_week_batch2_figures"
LOG = HERE / "202607_core_8_week_batch2_logs"
for folder in (OUT, LOG):
    folder.mkdir(parents=True, exist_ok=True)

ORDER = ["Cas9 WT", "core-2-8", "core-1-2", "core-704-6"]
SOURCE_TO_LABEL = {
    "cas9 WT": "Cas9 WT", "core 2-8": "core-2-8",
    "core-1-2": "core-1-2", "core 704-6": "core-704-6",
}
EARLIER_TO_LABEL = {
    "WT": "Cas9 WT", "2-8": "core-2-8",
    "2-1": "core-1-2", "704-6": "core-704-6",
}
COLORS = {
    "Cas9 WT": "#A6B0BA", "core-2-8": "#536FA8",
    "core-1-2": "#31539A", "core-704-6": "#287FB8",
}
BAR_OUTLINE = "#6F7780"
Y_SCALE = 1000
BOOTSTRAP_ITERATIONS = 20000
ROS_SCRIPT = HERE.parent / "ROS_timecourse" / "202607_core_ros_age_timecourse.py"

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": 17, "axes.titlesize": 19, "axes.labelsize": 18,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15,
    "figure.dpi": 120, "savefig.dpi": 300,
})


def load_8wk():
    data = pd.read_excel(SOURCE_8WK, sheet_name="Batch2", header=1)
    data = data.rename(columns={
        "Line": "source_genotype", "p19?": "p19",
        "Background": "background", "Signal": "gfp_signal",
    })
    data = data[data["source_genotype"].isin(SOURCE_TO_LABEL)].copy()
    data["genotype"] = data["source_genotype"].map(SOURCE_TO_LABEL)
    data["normalised_gfp"] = data["gfp_signal"] - data["background"]
    data["experiment"] = "8 weeks"
    data["replicate"] = data.groupby(["p19", "genotype"]).cumcount() + 1
    return data[[
        "experiment", "p19", "genotype", "replicate",
        "background", "gfp_signal", "normalised_gfp",
    ]]


def load_earlier():
    data = pd.read_excel(SOURCE_EARLIER, sheet_name="Sheet1", header=4).iloc[:, 1:5]
    data.columns = ["source_genotype", "background", "gfp_signal", "normalised_gfp"]
    data = data[data["source_genotype"].isin(EARLIER_TO_LABEL)].copy()
    data["genotype"] = data["source_genotype"].map(EARLIER_TO_LABEL)
    data["experiment"] = "Earlier"
    data["p19"] = "p19"
    data["replicate"] = data.groupby("genotype").cumcount() + 1
    return data[[
        "experiment", "p19", "genotype", "replicate",
        "background", "gfp_signal", "normalised_gfp",
    ]]


def cv(values):
    values = np.asarray(values, dtype=float)
    return np.std(values, ddof=1) / np.mean(values) if len(values) > 1 and np.mean(values) else np.nan


def bootstrap_cv_ratio(earlier, older, rng):
    ratios = np.empty(BOOTSTRAP_ITERATIONS)
    for index in range(BOOTSTRAP_ITERATIONS):
        a = rng.choice(earlier, len(earlier), replace=True)
        b = rng.choice(older, len(older), replace=True)
        ratios[index] = cv(b) / cv(a) if cv(a) > 0 else np.nan
    ratios = ratios[np.isfinite(ratios)]
    return np.quantile(ratios, [0.025, 0.5, 0.975])


def endpoint_tables(data):
    summary = (
        data.groupby(["p19", "genotype"])["normalised_gfp"]
        .agg(n="size", mean="mean", sd="std", sem="sem").reset_index()
    )
    tests = []
    for p19 in ["no p19", "p19"]:
        subset = data[data["p19"].eq(p19)]
        wt = subset.loc[subset["genotype"].eq("Cas9 WT"), "normalised_gfp"].to_numpy()
        groups = [
            subset.loc[subset["genotype"].eq(genotype), "normalised_gfp"].to_numpy()
            for genotype in ORDER[1:]
        ]
        result = stats.dunnett(*groups, control=wt, alternative="two-sided", random_state=20260716)
        for genotype, statistic, p_value in zip(ORDER[1:], result.statistic, result.pvalue):
            tests.append({
                "p19": p19, "comparison": f"{genotype} vs Cas9 WT",
                "test": "two-sided Dunnett test", "statistic": statistic,
                "p_adjusted": p_value,
            })
    return summary, pd.DataFrame(tests)


def variability_tables(earlier, older):
    matched = pd.concat([earlier, older[older["p19"].eq("p19")]], ignore_index=True)
    summary_rows, test_rows = [], []
    rng = np.random.default_rng(20260716)
    for genotype in ORDER:
        a = matched.loc[
            matched["experiment"].eq("Earlier") & matched["genotype"].eq(genotype),
            "normalised_gfp",
        ].to_numpy()
        b = matched.loc[
            matched["experiment"].eq("8 weeks") & matched["genotype"].eq(genotype),
            "normalised_gfp",
        ].to_numpy()
        low, median, high = bootstrap_cv_ratio(a, b, rng)
        # Test proportional rather than raw-scale dispersion because the two
        # experiments differ greatly in mean signal.
        fligner = stats.fligner(
            np.log(np.clip(a, 1, None)), np.log(np.clip(b, 1, None)), center="median"
        )
        summary_rows.extend([
            {"experiment": "Earlier", "genotype": genotype, "n": len(a),
             "mean": a.mean(), "sd": a.std(ddof=1), "cv": cv(a)},
            {"experiment": "8 weeks", "genotype": genotype, "n": len(b),
             "mean": b.mean(), "sd": b.std(ddof=1), "cv": cv(b)},
        ])
        test_rows.append({
            "genotype": genotype,
            "test": "Fligner-Killeen test of equal log-scale variance",
            "fligner_statistic": fligner.statistic, "p_raw": fligner.pvalue,
            "cv_ratio_8wk_over_earlier": cv(b) / cv(a),
            "bootstrap_cv_ratio_median": median,
            "bootstrap_cv_ratio_ci_low": low,
            "bootstrap_cv_ratio_ci_high": high,
        })
    tests = pd.DataFrame(test_rows)
    tests["p_holm"] = multipletests(tests["p_raw"], method="holm")[1]

    # Brown-Forsythe-style interaction on log-scale absolute median deviations.
    matched = matched.copy()
    matched["log_gfp"] = np.log(matched["normalised_gfp"].clip(lower=1))
    matched["group_median"] = matched.groupby(
        ["experiment", "genotype"]
    )["log_gfp"].transform("median")
    matched["absolute_log_deviation"] = (matched["log_gfp"] - matched["group_median"]).abs()
    matched["line_class"] = np.where(matched["genotype"].eq("Cas9 WT"), "WT", "CORE")
    model = smf.ols(
        "absolute_log_deviation ~ C(experiment) * C(line_class) + C(genotype)",
        data=matched,
    ).fit(cov_type="HC3")
    interaction_term = "C(experiment)[T.Earlier]:C(line_class)[T.WT]"
    pooled = pd.DataFrame([{
        "analysis": "log-scale Brown-Forsythe interaction",
        "question": "Does the experiment-related dispersion change differ between CORE and WT?",
        "term": interaction_term,
        "coefficient": model.params.get(interaction_term, np.nan),
        "robust_se_hc3": model.bse.get(interaction_term, np.nan),
        "t": model.tvalues.get(interaction_term, np.nan),
        "p": model.pvalues.get(interaction_term, np.nan),
        "n": len(matched),
    }])
    return pd.DataFrame(summary_rows), tests, matched, pooled


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    apply_axis_style(ax)


def load_ros_module():
    spec = spec_from_file_location("core_ros_age_timecourse", ROS_SCRIPT)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def italicise_core_legend(legend):
    if legend is None:
        return
    for text in legend.get_texts():
        if text.get_text().lower().startswith("core"):
            text.set_fontstyle("italic")


def p_to_label(p_value):
    if p_value >= 0.05:
        return "ns"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def draw_endpoint(ax, data, summary, tests):
    p19_order = ["no p19", "p19"]
    centers = np.arange(2, dtype=float)
    offsets = np.linspace(-0.30, 0.30, len(ORDER))
    width = 0.19
    rng = np.random.default_rng(8)
    for offset, genotype in zip(offsets, ORDER):
        stats_df = summary[summary["genotype"].eq(genotype)].set_index("p19").loc[p19_order]
        positions = centers + offset
        ax.bar(
            positions, stats_df["mean"] / Y_SCALE, width=width,
            yerr=stats_df["sem"] / Y_SCALE, capsize=3,
            color=COLORS[genotype], alpha=0.84, edgecolor=BAR_OUTLINE, linewidth=0.8,
            error_kw={"elinewidth": 1.2, "ecolor": TEXT_COLOR, "capthick": 1.2},
            label=genotype,
        )
        for position, condition in zip(positions, p19_order):
            values = data[
                data["genotype"].eq(genotype) & data["p19"].eq(condition)
            ]["normalised_gfp"].to_numpy() / Y_SCALE
            ax.scatter(
                np.full(len(values), position) + rng.uniform(-0.025, 0.025, len(values)),
                values, s=25, color=COLORS[genotype], edgecolor="none", alpha=0.55, zorder=4,
            )
    data_top = data["normalised_gfp"].max() / Y_SCALE
    bracket_height = 0.55
    bracket_bases = {"no p19": 20.0, "p19": 27.5}
    bracket_steps = {"no p19": 6.0, "p19": 6.0}
    for condition_index, condition in enumerate(p19_order):
        condition_tests = tests[tests["p19"].eq(condition)]
        wt_x = centers[condition_index] + offsets[0]
        for level, genotype in enumerate(ORDER[1:]):
            row = condition_tests[
                condition_tests["comparison"].eq(f"{genotype} vs Cas9 WT")
            ].iloc[0]
            core_x = centers[condition_index] + offsets[ORDER.index(genotype)]
            y = bracket_bases[condition] + bracket_steps[condition] * level
            ax.plot(
                [wt_x, wt_x, core_x, core_x],
                [y, y + bracket_height, y + bracket_height, y],
                color=TEXT_COLOR, linewidth=1.1, clip_on=False,
            )
            ax.text(
                (wt_x + core_x) / 2,
                y + bracket_height + bracket_steps[condition] * 0.08,
                p_to_label(row["p_adjusted"]), ha="center", va="bottom",
                fontsize=15, color=TEXT_COLOR,
            )
    ax.set_xticks(centers, ["No p19", "p19"])
    ax.set_ylabel("Normalised GFP signal (×10³ a.u.)")
    legend = ax.legend(
        frameon=False, ncol=4, mode="expand", loc="lower left",
        bbox_to_anchor=(0.01, 1.015, 0.98, 0.01), borderaxespad=0.0,
        columnspacing=1.0, handletextpad=0.55,
    )
    italicise_core_legend(legend)
    style_axis(ax)
    ax.set_ylim(0, 49)


def draw_variability(ax, variability):
    x = np.arange(len(ORDER))
    offsets = {"Earlier": -0.10, "8 weeks": 0.10}
    markers = {"Earlier": "o", "8 weeks": "s"}
    for experiment in ["Earlier", "8 weeks"]:
        subset = variability[variability["experiment"].eq(experiment)].set_index("genotype").loc[ORDER]
        ax.scatter(
            x + offsets[experiment], subset["cv"] * 100, s=62,
            marker=markers[experiment], color=[COLORS[g] for g in ORDER],
            edgecolor=TEXT_COLOR, linewidth=0.7, label=experiment, zorder=4,
        )
    for index, genotype in enumerate(ORDER):
        values = variability[variability["genotype"].eq(genotype)].set_index("experiment")
        ax.plot(
            [index - 0.10, index + 0.10],
            [values.loc["Earlier", "cv"] * 100, values.loc["8 weeks", "cv"] * 100],
            color="#9AA3AC", linewidth=1.2, zorder=2,
        )
    ax.axhline(0, color=TEXT_COLOR, linewidth=0.8)
    ax.set_xticks(x, ["Cas9\nWT", "core-\n2-8", "core-\n1-2", "core-\n704-6"])
    for label in ax.get_xticklabels()[1:]:
        label.set_fontstyle("italic")
    ax.set_ylabel("Coefficient of variation (%)")
    ax.set_title("p19 variability across experiments")
    ax.legend(frameon=False, loc="upper left")
    style_axis(ax)


def make_figure(data, summary, endpoint_tests):
    ros = load_ros_module()
    _, _, ros_plant_curves, _ = ros.process_experiments()
    ros_curve_summary = ros.curve_summary(ros_plant_curves)
    fig, axes = plt.subplots(2, 1, figsize=(8.27, 11.69))
    ros.draw_curves(axes[0], ros_curve_summary, "csp22", show_legend=True)
    draw_endpoint(axes[1], data, summary, endpoint_tests)
    for ax in axes:
        ax.xaxis.label.set_fontsize(20)
        ax.yaxis.label.set_fontsize(20)
        ax.title.set_fontsize(21)
        ax.tick_params(axis="both", labelsize=18)
        legend = ax.get_legend()
        if legend is not None:
            for text in legend.get_texts():
                text.set_fontsize(17)
            if legend.get_title() is not None:
                legend.get_title().set_fontsize(17)
    fig.subplots_adjust(
        left=0.14, right=0.965, top=0.955, bottom=0.08, hspace=0.39
    )
    for ax, letter in zip(axes, "AB"):
        box = ax.get_position()
        fig.text(
            box.x0 - 0.055, box.y1 + 0.020, letter, fontsize=23,
            fontweight="bold", fontfamily="Helvetica Neue", color=TEXT_COLOR,
        )
    stem = "202607_core_8_week_ros_timecourse_and_gfp"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    older = load_8wk()
    summary, endpoint_tests = endpoint_tables(older)
    with pd.ExcelWriter(LOG / "core_8_week_batch2_processed_data_and_statistics.xlsx") as writer:
        older.to_excel(writer, sheet_name="8wk_batch2_data", index=False)
        summary.to_excel(writer, sheet_name="8wk_endpoint_summary", index=False)
        endpoint_tests.to_excel(writer, sheet_name="8wk_dunnett_tests", index=False)
    make_figure(older, summary, endpoint_tests)
    (LOG / "core_8_week_batch2_methods_note.txt").write_text(
        "Eight-week endpoint data were read from Batch2 of " + str(SOURCE_8WK) +
        ". core-705-1 was excluded. Panel A shows the csp22 ROS timecourse across plant "
        "ages 28–40 days, using the processed age-series ROS data and excluding the "
        "failed 42-day experiment. Curves show plant means ± SEM. Panel B shows mean ± "
        "SEM and all observations for "
        "Cas9 WT, core-2-8, core-1-2 and core-704-6 under no-p19 and p19 conditions; "
        "each CORE line was compared with Cas9 WT using two-sided Dunnett tests within "
        "condition. Normalised GFP was calculated as signal minus background.\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(endpoint_tests.to_string(index=False))


if __name__ == "__main__":
    main()
