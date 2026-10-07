"""Plot the 20260612 PRp27 catalytic-mutant GFP experiment."""

from pathlib import Path
import os

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from patsy import build_design_matrices
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests


HERE = Path(__file__).resolve().parent
SOURCE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/"
    "prp27/newCatMuts/20260612/20260612_newCatMut_gfp.xlsx"
)
FIG_DIR = HERE / "figures"
LOG_DIR = HERE / "logs"
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

TREATMENTS = ["EV", "PRp27 OE", "KM001", "EW011"]
TREATMENT_LABELS = [
    "EV",
    "PRp27 OE",
    r"PRp27 OE$^{\mathsf{E123Q}}$",
    r"PRp27 OE$^{\mathsf{4mut}}$",
]
PLANT_LINES = ["WT", "prp27-1", "prp27-2"]
FINAL_LINES = ["prp27-1", "prp27-2"]
COLORS = ["#4E79A7", "#F28E2B", "#59A14F", "#B07AA1"]
MARKERS = ["o", "^", "s", "D"]
X_POSITIONS = np.arange(len(TREATMENTS)) * 1.48
TEXT_COLOR = "#4A4A4A"


def load_data():
    frame = pd.read_excel(SOURCE, sheet_name="Sheet1", header=1)
    frame = frame.dropna(how="all").dropna(axis=1, how="all")
    frame = frame.rename(
        columns={
            "Repeat": "repeat",
            "Plant Line": "plant_line",
            "Pre-infiltration": "treatment",
            "Background": "background",
            "GFP": "gfp",
            "Normalised GFP": "normalised_gfp",
        }
    )
    frame = frame[list({"repeat", "plant_line", "treatment", "background", "gfp", "normalised_gfp"})]
    frame["repeat"] = pd.to_numeric(frame["repeat"]).astype(int)
    frame["normalised_gfp"] = pd.to_numeric(frame["normalised_gfp"])
    return frame


def representative_repeat(frame):
    """Choose the repeat closest to the across-repeat median treatment profile."""
    means = (
        frame[frame["plant_line"].isin(FINAL_LINES)]
        .groupby(["repeat", "plant_line", "treatment"])["normalised_gfp"]
        .mean()
        .unstack("treatment")
    )
    ratios = means[TREATMENTS[1:]].div(means["EV"], axis=0)
    log_ratios = np.log(ratios)
    target = log_ratios.groupby("plant_line").median()
    rows = []
    for repeat in sorted(frame["repeat"].unique()):
        repeat_profile = log_ratios.xs(repeat, level="repeat")
        distance = (repeat_profile - target).abs().to_numpy().mean()
        rows.append({"repeat": repeat, "mean_absolute_log_ratio_distance": distance})
    scores = pd.DataFrame(rows).sort_values(["mean_absolute_log_ratio_distance", "repeat"])
    return int(scores.iloc[0]["repeat"]), scores


def anova_comparisons_for_repeat(frame, repeat):
    """Fit one-way ANOVAs and estimate the five planned treatment contrasts."""
    rows = []
    anova_rows = []
    selected = frame[(frame["repeat"] == repeat) & frame["plant_line"].isin(FINAL_LINES)]
    comparison_specs = [
        ("PRp27 OE", "KM001"),
        ("PRp27 OE", "EW011"),
        ("EV", "PRp27 OE"),
        ("EV", "KM001"),
        ("EV", "EW011"),
    ]
    for plant_line in FINAL_LINES:
        line_data = selected[selected["plant_line"] == plant_line].copy()
        line_data["treatment"] = pd.Categorical(
            line_data["treatment"], categories=TREATMENTS, ordered=False
        )
        model = ols("normalised_gfp ~ C(treatment)", data=line_data).fit()
        anova_table = anova_lm(model, typ=2).reset_index(names="term")
        for _, anova_row in anova_table.iterrows():
            anova_rows.append(
                {
                    "repeat": repeat,
                    "plant_line": plant_line,
                    "term": anova_row["term"],
                    "sum_sq": anova_row["sum_sq"],
                    "df": anova_row["df"],
                    "F": anova_row.get("F", np.nan),
                    "p": anova_row.get("PR(>F)", np.nan),
                }
            )

        design_info = model.model.data.design_info
        family = []
        for left, right in comparison_specs:
            contrast_grid = pd.DataFrame(
                {
                    "treatment": pd.Categorical(
                        [left, right], categories=TREATMENTS, ordered=False
                    )
                }
            )
            design = np.asarray(
                build_design_matrices([design_info], contrast_grid)[0]
            )
            result = model.t_test(design[0] - design[1])
            family.append(
                {
                    "repeat": repeat,
                    "plant_line": plant_line,
                    "left_group": left,
                    "right_group": right,
                    "comparison": f"{left} vs {right}",
                    "test": "planned contrast from one-way ANOVA",
                    "estimate": float(np.asarray(result.effect).ravel()[0]),
                    "t": float(np.asarray(result.tvalue).ravel()[0]),
                    "df": float(model.df_resid),
                    "raw_p": float(np.asarray(result.pvalue).ravel()[0]),
                }
            )
        adjusted = multipletests([row["raw_p"] for row in family], method="holm")[1]
        for row, adjusted_p in zip(family, adjusted):
            row["holm_adjusted_p"] = adjusted_p
            row["summary"] = "***" if adjusted_p < 0.001 else "**" if adjusted_p < 0.01 else "*" if adjusted_p < 0.05 else "ns"
        rows.extend(family)
    return pd.DataFrame(rows), pd.DataFrame(anova_rows)


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=1.2, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=TEXT_COLOR)
    ax.set_axisbelow(True)


def draw_bars(ax, subset, title=None, show_ylabel=True):
    tops = []
    for index, treatment in enumerate(TREATMENTS):
        values = subset.loc[subset["treatment"] == treatment, "normalised_gfp"].to_numpy() / 1000
        mean = values.mean()
        sem = stats.sem(values) if len(values) > 1 else 0
        tops.append(max(values.max(), mean + sem))
        x_position = X_POSITIONS[index]
        ax.bar(x_position, mean, width=0.58, color=COLORS[index], alpha=0.60, edgecolor="#303030", linewidth=1.2, zorder=2)
        ax.errorbar(x_position, mean, yerr=sem, color="#303030", capsize=4, linewidth=1.2, zorder=3)
        jitter = np.linspace(-0.10, 0.10, len(values))
        ax.scatter(
            x_position + jitter,
            values,
            s=42,
            marker=MARKERS[index],
            color=COLORS[index],
            alpha=0.90,
            edgecolor="#000000",
            linewidth=0.45,
            zorder=4,
        )
    ax.set_xticks(X_POSITIONS)
    ax.set_xticklabels(TREATMENT_LABELS, rotation=25, ha="right")
    if show_ylabel:
        ax.set_ylabel("GFP signal intensity (x1000)")
    if title:
        ax.set_title(title, fontstyle="italic" if title.startswith("prp27") else "normal")
    style_axis(ax)
    return max(tops)


def add_comparison_line(ax, left, right, y, label, step):
    ax.plot([left, right], [y, y], color="#303030", linewidth=1.1, clip_on=False)
    ax.text(
        (left + right) / 2,
        y + step * 0.10,
        label,
        ha="center",
        va="bottom",
        fontsize=23 if "*" in label else 20,
        color="#303030",
    )


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Helvetica Neue",
        "mathtext.sf": "Helvetica Neue",
        "mathtext.it": "Helvetica Neue:italic",
        "mathtext.bf": "Helvetica Neue:bold",
        "font.size": 20,
        "axes.labelsize": 20,
        "axes.titlesize": 20,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "axes.edgecolor": TEXT_COLOR,
        "axes.labelcolor": TEXT_COLOR,
        "text.color": TEXT_COLOR,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

data = load_data()
data.to_csv(LOG_DIR / "new_cat_mut_all_values.csv", index=False)
_, repeat_scores = representative_repeat(data)
chosen_repeat = 1
repeat_scores.to_csv(LOG_DIR / "representative_repeat_scores.csv", index=False)

# Diagnostic: every plant line in every repeat.
fig, axes = plt.subplots(3, 3, figsize=(16.2, 13.6), sharey=True)
for row, repeat in enumerate(sorted(data["repeat"].unique())):
    for column, plant_line in enumerate(PLANT_LINES):
        subset = data[(data["repeat"] == repeat) & (data["plant_line"] == plant_line)]
        draw_bars(axes[row, column], subset, title=plant_line, show_ylabel=column == 0)
        if column == 0:
            axes[row, column].text(-0.28, 0.5, f"Repeat {repeat}", transform=axes[row, column].transAxes, rotation=90, ha="center", va="center", fontweight="bold")
fig.subplots_adjust(left=0.11, right=0.985, top=0.95, bottom=0.09, wspace=0.18, hspace=0.48)
for extension in ("png", "svg", "pdf"):
    fig.savefig(FIG_DIR / f"202609_new_cat_mut_all_lines_all_repeats_diagnostic.{extension}", dpi=300, bbox_inches="tight")
plt.close(fig)

# Final: the objectively selected representative repeat, mutant backgrounds only.
comparisons, anova_results = anova_comparisons_for_repeat(data, chosen_repeat)
comparisons.to_csv(LOG_DIR / "representative_repeat_anova_holm_comparisons.csv", index=False)
anova_results.to_csv(LOG_DIR / "representative_repeat_one_way_anova.csv", index=False)
fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.6), sharey=True)
panel_tops = []
for index, (ax, plant_line) in enumerate(zip(axes, FINAL_LINES)):
    subset = data[(data["repeat"] == chosen_repeat) & (data["plant_line"] == plant_line)]
    top = draw_bars(ax, subset, title=plant_line, show_ylabel=index == 0)
    panel_tops.append(top)

global_top = max(panel_tops)
step = max(global_top * 0.16, 3.4)
for ax, plant_line, top in zip(axes, FINAL_LINES, panel_tops):
    line_comparisons = comparisons[comparisons["plant_line"] == plant_line].reset_index(drop=True)
    subset = data[(data["repeat"] == chosen_repeat) & (data["plant_line"] == plant_line)]
    treatment_top = subset.loc[subset["treatment"] != "EV", "normalised_gfp"].max() / 1000
    for comp_index, comparison in line_comparisons.iterrows():
        left = X_POSITIONS[TREATMENTS.index(comparison["left_group"])]
        right = X_POSITIONS[TREATMENTS.index(comparison["right_group"])]
        if comparison["left_group"] == "PRp27 OE":
            middle_base = treatment_top + 0.27 * (top - treatment_top)
            y = middle_base + step * comp_index
        else:
            y = top + step * (0.55 + comp_index - 2)
        add_comparison_line(ax, left, right, y, comparison["summary"], step)
for ax in axes:
    ax.set_ylim(-0.5, global_top + step * 4.1)
fig.subplots_adjust(left=0.105, right=0.985, top=0.92, bottom=0.27, wspace=0.18)
for extension in ("png", "svg", "pdf"):
    fig.savefig(FIG_DIR / f"202609_new_cat_mut_representative_repeat_{chosen_repeat}.{extension}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Selected representative repeat: {chosen_repeat}")
print(repeat_scores.to_string(index=False))
print(f"Saved figures to: {FIG_DIR}")
