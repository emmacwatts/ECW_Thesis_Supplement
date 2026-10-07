# %% [markdown]
# # RNA / actin-normalised expression timecourse
#
# Run this file one `# %%` block at a time with VS Code's **Run Cell** action.
# The script reads the original Excel workbook, uses the `RNA` sheet, and treats
# `Actin normalised` as the response variable.

# %%
from itertools import combinations
from pathlib import Path
import os

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from patsy import build_design_matrices
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests
from data_config import resolve_timecourse_workbook
from figure_styles import (
    ANALYSIS_RCPARAMS,
    COLORS,
    ERRORBAR_COLOR,
    FONT_STACK,
    INDIVIDUAL_POINT_ALPHA,
    LETTER_FONT_SIZE,
    LETTER_LABEL_GAP_FRACTION,
    LETTER_LABEL_STEP_FRACTION,
    MARKERS,
    MEAN_MARKER_EDGE,
    MEAN_MARKER_LINEWIDTH,
    apply_axis_style,
    make_errorbars_black,
)

# ---------------- CONFIG: edit to match your other figures ----------------
if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
elif (Path.cwd() / "PRp27" / "PRp27_timecourse_202606").is_dir():
    HERE = Path.cwd() / "PRp27" / "PRp27_timecourse_202606"
elif (Path.cwd() / "PRp27_timecourse_202606").is_dir():
    HERE = Path.cwd() / "PRp27_timecourse_202606"
else:
    HERE = Path.cwd()

SHEET = "RNA"


SRC = resolve_timecourse_workbook(HERE, env_var="RNA_SRC")
SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202606_RNA"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
STATS_DIR = LOG_DIR
SAVE_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)
SAVE_FMTS = ["png", "pdf", "svg"]

LINES = ["WT", "PRp27#1", "PRp27#2"]
SHORT_LINE_LABELS = {"WT": "WT", "PRp27#1": "#1", "PRp27#2": "#2"}
DPIS = [0, 1, 3, 5]
ERRBAR = "SEM"  # "SD" or "SEM"
SIG_ALPHA = 0.05
YLABEL = "Actin-normalised RNA signal"
XLABEL = "dpi"
FONT_FAMILY = FONT_STACK

# Outlier handling: flag for manual review, do not remove automatically.
# The RNA sheet contains biological samples, but no repeat/batch column.
# Flags use a two-sided Grubbs test within line × dpi groups.
OUTLIER_ALPHA = 0.05
OUTLIER_GROUP_COLS = ["line", "dpi"]
MANUAL_EXCLUSIONS = [
    # Example after manual confirmation:
    # {"line": "PRp27#2", "dpi": 5, "bio_rep": 2,
    #  "reason": "documented technical failure"},
]

# --------------------------------------------------------------------------

plt.rcParams.update(ANALYSIS_RCPARAMS)

print(f"Setup complete. Using data file: {SRC}")
print(f"Figures will save to: {SAVE_DIR}")


# %%
def load_raw(src):
    """Load RNA sheet and standardise column names/types."""
    # The RNA sheet has a title row above the true headers.
    df = pd.read_excel(src, sheet_name=SHEET, header=2)
    df = df.rename(
        columns={
            "Sample Type": "line",
            "dpi": "dpi",
            "Background": "background",
            "GFP Signal": "gfp_signal",
            "Background Normalised": "background_norm",
            "Actin control": "actin_control",
            "Actin normalised": "actin_norm",
        }
    )
    keep = [
        "line",
        "dpi",
        "background",
        "gfp_signal",
        "background_norm",
        "actin_control",
        "actin_norm",
    ]
    df = df[keep].copy()
    df["line"] = pd.Categorical(df["line"], categories=LINES, ordered=False)
    df["dpi"] = df["dpi"].astype(int)
    df["actin_norm"] = pd.to_numeric(df["actin_norm"], errors="coerce")
    df["bio_rep"] = (
        df.groupby(["line", "dpi"], observed=True)
        .cumcount()
        .add(1)
    )
    return df.dropna(subset=["line", "dpi", "actin_norm"])


def grubbs_two_sided_flag(values):
    """Return the most extreme value and two-sided Grubbs-test p-value."""
    values = values.dropna()
    n = len(values)
    if n < 3:
        return None

    sd = values.std(ddof=1)
    if sd == 0 or pd.isna(sd):
        return None

    mean = values.mean()
    deviations = (values - mean).abs()
    candidate_index = deviations.idxmax()
    g_statistic = deviations.loc[candidate_index] / sd

    denominator = (n - 1) ** 2 - n * g_statistic**2
    if denominator <= 0:
        p_value = 0.0
    else:
        t_squared = g_statistic**2 * n * (n - 2) / denominator
        p_value = 2 * n * (1 - stats.t.cdf(np.sqrt(t_squared), df=n - 2))
        p_value = max(0.0, min(1.0, float(p_value)))

    return {
        "index": candidate_index,
        "n": n,
        "group_mean": mean,
        "group_sd": sd,
        "g_statistic": float(g_statistic),
        "p_value": p_value,
    }


def clean_data(raw):
    """Flag potential outliers and remove only manually confirmed exclusions."""
    screened = raw.copy()
    screened["flag_outlier"] = False
    screened["flag_method"] = ""
    screened["flag_group"] = ""
    screened["flag_group_n"] = np.nan
    screened["flag_group_mean"] = np.nan
    screened["flag_group_sd"] = np.nan
    screened["flag_statistic"] = np.nan
    screened["flag_p_value"] = np.nan
    screened["flag_alpha"] = OUTLIER_ALPHA
    screened["flag_reason"] = ""
    screened["manual_exclude"] = False
    screened["manual_exclusion_reason"] = ""

    for group_key, group in screened.groupby(OUTLIER_GROUP_COLS, observed=True):
        result = grubbs_two_sided_flag(group["actin_norm"])
        if result is None or result["p_value"] >= OUTLIER_ALPHA:
            continue

        idx = result["index"]
        group_label = ", ".join(
            f"{col}={value}" for col, value in zip(OUTLIER_GROUP_COLS, group_key)
        )
        screened.loc[idx, "flag_outlier"] = True
        screened.loc[idx, "flag_method"] = "two-sided Grubbs test"
        screened.loc[idx, "flag_group"] = group_label
        screened.loc[idx, "flag_group_n"] = result["n"]
        screened.loc[idx, "flag_group_mean"] = result["group_mean"]
        screened.loc[idx, "flag_group_sd"] = result["group_sd"]
        screened.loc[idx, "flag_statistic"] = result["g_statistic"]
        screened.loc[idx, "flag_p_value"] = result["p_value"]
        screened.loc[idx, "flag_reason"] = (
            "Most extreme value in line × dpi group; "
            f"two-sided Grubbs p = {result['p_value']:.4g} < {OUTLIER_ALPHA}"
        )

    for exclusion in MANUAL_EXCLUSIONS:
        mask = pd.Series(True, index=screened.index)
        for column in ["line", "dpi", "bio_rep"]:
            if column not in exclusion:
                raise ValueError(f"Manual exclusion missing required key: {column}")
            mask &= screened[column].astype(str) == str(exclusion[column])

        reason = exclusion.get("reason", "Manually confirmed exclusion")
        screened.loc[mask, "manual_exclude"] = True
        screened.loc[mask, "manual_exclusion_reason"] = reason

    review_log = screened[
        screened["flag_outlier"] | screened["manual_exclude"]
    ].copy()
    post = screened[~screened["manual_exclude"]].copy()
    return post, review_log


def errbars(values):
    sd = values.std(axis=1, ddof=1)
    return sd / np.sqrt(values.count(axis=1)) if ERRBAR == "SEM" else sd


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3, direction="out")
    apply_axis_style(ax)
    make_errorbars_black(ax)


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(
            SAVE_DIR / f"{output_stem}.{file_format}",
            dpi=300,
            bbox_inches="tight",
        )




def p_to_stars(p_value):
    if pd.isna(p_value) or p_value >= SIG_ALPHA:
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def contrast_label(contrast):
    left, right = contrast.split(" - ")
    return f"{SHORT_LINE_LABELS[left]}–{SHORT_LINE_LABELS[right]}"


def model_ready(data):
    ready = data.rename(columns={"actin_norm": "response"}).copy()
    ready["line"] = pd.Categorical(ready["line"], categories=LINES, ordered=False)
    ready["dpi"] = pd.Categorical(ready["dpi"].astype(str), ordered=False)
    return ready.dropna(subset=["line", "dpi", "response"])


def anova_rows(model, context):
    table = anova_lm(model, typ=3).reset_index(names="term")
    table["term"] = table["term"].replace(
        {
            "C(line, Sum)": "sample_type",
            "C(dpi, Sum)": "dpi",
            "C(line, Sum):C(dpi, Sum)": "sample_type:dpi",
        }
    )
    table = table[table["term"] != "Intercept"].copy()
    for key, value in context.items():
        table[key] = value
    return table


def add_fdr(table, group_cols):
    table = table.copy()
    table["p_bh_fdr"] = np.nan
    p_col = "PR(>F)"
    groupby_cols = group_cols[0] if len(group_cols) == 1 else group_cols
    for _, index in table.groupby(groupby_cols).groups.items():
        p_values = table.loc[index, p_col]
        valid = p_values.notna()
        if valid.any():
            _, adjusted, _, _ = multipletests(p_values[valid], method="fdr_bh")
            table.loc[p_values[valid].index, "p_bh_fdr"] = adjusted
    return table


def line_contrast_rows(model, context):
    """Pairwise sample-type contrasts from a one-way sample-type model."""
    design_info = model.model.data.design_info
    marginal_rows = {}

    for line in LINES:
        grid = pd.DataFrame(
            {"line": pd.Categorical([line], categories=LINES)}
        )
        design = np.asarray(build_design_matrices([design_info], grid)[0])
        marginal_rows[line] = design.mean(axis=0)

    rows = []
    p_values = []
    for left, right in combinations(LINES, 2):
        contrast = marginal_rows[left] - marginal_rows[right]
        test = model.t_test(contrast)
        rows.append(
            {
                **context,
                "contrast": f"{left} - {right}",
                "estimate": float(np.asarray(test.effect).ravel()[0]),
                "t": float(np.asarray(test.tvalue).ravel()[0]),
                "p_raw": float(np.asarray(test.pvalue).ravel()[0]),
            }
        )
        p_values.append(float(np.asarray(test.pvalue).ravel()[0]))

    if rows:
        _, adjusted, _, _ = multipletests(p_values, method="holm")
        for row, p_adj in zip(rows, adjusted):
            row["p_holm_within_family"] = float(p_adj)
            row["stars"] = p_to_stars(p_adj)

    return pd.DataFrame(rows)


def run_stats(cleaned):
    """Run global and dpi-specific stats; write reproducible stats tables."""
    ready = model_ready(cleaned)
    global_model = ols(
        "response ~ C(line, Sum) * C(dpi, Sum)",
        data=ready,
    ).fit()
    global_anova = add_fdr(
        anova_rows(global_model, {"analysis": "global_timecourse"}),
        ["term"],
    )

    point_anova = []
    point_contrasts = []

    for dpi, data in cleaned.groupby("dpi", observed=True):
        point_ready = model_ready(data)
        if point_ready["line"].nunique() < 2:
            continue

        model = ols("response ~ C(line, Sum)", data=point_ready).fit()
        context = {"analysis": "pointwise_dpi", "dpi": dpi}
        point_anova.append(anova_rows(model, context))
        point_contrasts.append(line_contrast_rows(model, context))

    point_anova = (
        add_fdr(pd.concat(point_anova, ignore_index=True), ["term"])
        if point_anova
        else pd.DataFrame()
    )
    point_contrasts = (
        pd.concat(point_contrasts, ignore_index=True)
        if point_contrasts
        else pd.DataFrame()
    )

    global_anova.to_excel(STATS_DIR / "rna_global_full_factorial_anova.xlsx", index=False)
    point_anova.to_excel(STATS_DIR / "rna_pointwise_anova.xlsx", index=False)
    point_contrasts.to_excel(
        STATS_DIR / "rna_pointwise_sample_type_pairwise_contrasts.xlsx",
        index=False,
    )

    summary = {
        "global_sample_type_fdr_lt_0.05": int(
            (
                (global_anova.get("term") == "sample_type")
                & (global_anova.get("p_bh_fdr") < 0.05)
            ).sum()
        ),
        "global_sample_type_by_dpi_fdr_lt_0.05": int(
            (
                (global_anova.get("term") == "sample_type:dpi")
                & (global_anova.get("p_bh_fdr") < 0.05)
            ).sum()
        ),
        "pointwise_sample_type_fdr_lt_0.05": int(
            (
                (point_anova.get("term") == "sample_type")
                & (point_anova.get("p_bh_fdr") < 0.05)
            ).sum()
        )
        if not point_anova.empty
        else 0,
    }
    return {
        "global_anova": global_anova,
        "point_anova": point_anova,
        "point_contrasts": point_contrasts,
        "summary": summary,
    }


def contrast_subset(contrasts, **filters):
    if contrasts is None or contrasts.empty:
        return pd.DataFrame()
    subset = contrasts.copy()
    for column, value in filters.items():
        subset = subset[subset[column] == value]
    return subset.copy()


def compact_letter_display(contrast_rows):
    """Assign compact letters from pairwise Holm-adjusted sample-type contrasts."""
    if contrast_rows is None or contrast_rows.empty:
        return {}, False

    p_lookup = {}
    any_significant = False
    for _, row in contrast_rows.iterrows():
        left, right = row["contrast"].split(" - ")
        p_value = row.get("p_holm_within_family", np.nan)
        significant = pd.notna(p_value) and p_value < SIG_ALPHA
        p_lookup[frozenset([left, right])] = significant
        any_significant = any_significant or significant

    if not any_significant:
        return {line: "a" for line in LINES}, False

    letters = "abc"
    masks = range(1, 2 ** len(letters))
    candidates = []
    for mask_values in __import__("itertools").product(masks, repeat=len(LINES)):
        assignment = {
            line: {letters[i] for i in range(len(letters)) if mask & (1 << i)}
            for line, mask in zip(LINES, mask_values)
        }
        valid = True
        for left, right in combinations(LINES, 2):
            share_letter = bool(assignment[left] & assignment[right])
            significant = p_lookup.get(frozenset([left, right]), False)
            if significant and share_letter:
                valid = False
                break
            if not significant and not share_letter:
                valid = False
                break
        if not valid:
            continue
        used = set().union(*assignment.values())
        label_lengths = sum(len(v) for v in assignment.values())
        candidates.append((len(used), label_lengths, mask_values, assignment))

    if not candidates:
        return {line: letters[i] for i, line in enumerate(LINES)}, True

    _, _, _, best = sorted(candidates, key=lambda x: (x[0], x[1], x[2]))[0]
    return {line: "".join(sorted(best[line])) for line in LINES}, True


def rna_mean_error(data, line):
    values = (
        data[data["line"] == line]
        .pivot_table(index="dpi", columns="bio_rep", values="actin_norm")
        .reindex(DPIS)
    )
    return values.mean(axis=1), errbars(values)


def annotate_timecourse_letters(ax, contrasts, data, force_show_dpis=()):
    """Add stacked compact-letter labels above significant dpi families."""
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1
    label_gap_y = LETTER_LABEL_GAP_FRACTION * y_range
    label_step_y = LETTER_LABEL_STEP_FRACTION * y_range
    max_label_y = y_max

    for dpi in DPIS:
        rows = contrast_subset(contrasts, dpi=dpi)
        letters, should_show = compact_letter_display(rows)
        if not should_show and dpi not in force_show_dpis:
            continue
        if not letters and dpi in force_show_dpis:
            letters = {line: "a" for line in LINES}

        line_positions = []
        for line in LINES:
            mean, error = rna_mean_error(data, line)
            if dpi not in mean.index or pd.isna(mean.loc[dpi]):
                continue
            point_values = data[(data["line"] == line) & (data["dpi"] == dpi)]["actin_norm"]
            error_extent = mean.loc[dpi] + (0 if pd.isna(error.loc[dpi]) else error.loc[dpi])
            point_extent = point_values.max()
            data_extent = error_extent if pd.isna(point_extent) else max(point_extent, error_extent)
            line_positions.append(
                {
                    "line": line,
                    "mean": mean.loc[dpi],
                    "data_extent": data_extent,
                }
            )

        if not line_positions:
            continue

        line_positions = sorted(line_positions, key=lambda item: item["mean"], reverse=True)
        stack_base_y = max(item["data_extent"] for item in line_positions) + label_gap_y
        for idx, item in enumerate(line_positions):
            y = stack_base_y + (len(line_positions) - idx - 1) * label_step_y
            line = item["line"]
            ax.text(
                dpi,
                y,
                letters[line],
                color=COLORS[line],
                fontsize=LETTER_FONT_SIZE,
                fontweight="bold",
                ha="center",
                va="center",
                zorder=6,
            )
            max_label_y = max(max_label_y, y + label_step_y)

    if max_label_y > y_max:
        ax.set_ylim(y_min, max_label_y + label_step_y)


def draw_timecourse(ax, data, response_col):
    """Draw mean ± error timecourse with individual biological points."""
    for line in LINES:
        values = (
            data[data["line"] == line]
            .pivot_table(index="dpi", columns="bio_rep", values=response_col)
            .reindex(DPIS)
        )
        mean = values.mean(axis=1)
        error = errbars(values)

        ax.fill_between(
            DPIS,
            mean - error,
            mean + error,
            color=COLORS[line],
            alpha=0.18,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            DPIS,
            mean,
            color=COLORS[line],
            marker=MARKERS[line],
            markerfacecolor=COLORS[line],
            markeredgecolor=MEAN_MARKER_EDGE,
            markeredgewidth=MEAN_MARKER_LINEWIDTH,
            markersize=6,
            linewidth=1.8,
            label=line,
            zorder=4,
        )

        for column in values.columns:
            ax.scatter(
                DPIS,
                values[column],
                marker=MARKERS[line],
                s=30,
                facecolors=COLORS[line],
                edgecolors="none",
                linewidths=0,
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=3,
            )

    ax.axhline(0, color="#bbbbbb", linewidth=0.6)
    ax.set_xticks(DPIS)
    ax.set_xlabel(XLABEL)
    style(ax)


def draw_grouped_bar_by_dpi(ax, data):
    """Draw grouped bars by dpi with individual biological points."""
    x = np.arange(len(DPIS))
    bar_width = 0.22
    offsets = {"WT": -bar_width, "PRp27#1": 0, "PRp27#2": bar_width}

    for line in LINES:
        values = (
            data[data["line"] == line]
            .pivot_table(index="dpi", columns="bio_rep", values="actin_norm")
            .reindex(DPIS)
        )
        mean = values.mean(axis=1)
        error = errbars(values)
        positions = x + offsets[line]

        ax.bar(
            positions,
            mean,
            width=bar_width * 0.85,
            color=COLORS[line],
            alpha=0.28,
            edgecolor=COLORS[line],
            linewidth=1.1,
            label=line,
            zorder=1,
        )
        ax.errorbar(
            positions,
            mean,
            yerr=error,
            color=ERRORBAR_COLOR,
            fmt="none",
            capsize=3,
            elinewidth=1.0,
            zorder=3,
        )

        for dpi_index, _dpi in enumerate(DPIS):
            point_values = values.iloc[dpi_index].dropna().to_numpy()
            if len(point_values) == 0:
                continue
            jitter = np.linspace(-bar_width * 0.22, bar_width * 0.22, len(point_values))
            ax.scatter(
                np.full(len(point_values), positions[dpi_index]) + jitter,
                point_values,
                marker=MARKERS[line],
                s=32,
                facecolors=COLORS[line],
                edgecolors="none",
                linewidths=0,
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=4,
            )

    ax.axhline(0, color="#bbbbbb", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([str(dpi) for dpi in DPIS])
    ax.set_xlabel(XLABEL)
    style(ax)


def annotate_grouped_bar_letters(ax, contrasts, data, force_show_dpis=()):
    """Add stacked compact-letter labels above grouped bar families."""
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1
    label_gap_y = LETTER_LABEL_GAP_FRACTION * y_range
    label_step_y = LETTER_LABEL_STEP_FRACTION * y_range
    max_label_y = y_max
    x = np.arange(len(DPIS))

    for dpi_index, dpi in enumerate(DPIS):
        rows = contrast_subset(contrasts, dpi=dpi)
        letters, should_show = compact_letter_display(rows)
        if not should_show and dpi not in force_show_dpis:
            continue
        if not letters and dpi in force_show_dpis:
            letters = {line: "a" for line in LINES}

        line_positions = []
        for line in LINES:
            values = data[(data["line"] == line) & (data["dpi"] == dpi)]["actin_norm"]
            if values.empty:
                continue
            mean = values.mean()
            error = values.sem() if ERRBAR == "SEM" else values.std(ddof=1)
            error_extent = mean + (0 if pd.isna(error) else error)
            point_extent = values.max()
            data_extent = error_extent if pd.isna(point_extent) else max(point_extent, error_extent)
            line_positions.append(
                {
                    "line": line,
                    "mean": mean,
                    "data_extent": data_extent,
                }
            )

        if not line_positions:
            continue

        line_positions = sorted(line_positions, key=lambda item: item["mean"], reverse=True)
        stack_base_y = max(item["data_extent"] for item in line_positions) + label_gap_y
        for idx, item in enumerate(line_positions):
            y = stack_base_y + (len(line_positions) - idx - 1) * label_step_y
            line = item["line"]
            ax.text(
                x[dpi_index],
                y,
                letters[line],
                color=COLORS[line],
                fontsize=LETTER_FONT_SIZE,
                fontweight="bold",
                ha="center",
                va="center",
                zorder=6,
            )
            max_label_y = max(max_label_y, y + label_step_y)

    if max_label_y > y_max:
        ax.set_ylim(y_min, max_label_y + label_step_y)


# %%
raw = load_raw(SRC)
clean, outlier_log = clean_data(raw)
stats_results = run_stats(clean)

raw.to_excel(STATS_DIR / "rna_raw_long.xlsx", index=False)
clean.to_excel(STATS_DIR / "rna_cleaned_long.xlsx", index=False)
outlier_log.to_excel(STATS_DIR / "rna_outlier_log.xlsx", index=False)

print(f"Raw RNA measurements before manual exclusion review: {len(raw)}")
print(f"Potential outlier flags for manual review: {int(outlier_log['flag_outlier'].sum()) if not outlier_log.empty else 0}")
print(f"Manually confirmed exclusions applied: {int(clean['manual_exclude'].sum()) if 'manual_exclude' in clean else 0}")
print(f"RNA measurements used for main plots/stats: {len(clean)}")
print(f"Log and stats Excel .xlsx files saved to: {STATS_DIR}")
print("Stats summary:", stats_results["summary"])
clean.head()


# %% [markdown]
# ## 0. Cleaning reference: pre-clean vs post-clean
#
# Left panel shows raw data; right panel shows post-clean data. Points are
# individual biological samples.

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 6.0), sharey=True)

for ax, title, dataset in zip(
    axes,
    ["Pre-clean raw data", "Post-clean data"],
    [raw, clean],
):
    draw_timecourse(ax, dataset, "actin_norm")
    ax.set_title(title, pad=10)

axes[0].set_ylabel(YLABEL)
axes[-1].legend(frameon=False, loc="upper right")
fig.suptitle(f"RNA cleaning reference (mean ± {ERRBAR})", fontsize=16, y=1.04)
fig.tight_layout(pad=1.8, w_pad=2.2, rect=[0, 0.32, 1, 1])
save(fig, "rna_cleaning_reference_pre_vs_post")
plt.show()


# %% [markdown]
# ## 1. RNA timecourse
#
# Final plot of actin-normalised RNA signal over time. Points are individual
# biological samples; shaded regions are mean ± SEM/SD.

# %%
fig, ax = plt.subplots(figsize=(8.0, 5.8))
draw_timecourse(ax, clean, "actin_norm")
annotate_timecourse_letters(ax, stats_results["point_contrasts"], clean)

ax.set_ylabel(YLABEL)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
fig.suptitle(f"RNA timecourse (mean ± {ERRBAR})", fontsize=16, y=1.04)
fig.tight_layout(pad=1.8, rect=[0, 0.32, 0.84, 1])
save(fig, "rna_timecourse")
plt.show()


# %% [markdown]
# ## 2. RNA grouped by dpi
#
# Bars show mean ± SEM/SD by sample type at each dpi. Points are individual
# biological samples.

# %%
fig, ax = plt.subplots(figsize=(9.2, 5.8))
draw_grouped_bar_by_dpi(ax, clean)
annotate_grouped_bar_letters(ax, stats_results["point_contrasts"], clean)

ax.set_ylabel(YLABEL)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
fig.suptitle(f"RNA by dpi (mean ± {ERRBAR})", fontsize=16, y=1.04)
fig.tight_layout(pad=1.8, rect=[0, 0.32, 0.84, 1])
save(fig, "rna_by_dpi_grouped_bars")
plt.show()


# %% [markdown]
# ---
#
# Figures are saved in `PRp27/PRp27_timecourse_202606/202606_RNA_figures/` and logs/statistical tables are saved in
# `PRp27/PRp27_timecourse_202606/202606_RNA_logs/`.
