# %% [markdown]
# # Colony count timecourse
#
# Run this file one `# %%` block at a time with VS Code's **Run Cell** action.
# The script reads the `ColonyCount` sheet from the original timecourse workbook,
# averages technical counts to biological-replicate means, exports cleaning/stat
# logs, and plots CFU/cm^2 leaf tissue over time.

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
    FONT_STACK,
    INDIVIDUAL_POINT_ALPHA,
    LETTER_FONT_SIZE,
    LETTER_LABEL_GAP_FRACTION,
    LETTER_LABEL_STEP_FRACTION,
    LOG_LETTER_GAP,
    LOG_LETTER_STEP,
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

SHEET = "ColonyCount"


SRC = resolve_timecourse_workbook(HERE, env_var="COLONY_SRC")
SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202606_ColonyCount"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
STATS_DIR = LOG_DIR
SAVE_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)
SAVE_FMTS = ["png", "pdf", "svg"]

LINES = ["WT", "PRp27#1", "PRp27#2"]
SHORT_LINE_LABELS = {"WT": "WT", "PRp27#1": "#1", "PRp27#2": "#2"}
DPIS = [0, 1, 3, 5]
MAIN_DPIS = [0, 1, 3]
BATCHES = [1, 2, 3]
ERRBAR = "SEM"  # "SD" or "SEM"
SIG_ALPHA = 0.05
OUTLIER_ALPHA = 0.01
OUTLIER_MIN_N = 6
YLABEL = r"CFU/cm$^2$ leaf tissue"
XLABEL = "dpi"
FONT_FAMILY = FONT_STACK

# Outlier handling: flag for manual review, do not remove automatically.
# Colony counts vary more than fluorescence/RNA readouts, so the default flagging
# is intentionally conservative: biological-replicate means are screened within
# sample type x dpi groups with alpha = 0.01, and only groups with at least 6
# biological-replicate means are tested.
MANUAL_EXCLUSIONS = [
    # Example after manual confirmation:
    # {"line": "PRp27#2", "batch": 2, "bio_rep": 1, "dpi": 3,
    #  "reason": "documented plating/counting failure"},
]

# --------------------------------------------------------------------------

plt.rcParams.update(ANALYSIS_RCPARAMS)

print(f"Setup complete. Using data file: {SRC}")
print(f"Figures will save to: {SAVE_DIR}")


# %%
def load_raw(src):
    """Load the ColonyCount main table and standardise column names/types."""
    df = pd.read_excel(src, sheet_name=SHEET, header=1, usecols="B:J")
    df = df.rename(
        columns={
            "Sample Type": "line",
            "Batch": "batch",
            "Tech Rep": "tech_rep",
            "Bio rep": "bio_rep",
            "dpi": "dpi",
            "Spot number": "spot_number",
            "Colony count": "colony_count",
            "Multiplication factor": "multiplication_factor",
            "Colony count (CFU/cm^2  leaf tissue)": "cfu_cm2",
        }
    )
    df = df.dropna(how="all").copy()

    numeric_cols = ["tech_rep", "dpi", "spot_number", "colony_count", "multiplication_factor", "cfu_cm2"]
    for column in numeric_cols:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    mq_controls = df[df["line"] == "MQ"].copy()
    analysis = df[df["line"].isin(LINES)].copy()
    analysis["line"] = pd.Categorical(analysis["line"], categories=LINES, ordered=False)
    analysis["batch"] = pd.to_numeric(analysis["batch"], errors="coerce").astype("Int64")
    analysis["bio_rep"] = pd.to_numeric(analysis["bio_rep"], errors="coerce").astype("Int64")
    analysis["tech_rep"] = analysis["tech_rep"].astype("Int64")
    analysis["dpi"] = analysis["dpi"].astype("Int64")
    analysis["spot_number"] = analysis["spot_number"].astype("Int64")

    return (
        analysis.dropna(subset=["line", "batch", "bio_rep", "dpi", "cfu_cm2"]),
        mq_controls,
    )


def biological_means(raw):
    """Average technical colony counts to one row per biological replicate."""
    grouped = (
        raw.groupby(["line", "batch", "bio_rep", "dpi"], observed=True)
        .agg(
            cfu_cm2=("cfu_cm2", "mean"),
            colony_count_mean=("colony_count", "mean"),
            colony_count_min=("colony_count", "min"),
            colony_count_max=("colony_count", "max"),
            n_tech=("tech_rep", "count"),
            spot_number_min=("spot_number", "min"),
            spot_number_max=("spot_number", "max"),
        )
        .reset_index()
    )
    grouped["line"] = pd.Categorical(grouped["line"], categories=LINES, ordered=False)
    grouped["dpi"] = grouped["dpi"].astype(int)
    grouped["batch"] = grouped["batch"].astype(int)
    grouped["bio_rep"] = grouped["bio_rep"].astype(int)
    grouped["log10_cfu_cm2"] = np.log10(grouped["cfu_cm2"] + 1)
    return grouped


def grubbs_two_sided_flag(values):
    """Return the most extreme value and two-sided Grubbs-test p-value."""
    values = values.dropna()
    n = len(values)
    if n < OUTLIER_MIN_N:
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


def clean_data(bio_means):
    """Flag potential outliers and remove only manually confirmed exclusions."""
    screened = bio_means.copy()
    screened["flag_outlier"] = False
    screened["flag_method"] = ""
    screened["flag_group"] = ""
    screened["flag_group_n"] = np.nan
    screened["flag_group_mean"] = np.nan
    screened["flag_group_sd"] = np.nan
    screened["flag_statistic"] = np.nan
    screened["flag_p_value"] = np.nan
    screened["flag_alpha"] = OUTLIER_ALPHA
    screened["manual_exclude"] = False
    screened["manual_exclusion_reason"] = ""

    for group_key, group in screened.groupby(["line", "dpi"], observed=True):
        result = grubbs_two_sided_flag(group["log10_cfu_cm2"])
        if result is None or result["p_value"] >= OUTLIER_ALPHA:
            continue

        idx = result["index"]
        group_label = f"line={group_key[0]}, dpi={group_key[1]}"
        screened.loc[idx, "flag_outlier"] = True
        screened.loc[idx, "flag_method"] = "two-sided Grubbs test on log10 biological means"
        screened.loc[idx, "flag_group"] = group_label
        screened.loc[idx, "flag_group_n"] = result["n"]
        screened.loc[idx, "flag_group_mean"] = result["group_mean"]
        screened.loc[idx, "flag_group_sd"] = result["group_sd"]
        screened.loc[idx, "flag_statistic"] = result["g_statistic"]
        screened.loc[idx, "flag_p_value"] = result["p_value"]

    for exclusion in MANUAL_EXCLUSIONS:
        mask = pd.Series(True, index=screened.index)
        for column in ["line", "batch", "bio_rep", "dpi"]:
            if column not in exclusion:
                raise ValueError(f"Manual exclusion missing required key: {column}")
            mask &= screened[column].astype(str) == str(exclusion[column])

        reason = exclusion.get("reason", "Manually confirmed exclusion")
        screened.loc[mask, "manual_exclude"] = True
        screened.loc[mask, "manual_exclusion_reason"] = reason

    review_log = screened[
        screened["flag_outlier"] | screened["manual_exclude"]
    ].copy()
    clean = screened[~screened["manual_exclude"]].copy()
    return clean, review_log, screened


def errbars(values):
    """Calculate row-wise SD or SEM."""
    sd = values.std(axis=1, ddof=1)
    return sd / np.sqrt(values.count(axis=1)) if ERRBAR == "SEM" else sd


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3)
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


def add_fdr(table, group_cols):
    if table.empty:
        return table
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


def colony_ready(data):
    ready = data.copy()
    ready["line"] = pd.Categorical(ready["line"], categories=LINES, ordered=False)
    ready["dpi"] = pd.Categorical(ready["dpi"].astype(str), ordered=False)
    ready["batch"] = pd.Categorical(ready["batch"].astype(str), ordered=False)
    ready["response"] = ready["log10_cfu_cm2"]
    return ready.dropna(subset=["line", "dpi", "batch", "response"])


def clean_anova_table(table):
    return (
        table.reset_index()
        .rename(columns={"index": "term"})
        .replace(
            {
                "C(line, Sum)": "sample_type",
                "C(dpi, Sum)": "dpi",
                "C(batch, Sum)": "batch",
                "C(line, Sum):C(dpi, Sum)": "sample_type:dpi",
            }
        )
    )


def line_contrast_rows(model, ready, context, include_batch=True):
    design_info = model.model.data.design_info
    present_lines = [line for line in LINES if line in set(ready["line"].astype(str))]
    if include_batch:
        batch_levels = list(ready["batch"].cat.categories)
    else:
        batch_levels = [None]

    marginal_rows = {}
    for line in present_lines:
        new_data = {
            "line": pd.Categorical([line] * len(batch_levels), categories=LINES),
        }
        if include_batch:
            new_data["batch"] = pd.Categorical(batch_levels, categories=batch_levels)
        design = build_design_matrices([design_info], pd.DataFrame(new_data))[0]
        marginal_rows[line] = design.mean(axis=0)

    rows = []
    p_values = []
    for left, right in combinations(present_lines, 2):
        contrast = marginal_rows[left] - marginal_rows[right]
        test = model.t_test(contrast)
        row = {
            **context,
            "contrast": f"{left} - {right}",
            "estimate_log10": float(np.asarray(test.effect).ravel()[0]),
            "t": float(np.asarray(test.tvalue).ravel()[0]),
            "p_raw": float(np.asarray(test.pvalue).ravel()[0]),
        }
        rows.append(row)
        p_values.append(row["p_raw"])

    if rows:
        _, adjusted, _, _ = multipletests(p_values, method="holm")
        for row, p_adj in zip(rows, adjusted):
            row["p_holm_within_family"] = float(p_adj)
            row["stars"] = p_to_stars(p_adj)

    return pd.DataFrame(rows)


def run_stats(clean):
    """Export global and pointwise ANOVA/contrast tables."""
    ready_all = colony_ready(clean)
    global_ready = ready_all[ready_all["dpi"].isin(["0", "1", "3"])].copy()

    global_model = ols(
        "response ~ C(line, Sum) * C(dpi, Sum) + C(batch, Sum)",
        data=global_ready,
    ).fit()
    global_anova = clean_anova_table(anova_lm(global_model, typ=2))
    global_anova["family"] = "global_0_1_3_dpi"
    global_anova = add_fdr(global_anova, ["family"])

    point_anova_rows = []
    point_contrasts = []
    for dpi in DPIS:
        point_ready = ready_all[ready_all["dpi"] == str(dpi)].copy()
        present_lines = point_ready["line"].dropna().astype(str).nunique()
        if point_ready.empty or present_lines < 2:
            continue

        include_batch = point_ready["batch"].nunique() > 1
        formula = (
            "response ~ C(line, Sum) + C(batch, Sum)"
            if include_batch
            else "response ~ C(line, Sum)"
        )
        model = ols(formula, data=point_ready).fit()
        point_anova = clean_anova_table(anova_lm(model, typ=2))
        point_anova["dpi"] = dpi
        point_anova["model"] = formula
        point_anova_rows.append(point_anova)

        point_contrasts.append(
            line_contrast_rows(
                model,
                point_ready,
                {"dpi": dpi, "model": formula},
                include_batch=include_batch,
            )
        )

    point_anova = (
        pd.concat(point_anova_rows, ignore_index=True)
        if point_anova_rows
        else pd.DataFrame()
    )
    point_anova = add_fdr(point_anova, ["dpi"]) if not point_anova.empty else point_anova
    point_contrasts = (
        pd.concat(point_contrasts, ignore_index=True)
        if point_contrasts
        else pd.DataFrame()
    )

    global_anova.to_excel(STATS_DIR / "colony_global_anova_0_1_3dpi.xlsx", index=False)
    point_anova.to_excel(STATS_DIR / "colony_pointwise_anova.xlsx", index=False)
    point_contrasts.to_excel(
        STATS_DIR / "colony_pointwise_sample_type_pairwise_contrasts.xlsx",
        index=False,
    )

    return {
        "global_anova": global_anova,
        "point_anova": point_anova,
        "point_contrasts": point_contrasts,
        "summary": {
            "global_sample_type_fdr_lt_0.05": int(
                (
                    (global_anova["term"] == "sample_type")
                    & (global_anova["p_bh_fdr"] < SIG_ALPHA)
                ).sum()
            ),
            "global_sample_type_by_dpi_fdr_lt_0.05": int(
                (
                    (global_anova["term"] == "sample_type:dpi")
                    & (global_anova["p_bh_fdr"] < SIG_ALPHA)
                ).sum()
            ),
            "pointwise_sample_type_fdr_lt_0.05": int(
                (
                    (point_anova.get("term", pd.Series(dtype=str)) == "sample_type")
                    & (point_anova.get("p_bh_fdr", pd.Series(dtype=float)) < SIG_ALPHA)
                ).sum()
            )
            if not point_anova.empty
            else 0,
        },
    }


def contrast_subset(contrasts, **filters):
    if contrasts is None or contrasts.empty:
        return pd.DataFrame()
    subset = contrasts.copy()
    for column, value in filters.items():
        subset = subset[subset[column] == value]
    return subset.copy()


def compact_letter_display(contrast_rows, present_lines):
    """Assign compact letters from pairwise Holm-adjusted sample-type contrasts."""
    present_lines = [line for line in LINES if line in present_lines]
    if contrast_rows is None or contrast_rows.empty or len(present_lines) < 2:
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
        return {line: "a" for line in present_lines}, False

    letters = "abc"
    masks = range(1, 2 ** len(letters))
    candidates = []
    for mask_values in __import__("itertools").product(masks, repeat=len(present_lines)):
        assignment = {
            line: {letters[i] for i in range(len(letters)) if mask & (1 << i)}
            for line, mask in zip(present_lines, mask_values)
        }
        valid = True
        for left, right in combinations(present_lines, 2):
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
        return {line: letters[i] for i, line in enumerate(present_lines)}, True

    _, _, _, best = sorted(candidates, key=lambda x: (x[0], x[1], x[2]))[0]
    return {line: "".join(sorted(best[line])) for line in present_lines}, True


def draw_timecourse(ax, data, dpi_values=MAIN_DPIS, lines=LINES):
    """Draw mean ± error colony-count timecourse with biological points."""
    x = np.asarray(dpi_values, dtype=float)
    positive_values = []
    for line in lines:
        values = (
            data[data["line"] == line]
            .pivot_table(index="dpi", columns=["batch", "bio_rep"], values="cfu_cm2")
            .reindex(dpi_values)
        )
        mean = values.mean(axis=1).astype(float).to_numpy()
        error = errbars(values).astype(float).to_numpy()
        valid = ~np.isnan(mean)
        positive_values.extend(values.to_numpy(dtype=float)[values.to_numpy(dtype=float) > 0])

        if valid.any():
            fill_valid = valid & ~np.isnan(error)
            ax.fill_between(
                x[fill_valid],
                np.maximum(mean[fill_valid] - error[fill_valid], 1),
                mean[fill_valid] + error[fill_valid],
                color=COLORS[line],
                alpha=0.18,
                linewidth=0,
                zorder=1,
            )
            ax.plot(
                x[valid],
                mean[valid],
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
            point_values = values[column].astype(float).to_numpy()
            point_valid = ~np.isnan(point_values)
            ax.scatter(
                x[point_valid],
                point_values[point_valid],
                marker=MARKERS[line],
                s=24,
                facecolors=COLORS[line],
                edgecolors="none",
                linewidths=0,
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=3,
            )

    if positive_values:
        min_positive = min(positive_values)
        max_positive = max(positive_values)
        ax.set_yscale("log")
        ax.set_ylim(min_positive / 3, max_positive * 4)
    ax.set_xticks(dpi_values)
    ax.set_xlabel(XLABEL)
    style(ax)


def draw_batch_timecourse(ax, data, batch, dpi_values=MAIN_DPIS):
    """Draw one batch diagnostic timecourse."""
    panel = data[data["batch"] == batch]
    draw_timecourse(ax, panel, dpi_values=dpi_values)
    ax.set_title(f"Batch {batch}", pad=10)


def annotate_timecourse_letters(
    ax,
    contrasts,
    data,
    dpi_values=MAIN_DPIS,
    lines=LINES,
    log_stack_gap=LOG_LETTER_GAP,
    log_stack_step=LOG_LETTER_STEP,
    force_show_dpis=(),
):
    """Add stacked compact-letter labels above significant dpi families."""
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1
    label_gap_y = LETTER_LABEL_GAP_FRACTION * y_range
    label_step_y = LETTER_LABEL_STEP_FRACTION * y_range
    max_label_y = y_max

    for dpi in dpi_values:
        rows = contrast_subset(contrasts, dpi=dpi)
        present_lines = [
            line
            for line in data[data["dpi"] == dpi]["line"].dropna().astype(str).unique()
            if line in lines
        ]
        letters, should_show = compact_letter_display(rows, present_lines)
        if not should_show and dpi not in force_show_dpis:
            continue
        if not letters and dpi in force_show_dpis:
            letters = {line: "a" for line in present_lines}

        line_positions = []
        for line in lines:
            values = data[(data["line"] == line) & (data["dpi"] == dpi)]["cfu_cm2"]
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
        stack_base_y = max(item["data_extent"] for item in line_positions)
        if ax.get_yscale() == "log":
            stack_base_y *= log_stack_gap
        else:
            stack_base_y += label_gap_y
        for idx, item in enumerate(line_positions):
            if ax.get_yscale() == "log":
                y = stack_base_y * (log_stack_step ** (len(line_positions) - idx - 1))
            else:
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
                zorder=7,
            )
            if ax.get_yscale() == "log":
                max_label_y = max(max_label_y, y * log_stack_step)
            else:
                max_label_y = max(max_label_y, y + label_step_y)

    if max_label_y > y_max:
        if ax.get_yscale() == "log":
            ax.set_ylim(y_min, max_label_y)
        else:
            ax.set_ylim(y_min, max_label_y + label_step_y)


# %%
raw, mq_controls = load_raw(SRC)
bio_mean_pre_clean = biological_means(raw)
clean, outlier_log, screened = clean_data(bio_mean_pre_clean)
stats_results = run_stats(clean)

raw.to_excel(STATS_DIR / "colony_raw_technical_rows.xlsx", index=False)
mq_controls.to_excel(STATS_DIR / "colony_mq_water_control_rows.xlsx", index=False)
bio_mean_pre_clean.to_excel(STATS_DIR / "colony_biological_means_pre_clean.xlsx", index=False)
screened.to_excel(STATS_DIR / "colony_outlier_screened_all_biological_means.xlsx", index=False)
clean.to_excel(STATS_DIR / "colony_cleaned_biological_means.xlsx", index=False)
outlier_log.to_excel(STATS_DIR / "colony_outlier_log.xlsx", index=False)

print(f"Raw technical colony-count rows before MQ/control removal: {len(raw) + len(mq_controls)}")
print(f"MQ water-control rows excluded from analysis plots/stats: {len(mq_controls)}")
print(f"Technical colony-count rows used for biological means: {len(raw)}")
print(f"Biological-replicate means before manual exclusion review: {len(bio_mean_pre_clean)}")
print(f"Potential outlier flags for manual review: {int(outlier_log['flag_outlier'].sum()) if not outlier_log.empty else 0}")
print(f"Manually confirmed exclusions applied: {int(screened['manual_exclude'].sum()) if 'manual_exclude' in screened else 0}")
print(f"Biological-replicate means used for plots/stats: {len(clean)}")
print(f"Log and stats Excel .xlsx files saved to: {STATS_DIR}")
print("Stats summary:", stats_results["summary"])
clean.head()


# %% [markdown]
# ## 0. Cleaning reference: pre-clean vs post-clean
#
# Biological-replicate means are shown before and after any manually confirmed
# exclusions. By default, no flagged values are removed automatically.

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), sharey=True)

for ax, title, dataset in zip(
    axes,
    ["Pre-clean biological means", "Post-clean biological means"],
    [bio_mean_pre_clean, clean],
):
    draw_timecourse(ax, dataset)
    ax.set_title(title, pad=10)

axes[0].set_ylabel(YLABEL)
axes[-1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
fig.suptitle(f"Colony count cleaning reference (mean ± {ERRBAR})", fontsize=16, y=1.04)
fig.subplots_adjust(left=0.08, right=0.84, top=0.84, bottom=0.34, wspace=0.35)
save(fig, "colony_cleaning_reference_pre_vs_post")
plt.show()


# %% [markdown]
# ## 1. Batch-separated diagnostics
#
# Rows show biological-replicate means before and after any manually confirmed
# exclusions; columns show independent batches.

# %%
fig, axes = plt.subplots(
    2,
    len(BATCHES),
    figsize=(14.5, 8.6),
    sharex=True,
    sharey=True,
)

for row, (row_label, dataset) in enumerate(
    [("Pre-clean", bio_mean_pre_clean), ("Post-clean", clean)]
):
    for col, batch in enumerate(BATCHES):
        ax = axes[row, col]
        draw_batch_timecourse(ax, dataset, batch)
        if col == 0:
            ax.set_ylabel(f"{row_label}\n{YLABEL}")
        if row == 0:
            ax.set_xlabel("")

axes[0, -1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
fig.suptitle(f"Colony count by batch (mean ± {ERRBAR})", fontsize=16, y=1.02)
fig.subplots_adjust(left=0.08, right=0.86, top=0.86, bottom=0.30, wspace=0.35, hspace=0.42)
save(fig, "colony_batch_diagnostics_pre_vs_post")
plt.show()


# %% [markdown]
# ## 2. Final colony-count timecourse
#
# Lines show sample-type means across biological-replicate means. Shaded regions
# are mean ± SEM/SD; points are individual biological-replicate means after
# technical averaging and any manually confirmed exclusions.

# %%
fig, ax = plt.subplots(figsize=(8.0, 5.8))
draw_timecourse(ax, clean, dpi_values=MAIN_DPIS)
annotate_timecourse_letters(ax, stats_results["point_contrasts"], clean, dpi_values=MAIN_DPIS)

ax.set_ylabel(YLABEL)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
fig.suptitle(f"Colony count timecourse to 3 dpi (mean ± {ERRBAR})", fontsize=16, y=1.04)
fig.subplots_adjust(left=0.14, right=0.80, top=0.84, bottom=0.34)
save(fig, "colony_final_timecourse")
plt.show()


# %% [markdown]
# ## 3. Complete 5 dpi sets
#
# Only line/batch sets with biological means at every dpi in 0, 1, 3, and 5 are
# shown here.

# %%
complete_set_keys = []
for (line, batch), group in clean.groupby(["line", "batch"], observed=True):
    if set(group["dpi"].astype(int)) >= set(DPIS):
        complete_set_keys.append((line, batch))

complete_set = pd.concat(
    [
        clean[(clean["line"] == line) & (clean["batch"] == batch)]
        for line, batch in complete_set_keys
    ],
    ignore_index=True,
) if complete_set_keys else clean.iloc[0:0].copy()
complete_lines = [line for line in LINES if line in complete_set["line"].astype(str).unique()]

fig, ax = plt.subplots(figsize=(8.0, 5.8))
draw_timecourse(ax, complete_set, dpi_values=DPIS, lines=complete_lines)
annotate_timecourse_letters(
    ax,
    stats_results["point_contrasts"],
    complete_set,
    dpi_values=DPIS,
    lines=complete_lines,
)

ax.set_ylabel(YLABEL)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
fig.suptitle(f"Complete 5 dpi colony-count sets (mean ± {ERRBAR})", fontsize=16, y=1.04)
fig.subplots_adjust(left=0.14, right=0.80, top=0.84, bottom=0.34)
save(fig, "colony_complete_5dpi_sets")
plt.show()


# %% [markdown]
# ---
#
# Figures are saved in `PRp27/PRp27_timecourse_202606/202606_ColonyCount_figures/` and logs/statistical
# tables are saved in `PRp27/PRp27_timecourse_202606/202606_ColonyCount_logs/`.
