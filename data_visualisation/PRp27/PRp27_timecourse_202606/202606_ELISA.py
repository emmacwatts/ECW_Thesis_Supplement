# %% [markdown]
# # ELISA dilution curves — individual, merged, and across-timepoint plots
#
# If this page is not working - check the filepath and sheet names!
# Run this file one `# %%` block at a time with VS Code's **Run Cell** action.
# It reproduces the notebook workflow without embedding outputs in the file.
#
# ELISA exclusions are applied automatically because they are based on explicit
# technical-QC criteria: failed duplicate agreement and failed/background-level
# dilution curves. Every exclusion is written to the logs folder. If a row looks
# biologically plausible rather than technically failed, review the log and revise
# the pre-specified QC rule before using the plot for final reporting.

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
from patsy import build_design_matrices
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests
from figure_styles import (
    ANALYSIS_RCPARAMS,
    COLORS,
    DPI_COLORS,
    ERRORBAR_COLOR,
    FONT_STACK,
    INDIVIDUAL_POINT_ALPHA,
    LINESTYLES,
    MARKER_FACES,
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
SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/prp27/prp_TimeCourse2/COVA/202606_COVAELISA.xlsx"
)
SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202606_ELISA"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
STATS_DIR = LOG_DIR
SAVE_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)
SAVE_FMTS = ["png", "pdf", "svg"]
USE_LOGS_FOR_PLOTS = os.environ.get("ELISA_USE_LOGS", "0") == "1"

ERRBAR = "SEM"  # "SD" or "SEM"
YLABEL = r"Normalised A$_{450}$"
LINES = ["WT", "PRp27#1", "PRp27#2"]
DPIS = [0, 1, 3, 5]
DILS = [4, 40, 400, 4000]
FONT_FAMILY = FONT_STACK

# Workbook sheets
RAW_SHEET = "Orig_Raw"

# Cleaning toggles; these match the transparent in-script QC rules.
RULE_A = True
A_ABSDIFF = 0.10
A_CV = 20.0
RULE_B = True
B_BASELINE = 0.05
B_SIBLING = 0.15

BASE_DESCRIPTOR = (
    "Data source and cleaning: COVA ELISA data were read from the Orig_Raw sheet "
    "of the source workbook. Each measurement was background-normalised in the "
    "workbook, then reciprocal dilution was calculated as 1 / Dilution. Technical "
    "duplicate wells were averaged to give one value for each sample type, dpi, "
    "experimental repeat, biological replicate, and reciprocal dilution. Cleaning "
    "was performed in this script before plotting using explicit technical-QC "
    "rules rather than exploratory biological outlier trimming. Rule A removed "
    "individual measurements when the two technical duplicate wells disagreed by "
    f"more than {A_ABSDIFF:.2f} OD units and the raw-signal duplicate CV was greater "
    f"than {A_CV:.0f}%. Rule B removed complete biological-replicate dilution curves "
    f"at 3 or 5 dpi when the top-concentration measurement was below {B_BASELINE:.2f} "
    "while sibling biological replicates in the same sample type/dpi/repeat were "
    f"positive above {B_SIBLING:.2f}, consistent with a failed or background-level "
    "sample. These ELISA exclusions are applied automatically because they are tied "
    "to technical replicate disagreement or failed dilution-curve behaviour; the "
    "full exclusion log is exported to PRp27/PRp27_timecourse_202606/202606_ELISA_logs/. If a logged exclusion "
    "appears biologically plausible rather than technically failed, it should be "
    "reviewed manually and the pre-specified QC rule should be revised before final "
    "reporting. No additional hidden filtering is applied."
)

STATS_DESCRIPTOR = (
    f"Summary statistics: curves and bars show mean ± {ERRBAR}. SEM is used here "
    "to show uncertainty in the estimated mean; individual biological or repeat "
    "points are overlaid so the spread of the underlying observations remains "
    "visible. The plotted group means are calculated only from the post-clean "
    "dataset unless the panel explicitly states pre-clean. Inferential testing is "
    "exported to PRp27/PRp27_timecourse_202606/202606_ELISA_logs/. AUC statistics use biological-replicate "
    "AUC values from log10 reciprocal dilution and fit full-factorial OLS ANOVA "
    "models with sample type, dpi, and experimental repeat. Pointwise AUC tests "
    "compare sample types within each dpi while blocking by repeat. Dilution-curve "
    "pointwise tests compare sample types within each dpi and reciprocal dilution "
    "while blocking by repeat. Pairwise sample-type contrasts are Holm-adjusted "
    "within each dpi or dpi/dilution family, and ANOVA p-values also include "
    "Benjamini-Hochberg FDR-adjusted values. The intended biological comparison is "
    "between sample types within each dpi, not between days."
)

# --------------------------------------------------------------------------

plt.rcParams.update(ANALYSIS_RCPARAMS)


# %%
def load_raw_tech_means(src):
    """Average technical duplicates from the original raw sheet.

    This returns the pre-cleaning long format, so it can be compared with the
    in-script cleaned dataset.
    """
    df = pd.read_excel(src, sheet_name=RAW_SHEET)
    df["recip"] = (1 / df["Dilution"]).round().astype(int)
    keys = ["Sample Type", "dpi", "Repeat", "Replicate", "recip"]

    tech = df.pivot_table(
        index=keys,
        columns="Tech Rep",
        values="Normalised Signal",
    ).reset_index()
    tech.columns = [*keys, "t1", "t2"]

    tech["mean"] = tech[["t1", "t2"]].mean(axis=1)
    return tech[[*keys, "mean"]].copy()


def append_qc_reason(existing, new_reason):
    if existing:
        return f"{existing}; {new_reason}"
    return new_reason


def load_clean(src):
    """Average technical duplicates, apply documented QC removals, and log them."""
    df = pd.read_excel(src, sheet_name=RAW_SHEET)
    df["recip"] = (1 / df["Dilution"]).round().astype(int)
    keys = ["Sample Type", "dpi", "Repeat", "Replicate", "recip"]

    tech = df.pivot_table(
        index=keys,
        columns="Tech Rep",
        values="Normalised Signal",
    ).reset_index()
    tech.columns = [*keys, "t1", "t2"]

    signal = df.pivot_table(
        index=keys,
        columns="Tech Rep",
        values="Signal",
    ).reset_index()

    tech["mean"] = tech[["t1", "t2"]].mean(axis=1)
    tech["raw_signal_t1"] = signal[1].values
    tech["raw_signal_t2"] = signal[2].values
    tech["absdiff"] = (tech["t1"] - tech["t2"]).abs()
    tech["cv"] = (
        signal[[1, 2]].std(axis=1, ddof=0)
        / signal[[1, 2]].mean(axis=1)
        * 100
    ).values
    tech["drop"] = False
    tech["qc_rule"] = ""
    tech["qc_reason"] = ""
    tech["rule_a_absdiff_threshold"] = A_ABSDIFF
    tech["rule_a_cv_threshold"] = A_CV
    tech["rule_b_baseline_threshold"] = B_BASELINE
    tech["rule_b_sibling_threshold"] = B_SIBLING
    tech["rule_b_top_dilution_mean"] = np.nan
    tech["rule_b_sibling_max_top_dilution"] = np.nan

    if RULE_A:
        disagreement = (
            (tech["absdiff"] > A_ABSDIFF)
            & (tech["cv"] > A_CV)
        )
        reason = (
            f"Rule A: technical duplicate disagreement; absdiff > {A_ABSDIFF} "
            f"and raw-signal duplicate CV > {A_CV}%"
        )
        tech.loc[disagreement, "drop"] = True
        tech.loc[disagreement, "qc_rule"] = tech.loc[disagreement, "qc_rule"].apply(
            lambda x: append_qc_reason(x, "Rule A")
        )
        tech.loc[disagreement, "qc_reason"] = tech.loc[
            disagreement, "qc_reason"
        ].apply(lambda x: append_qc_reason(x, reason))

    if RULE_B:
        top_dilution = tech[tech["recip"] == 4]
        groups = top_dilution.groupby(["Sample Type", "dpi", "Repeat"])

        for (line, dpi, repeat), group in groups:
            sibling_max = group["mean"].max()
            if dpi not in (3, 5) or sibling_max <= B_SIBLING:
                continue

            for _, row in group.iterrows():
                if row["mean"] >= B_BASELINE:
                    continue

                dead_well = (
                    (tech["Sample Type"] == line)
                    & (tech["dpi"] == dpi)
                    & (tech["Repeat"] == repeat)
                    & (tech["Replicate"] == row["Replicate"])
                )
                reason = (
                    "Rule B: complete biological-replicate dilution curve removed; "
                    f"top-concentration mean {row['mean']:.4g} < {B_BASELINE} while "
                    f"sibling max {sibling_max:.4g} > {B_SIBLING}"
                )
                tech.loc[dead_well, "drop"] = True
                tech.loc[dead_well, "qc_rule"] = tech.loc[dead_well, "qc_rule"].apply(
                    lambda x: append_qc_reason(x, "Rule B")
                )
                tech.loc[dead_well, "qc_reason"] = tech.loc[
                    dead_well, "qc_reason"
                ].apply(lambda x: append_qc_reason(x, reason))
                tech.loc[dead_well, "rule_b_top_dilution_mean"] = row["mean"]
                tech.loc[dead_well, "rule_b_sibling_max_top_dilution"] = sibling_max

    qc_log = tech[tech["drop"]].copy()
    cleaned = tech[~tech["drop"]].copy()
    return cleaned, qc_log, tech


def repeat_means(cleaned):
    """Return one mean value per experimental repeat and dilution."""
    return (
        cleaned.groupby(
            ["Sample Type", "dpi", "Repeat", "recip"],
            as_index=False,
        )["mean"]
        .mean()
        .rename(columns={"mean": "val"})
    )


def curve_auc(cleaned):
    """Calculate one AUC per biological replicate over log10 dilution."""
    rows = []
    group_cols = ["Sample Type", "dpi", "Repeat", "Replicate"]

    for keys, group in cleaned.groupby(group_cols):
        group = group.sort_values("recip")
        if group["recip"].nunique() < 2:
            continue

        auc = np.trapezoid(
            group["mean"].to_numpy(),
            x=np.log10(group["recip"].to_numpy()),
        )
        rows.append(
            {
                "Sample Type": keys[0],
                "dpi": keys[1],
                "Repeat": keys[2],
                "Replicate": keys[3],
                "auc": auc,
            }
        )

    return pd.DataFrame(rows)


def model_ready(data, response_col):
    """Prepare ELISA data for statsmodels formulas."""
    ready = data.rename(
        columns={
            "Sample Type": "sample_type",
            "Repeat": "repeat",
            response_col: "response",
        }
    ).copy()
    ready["sample_type"] = pd.Categorical(
        ready["sample_type"], categories=LINES, ordered=False
    )
    ready["dpi"] = pd.Categorical(ready["dpi"].astype(str), ordered=False)
    ready["repeat"] = pd.Categorical(ready["repeat"].astype(str), ordered=False)
    if "recip" in ready:
        ready["recip"] = pd.Categorical(ready["recip"].astype(str), ordered=False)
    return ready.dropna(subset=["sample_type", "dpi", "repeat", "response"])


def anova_rows(model, context):
    table = anova_lm(model, typ=3).reset_index(names="term")
    table["term"] = table["term"].replace(
        {
            "C(sample_type, Sum)": "sample_type",
            "C(dpi, Sum)": "dpi",
            "C(repeat, Sum)": "repeat",
            "C(recip, Sum)": "reciprocal_dilution",
            "C(sample_type, Sum):C(dpi, Sum)": "sample_type:dpi",
            "C(sample_type, Sum):C(repeat, Sum)": "sample_type:repeat",
            "C(dpi, Sum):C(repeat, Sum)": "dpi:repeat",
            "C(sample_type, Sum):C(dpi, Sum):C(repeat, Sum)": (
                "sample_type:dpi:repeat"
            ),
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


def line_contrast_rows(model, ready, context):
    """Pairwise sample-type contrasts averaged over repeat blocks."""
    repeat_levels = list(ready["repeat"].cat.categories)
    design_info = model.model.data.design_info
    marginal_rows = {}

    for line in LINES:
        grid = pd.DataFrame(
            {
                "sample_type": pd.Categorical([line] * len(repeat_levels), categories=LINES),
                "repeat": pd.Categorical(repeat_levels, categories=repeat_levels),
            }
        )
        design = np.asarray(build_design_matrices([design_info], grid)[0])
        marginal_rows[line] = design.mean(axis=0)

    rows = []
    p_values = []
    for left, right in combinations(LINES, 2):
        contrast = marginal_rows[left] - marginal_rows[right]
        test = model.t_test(contrast)
        p_value = float(np.asarray(test.pvalue).ravel()[0])
        rows.append(
            {
                **context,
                "contrast": f"{left} - {right}",
                "estimate": float(np.asarray(test.effect).ravel()[0]),
                "t": float(np.asarray(test.tvalue).ravel()[0]),
                "p_raw": p_value,
            }
        )
        p_values.append(p_value)

    if rows:
        _, adjusted, _, _ = multipletests(p_values, method="holm")
        for row, p_adj in zip(rows, adjusted):
            row["p_holm_within_family"] = float(p_adj)
            row["stars"] = p_to_stars(p_adj)

    return pd.DataFrame(rows)


def p_to_stars(p_value):
    if pd.isna(p_value) or p_value >= 0.05:
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def run_elisa_stats(cleaned):
    """Export AUC-level and dilution-point model statistics."""
    auc = curve_auc(cleaned)
    auc.to_excel(STATS_DIR / "elisa_auc_values.xlsx", index=False)

    global_anova = []
    auc_point_anova = []
    auc_point_contrasts = []
    dilution_point_anova = []
    dilution_point_contrasts = []

    auc_ready = model_ready(auc, "auc")
    if (
        auc_ready["sample_type"].nunique() >= 2
        and auc_ready["dpi"].nunique() >= 2
        and auc_ready["repeat"].nunique() >= 2
    ):
        model = ols(
            "response ~ C(sample_type, Sum) * C(dpi, Sum) * C(repeat, Sum)",
            data=auc_ready,
        ).fit()
        global_anova.append(anova_rows(model, {"analysis": "auc_global"}))

    for dpi, data in auc.groupby("dpi"):
        ready = model_ready(data, "auc")
        if ready["sample_type"].nunique() < 2 or ready["repeat"].nunique() < 2:
            continue
        model = ols("response ~ C(sample_type, Sum) + C(repeat, Sum)", data=ready).fit()
        context = {"analysis": "auc_pointwise_dpi", "dpi": dpi}
        auc_point_anova.append(anova_rows(model, context))
        auc_point_contrasts.append(line_contrast_rows(model, ready, context))

    for (dpi, recip), data in cleaned.groupby(["dpi", "recip"]):
        ready = model_ready(data, "mean")
        if ready["sample_type"].nunique() < 2 or ready["repeat"].nunique() < 2:
            continue
        model = ols("response ~ C(sample_type, Sum) + C(repeat, Sum)", data=ready).fit()
        context = {
            "analysis": "dilution_pointwise_dpi_recip",
            "dpi": dpi,
            "reciprocal_dilution": recip,
        }
        dilution_point_anova.append(anova_rows(model, context))
        dilution_point_contrasts.append(line_contrast_rows(model, ready, context))

    tables = {
        "elisa_auc_global_full_factorial_anova.xlsx": add_fdr(
            pd.concat(global_anova, ignore_index=True), ["term"]
        )
        if global_anova
        else pd.DataFrame(),
        "elisa_auc_pointwise_anova.xlsx": add_fdr(
            pd.concat(auc_point_anova, ignore_index=True), ["term"]
        )
        if auc_point_anova
        else pd.DataFrame(),
        "elisa_auc_pointwise_sample_type_pairwise_contrasts.xlsx": pd.concat(
            auc_point_contrasts, ignore_index=True
        )
        if auc_point_contrasts
        else pd.DataFrame(),
        "elisa_dilution_pointwise_anova.xlsx": add_fdr(
            pd.concat(dilution_point_anova, ignore_index=True), ["term"]
        )
        if dilution_point_anova
        else pd.DataFrame(),
        "elisa_dilution_pointwise_sample_type_pairwise_contrasts.xlsx": pd.concat(
            dilution_point_contrasts, ignore_index=True
        )
        if dilution_point_contrasts
        else pd.DataFrame(),
    }

    for filename, table in tables.items():
        table.to_excel(STATS_DIR / filename, index=False)

    return {"auc": auc, **tables}


def errbars(values):
    """Calculate row-wise SD or SEM."""
    sd = values.std(axis=1, ddof=1)
    return sd / np.sqrt(values.count(axis=1)) if ERRBAR == "SEM" else sd


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3)
    apply_axis_style(ax)
    make_errorbars_black(ax)


def logx(ax):
    ax.set_xscale("log")
    ax.set_xticks(DILS)
    ax.set_xticklabels(DILS)


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(
            SAVE_DIR / f"{output_stem}.{file_format}",
            dpi=300,
            bbox_inches="tight",
        )




def draw_summary_curve(ax, data, dpi, line):
    """Draw mean ± error and biological points for one line at one dpi."""
    values = (
        data[(data["Sample Type"] == line) & (data["dpi"] == dpi)]
        .pivot_table(
            index="recip",
            columns=["Repeat", "Replicate"],
            values="mean",
        )
        .reindex(DILS)
    )
    mean = values.mean(axis=1)
    error = errbars(values)

    ax.fill_between(
        DILS,
        mean - error,
        mean + error,
        color=COLORS[line],
        alpha=0.18,
        linewidth=0,
        zorder=1,
    )
    ax.plot(
        DILS,
        mean,
        color=COLORS[line],
        linestyle=LINESTYLES[line],
        marker=MARKERS[line],
        markerfacecolor=MARKER_FACES[line],
        markeredgecolor=MEAN_MARKER_EDGE,
        markeredgewidth=MEAN_MARKER_LINEWIDTH,
        markersize=5,
        linewidth=1.35,
        label=line,
        zorder=4,
    )

    for column in values.columns:
        ax.scatter(
            DILS,
            values[column],
            marker=MARKERS[line],
            s=15,
            facecolors=MARKER_FACES[line],
            edgecolors="none",
            linewidths=0,
            alpha=INDIVIDUAL_POINT_ALPHA,
            zorder=3,
        )


def draw_auc_by_dpi(ax, auc_data):
    """Draw AUC grouped bars by dpi with biological-replicate points."""
    x = np.arange(len(DPIS))
    bar_width = 0.22
    offsets = {"WT": -bar_width, "PRp27#1": 0, "PRp27#2": bar_width}

    for line in LINES:
        line_data = auc_data[auc_data["Sample Type"] == line]
        means = line_data.groupby("dpi")["auc"].mean().reindex(DPIS)
        errors = (
            line_data.groupby("dpi")["auc"].sem().reindex(DPIS)
            if ERRBAR == "SEM"
            else line_data.groupby("dpi")["auc"].std(ddof=1).reindex(DPIS)
        )
        positions = x + offsets[line]

        ax.bar(
            positions,
            means,
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
            means,
            yerr=errors,
            color=ERRORBAR_COLOR,
            fmt="none",
            capsize=3,
            elinewidth=1.0,
            zorder=3,
        )

        for dpi_index, dpi in enumerate(DPIS):
            values = line_data[line_data["dpi"] == dpi]["auc"].to_numpy()
            if len(values) == 0:
                continue
            jitter = np.linspace(-bar_width * 0.22, bar_width * 0.22, len(values))
            ax.scatter(
                np.full(len(values), positions[dpi_index]) + jitter,
                values,
                marker=MARKERS[line],
                s=33,
                facecolors=MARKER_FACES[line],
                edgecolors="none",
                linewidths=0,
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=4,
            )

    ax.axhline(0, color="#bbbbbb", linewidth=0.6)
    style(ax)
    ax.set_xticks(x)
    ax.set_xticklabels([str(dpi) for dpi in DPIS])
    ax.set_xlabel("dpi")


# %%
if USE_LOGS_FOR_PLOTS:
    raw = pd.read_excel(STATS_DIR / "elisa_raw_tech_mean_pre_clean.xlsx")
    clean = pd.read_excel(STATS_DIR / "elisa_cleaned_long.xlsx")
    qc_screened = pd.read_excel(STATS_DIR / "elisa_qc_screened_all_rows.xlsx")
    qc_log = pd.read_excel(STATS_DIR / "elisa_qc_exclusion_log.xlsx")
    rm = repeat_means(clean)
else:
    raw = load_raw_tech_means(SRC)
    clean, qc_log, qc_screened = load_clean(SRC)
    rm = repeat_means(clean)
    run_elisa_stats(clean)

raw.to_excel(STATS_DIR / "elisa_raw_tech_mean_pre_clean.xlsx", index=False)
clean.to_excel(STATS_DIR / "elisa_cleaned_long.xlsx", index=False)
qc_screened.to_excel(STATS_DIR / "elisa_qc_screened_all_rows.xlsx", index=False)
qc_log.to_excel(STATS_DIR / "elisa_qc_exclusion_log.xlsx", index=False)

qc_summary = (
    qc_log.groupby(["qc_rule", "Sample Type", "dpi", "Repeat"], dropna=False)
    .size()
    .reset_index(name="n_rows_excluded")
)
qc_summary.to_excel(STATS_DIR / "elisa_qc_exclusion_summary.xlsx", index=False)

print(f"Raw tech-mean measurements before cleaning: {len(raw)}")
print(f"In-script cleaned measurements used for main plots: {len(clean)}")
print(f"Measurements excluded by in-script cleaning: {len(qc_log)}")
print(f"ELISA logs and stats Excel .xlsx files saved to: {STATS_DIR}")
clean.head()


# %% [markdown]
# ## 0. Cleaning reference: pre-clean vs post-clean
#
# The main plots below use the transparent cleaning rules in this script. This
# reference figure shows the raw tech-averaged data before outlier exclusion
# versus the cleaned data used for plotting.

# %%
fig, axes = plt.subplots(
    2,
    len(DPIS),
    figsize=(15, 9.2),
    sharex=True,
    sharey=True,
)

for row, (dataset_name, dataset) in enumerate(
    [("Pre-clean raw tech-mean", raw), ("Post-clean in-script data", clean)]
):
    for column, dpi in enumerate(DPIS):
        ax = axes[row, column]

        for line in LINES:
            draw_summary_curve(ax, dataset, dpi, line)

        ax.axhline(0, color="#bbbbbb", linewidth=0.7)
        logx(ax)
        style(ax)
        ax.set_title(f"{dpi} dpi", pad=10)

        if column == 0:
            ax.set_ylabel(f"{dataset_name}\n{YLABEL}")
        if row == 1:
            ax.set_xlabel("Reciprocal dilution")

axes[0, -1].legend(frameon=False, loc="upper right")
fig.suptitle(
    f"Cleaning reference: before vs after outlier exclusion (mean ± {ERRBAR})",
    fontsize=14,
    y=1.02,
)
fig.tight_layout(pad=1.8, w_pad=2.2, h_pad=2.4, rect=[0, 0.30, 1, 1])
save(fig, "cleaning_reference_pre_vs_post")
plt.show()


# %% [markdown]
# ## 1. Per timepoint × repeat
#
# Rows are dpi and columns are experimental repeats. Points are biological
# samples; the line and translucent zone show mean ± SD/SEM.

# %%
reps = sorted(clean["Repeat"].unique())

prism_color = COLORS
prism_line = LINESTYLES
prism_face = MARKER_FACES

fig, axes = plt.subplots(
    len(DPIS),
    len(reps),
    figsize=(3.2 * len(reps), 2.9 * len(DPIS)),
    sharex=True,
    sharey=True,
)

for row, dpi in enumerate(DPIS):
    for column, repeat in enumerate(reps):
        ax = axes[row, column]

        for line in LINES:
            values = (
                clean[
                    (clean["Sample Type"] == line)
                    & (clean["dpi"] == dpi)
                    & (clean["Repeat"] == repeat)
                ]
                .pivot_table(
                    index="recip",
                    columns="Replicate",
                    values="mean",
                )
                .reindex(DILS)
            )
            mean = values.mean(axis=1)
            error = errbars(values)

            ax.fill_between(
                DILS,
                mean - error,
                mean + error,
                color=prism_color[line],
                alpha=0.18,
                linewidth=0,
                zorder=1,
            )
            ax.plot(
                DILS,
                mean,
                color=prism_color[line],
                linestyle=prism_line[line],
                marker=MARKERS[line],
                markerfacecolor=prism_face[line],
                markeredgecolor=MEAN_MARKER_EDGE,
                markeredgewidth=MEAN_MARKER_LINEWIDTH,
                markersize=5,
                linewidth=1.35,
                label=line,
                zorder=4,
            )

            for replicate in values.columns:
                ax.scatter(
                    DILS,
                    values[replicate],
                    marker=MARKERS[line],
                    s=15,
                    facecolors=prism_face[line],
                    edgecolors="none",
                    linewidths=0,
                    alpha=INDIVIDUAL_POINT_ALPHA,
                    zorder=3,
                )

        ax.axhline(0, color="#b3b3b3", linewidth=0.6)
        logx(ax)
        style(ax)
        ax.tick_params(direction="out")

        if row == 0:
            ax.set_title(f"Repeat {repeat}", pad=10)
        if column == 0:
            ax.set_ylabel(f"{dpi} dpi\n{YLABEL}")
        if row == len(DPIS) - 1:
            ax.set_xlabel("Reciprocal dilution")

axes[0, -1].legend(
    frameon=False,
    handlelength=2.4,
    borderaxespad=0.4,
    loc="upper right",
)
fig.suptitle(f"Dilution curves (mean ± {ERRBAR})", fontsize=14, y=1.015)
fig.tight_layout(pad=1.8, w_pad=2.0, h_pad=2.4, rect=[0, 0.17, 1, 1])
save(fig, "curves_per_timepoint_repeat")
plt.show()


# %% [markdown]
# ### Optional: export each timepoint × repeat panel separately
#
# Change the switch below to `True` to write 12 standalone figures.

# %%
EXPORT_INDIVIDUAL_PANELS = False

if EXPORT_INDIVIDUAL_PANELS:
    for dpi in DPIS:
        for repeat in reps:
            fig, ax = plt.subplots(figsize=(4.5, 5.0))

            for line in LINES:
                values = (
                    clean[
                        (clean["Sample Type"] == line)
                        & (clean["dpi"] == dpi)
                        & (clean["Repeat"] == repeat)
                    ]
                    .pivot_table(
                        index="recip",
                        columns="Replicate",
                        values="mean",
                    )
                    .reindex(DILS)
                )
                ax.plot(
                    DILS,
                    values.mean(axis=1),
                    color=COLORS[line],
                    marker=MARKERS[line],
                    markersize=5,
                    linewidth=1.7,
                    label=line,
                )
                for replicate in values.columns:
                    ax.scatter(
                        DILS,
                        values[replicate],
                        color=COLORS[line],
                        s=14,
                        alpha=INDIVIDUAL_POINT_ALPHA,
                    )

            ax.axhline(0, color="#cccccc", linewidth=0.6)
            logx(ax)
            style(ax)
            ax.set_title(f"{dpi} dpi, repeat {repeat}")
            ax.set_xlabel("Reciprocal dilution")
            ax.set_ylabel(YLABEL)
            ax.legend(frameon=False)
            fig.tight_layout(pad=1.6, rect=[0, 0.40, 1, 1])
            save(fig, f"panel_dpi{dpi}_rep{repeat}")
            plt.close(fig)


# %% [markdown]
# ## 2. Per timepoint, repeats merged
#
# One panel per dpi. The line is the mean of the three experimental-repeat
# means with SD/SEM error bars (n = 3). Points are individual repeat means.

# %%
fig, axes = plt.subplots(1, 4, figsize=(15, 6.4), sharey=True)

for ax, dpi in zip(axes, DPIS):
    for line in LINES:
        values = (
            rm[(rm["Sample Type"] == line) & (rm["dpi"] == dpi)]
            .pivot_table(
                index="recip",
                columns="Repeat",
                values="val",
            )
            .reindex(DILS)
        )

        mean = values.mean(axis=1)
        error = errbars(values)

        ax.fill_between(
            DILS,
            mean - error,
            mean + error,
            color=COLORS[line],
            alpha=0.18,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            DILS,
            mean,
            color=COLORS[line],
            linestyle=LINESTYLES[line],
            marker=MARKERS[line],
            markerfacecolor=MARKER_FACES[line],
            markeredgecolor=MEAN_MARKER_EDGE,
            markeredgewidth=MEAN_MARKER_LINEWIDTH,
            markersize=5,
            linewidth=1.35,
            label=line,
            zorder=4,
        )

        for repeat in values.columns:
            ax.scatter(
                DILS,
                values[repeat],
                marker=MARKERS[line],
                s=15,
                facecolors=MARKER_FACES[line],
                edgecolors="none",
                linewidths=0,
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=3,
            )

    ax.axhline(0, color="#bbbbbb", linewidth=0.6)
    logx(ax)
    style(ax)
    ax.set_title(f"{dpi} dpi", pad=10)
    ax.set_xlabel("Reciprocal dilution")

axes[0].set_ylabel(YLABEL)
axes[-1].legend(frameon=False)
fig.suptitle(
    f"Dilution curves by day, repeats merged (mean ± {ERRBAR}, n = 3)",
    fontsize=14,
    y=1.06,
)
fig.tight_layout(pad=1.8, w_pad=2.2, rect=[0, 0.38, 1, 1])
save(fig, "curves_by_dpi_merged")
plt.show()


# %% [markdown]
# ## 3. AUC by dpi
#
# AUC summarises each biological replicate's dilution curve across the log10
# reciprocal-dilution series. Bars show mean ± SEM/SD and points show individual
# biological-replicate AUC values.

# %%
auc = curve_auc(clean)

fig, ax = plt.subplots(figsize=(8.5, 6.8))
draw_auc_by_dpi(ax, auc)
ax.set_ylabel(r"AUC of normalised A$_{450}$")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
fig.suptitle(f"COVA AUC (mean ± {ERRBAR})", fontsize=14, y=1.02)
fig.tight_layout(pad=1.8, rect=[0, 0.40, 0.86, 1])
save(fig, "auc_by_dpi")
plt.show()


# %% [markdown]
# ## 4. Across timepoints, by line
#
# One panel per line, with the dilution curve at each dpi overlaid
# (mean ± SD/SEM, n = 3).

# %%
dpi_colors = [DPI_COLORS[dpi] for dpi in DPIS]
fig, axes = plt.subplots(1, 3, figsize=(13, 6.4), sharey=True)

for ax, line in zip(axes, LINES):
    for dpi, color in zip(DPIS, dpi_colors):
        values = (
            rm[(rm["Sample Type"] == line) & (rm["dpi"] == dpi)]
            .pivot_table(
                index="recip",
                columns="Repeat",
                values="val",
            )
            .reindex(DILS)
        )

        mean = values.mean(axis=1)
        error = errbars(values)

        ax.fill_between(
            DILS,
            mean - error,
            mean + error,
            color=color,
            alpha=0.18,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            DILS,
            mean,
            color=color,
            marker="o",
            markerfacecolor=color,
            markeredgecolor=MEAN_MARKER_EDGE,
            markeredgewidth=MEAN_MARKER_LINEWIDTH,
            markersize=5,
            linewidth=1.35,
            label=f"{dpi} dpi",
            zorder=4,
        )

    ax.axhline(0, color="#bbbbbb", linewidth=0.6)
    logx(ax)
    style(ax)
    ax.set_title(line, pad=10)
    ax.set_xlabel("Reciprocal dilution")

axes[0].set_ylabel(YLABEL)
axes[-1].legend(frameon=False, title="dpi")
fig.suptitle(
    f"Dilution curves across timepoints, by line "
    f"(mean ± {ERRBAR}, n = 3)",
    fontsize=14,
    y=1.06,
)
fig.tight_layout(pad=1.8, w_pad=2.2, rect=[0, 0.38, 1, 1])
save(fig, "curves_across_timepoints")
plt.show()


# %% [markdown]
# ---
#
# Figures are saved in `PRp27/PRp27_timecourse_202606/202606_ELISA_figures/` as PNG, PDF, and SVG.
# Logs/statistical tables are saved in `PRp27/PRp27_timecourse_202606/202606_ELISA_logs/`. SVG and PDF
# retain editable text for later adjustment in Inkscape or Illustrator.
