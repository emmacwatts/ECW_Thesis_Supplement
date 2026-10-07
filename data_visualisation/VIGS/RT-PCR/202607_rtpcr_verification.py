# %% [markdown]
# # RT-PCR VIGS verification plots
#
# Parse the manually arranged RT-PCR verification workbook and make one compiled
# panel figure showing WT versus matched VIGS line band-intensity means.

# %%
from pathlib import Path
import os
import re
import sys

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
else:
    HERE = Path.cwd() / "VIGS" / "RT-PCR"

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))

from figure_styles import (  # noqa: E402
    ANALYSIS_RCPARAMS,
    ERRORBAR_COLOR,
    FONT_STACK,
    TEXT_COLOR,
    apply_axis_style,
)

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_rtpcr_verification"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
for folder in [SAVE_DIR, LOG_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "4. Immunity Silencing/VIGS/FurtherVIGS/rt-PCR_verification/"
    "RT-PCR-202606/RT-PCR_verif.xlsx"
)
NEW_SRC = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "4. Immunity Silencing/VIGS/FurtherVIGS/rt-PCR_verification/202609/202609_RTVIGS.xlsx"
)

SAVE_FMTS = ["png", "svg"]
SIG_ALPHA = 0.05
WT_COLOR = "#6D6D6D"
VIGS_COLOR = "#3E6E86"
BAR_ALPHA = 0.72
POINT_ALPHA = 0.82
BAR_WIDTH = 0.58

TARGET_DISPLAY = {
    "EDS5": "EDS5",
    "EAH": "EAH",
    "EAS": "EAS",
    "EIN3": "EIN3",
    "ETR1": "ETR1",
    "JAZ": "JAZ",
    "LYM3": "LYM3",
    "PEROX72": "PEROX72",
    "PYL8": "PYL8",
    "SQS": "SQS",
}

PANEL_ORDER = ["EAH", "EAS", "EDS5", "EIN3", "ETR1", "JAZ", "LYM3", "PEROX72", "PYL8", "SQS"]

STATS_METHODS = (
    "Stats/processing: original target-band data were parsed from RT-PCR_verif.xlsx and "
    "Actin plus optimized target-band data from 202609_RTVIGS.xlsx. Every target-band AUC "
    "was divided by the matching sample line's mean of three Actin AUC readings; GUS was "
    "interpreted as the control corresponding to the original WT group. Optimized readings "
    "replace LYM14/NbL14g16690.1, LYM04/NbL04g07070.1, EIN3/"
    "NbL07g03540.1 + NbL17g06220.1, both PYL8 "
    "assays (PYL03/NbL03g12610.1 and PYL8_03/NbL13g12750.1), SQS17/NbL17g01610.1 only, "
    "JAZ/NbL16g20890.1, and EAS/NbL17g23470.1. The EIN3 assay is recorded as targeting "
    "both supplied EIN3 genes; the first supplied ID was assigned to the single JAZ assay "
    "because the workbook does not distinguish its paralogues. "
    "Original readings are retained for PEROX72, EAH, ETR1, EDS5, and SQS09, and are also "
    "Actin-normalized. Optimized blocks with two target readings remain n=2. "
    "Original source data were parsed from repeated tables in RT-PCR_verif.xlsx "
    "using rows headed Sample Type, Gene ID, and Band intensity (AUC). Side notes and "
    "non-numeric band-intensity cells were ignored. Each target-gene comparison was "
    "summarised as mean +/- SEM from biological replicate band-intensity AUC values, "
    "with individual replicate values overlaid. For targets with multiple homologous "
    "genes/primer sets, WT and target VIGS bars are plotted as separate paired groups "
    "and the Gene ID is shown under the corresponding group; single-gene panels omit "
    "the Gene ID label. WT was compared with the matching TRV2::target line using a "
    "directional Welch two-sample t-test for each target-gene set, testing the "
    "pre-specified verification hypothesis that target VIGS band intensity is lower "
    "than WT. This directional test was used as the primary/plotted comparison "
    "because the biological question is knockdown rather than any difference in "
    "either direction, while Welch's test retains unequal-variance handling for low "
    "replicate counts. Two-sided Welch and two-sided equal-variance t-test p-values "
    "are also saved as sensitivity checks. Raw primary p-values, Holm-adjusted "
    "p-values within each target, and Holm-adjusted p-values across all target-gene "
    "tests are saved in the stats table; the plotted significance stars use the "
    "within-target Holm-adjusted primary p-values. No outlier removal or "
    "transformation or outlier removal was applied beyond Actin normalization."
)

REPLACEMENT_BLOCKS = [
    {"target": "EIN3", "gene_key": "NbL07g03540.1, NbL17g06220.1", "start": 4, "end": 7},
    {"target": "LYM3", "gene_key": "NbL14g16690.1", "start": 17, "end": 20},
    {"target": "PYL8", "gene_key": "NbL03g12610.1", "start": 23, "end": 26},
    {"target": "SQS", "gene_key": "NbL17g01610.1", "start": 29, "end": 34},
    {"target": "JAZ", "gene_key": "NbL16g20890.1", "start": 37, "end": 42},
    {"target": "LYM3", "gene_key": "NbL04g07070.1", "start": 45, "end": 48},
    {"target": "PYL8", "gene_key": "NbL13g12750.1", "start": 50, "end": 55},
    {"target": "EAS", "gene_key": "NbL17g23470.1", "start": 58, "end": 61},
]
ACTIN_LABEL_FOR_TARGET = {
    "LYM3": "LYM3", "PEROX72": "PEROX72", "EAH": "EAH", "EIN3": "EIN3",
    "PYL8": "PYL8", "ETR1": "ETR1", "SQS": "SQS", "EAS": "EAS",
    "JAZ": "JAZ", "EDS5": "EDS5",
}

RTPCR_RCPARAMS = {
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 32.5,
    "axes.titlesize": 34,
    "axes.labelsize": 32.5,
    "xtick.labelsize": 31.5,
    "ytick.labelsize": 31.5,
    "legend.fontsize": 32,
    "figure.dpi": 120,
    "savefig.dpi": 300,
}
plt.rcParams.update(RTPCR_RCPARAMS)


# %%
def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def clean_text(value):
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def is_numeric(value):
    return pd.notna(pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0])


def clean_target_key(value):
    text = clean_text(value)
    text = re.sub(r"\([^)]*\)", "", text)
    text = re.sub(r"^TRV2::", "", text, flags=re.IGNORECASE)
    text = text.replace("Peroxidase-72", "PEROX72")
    text = text.replace("Peroxidase 72", "PEROX72")
    text = re.sub(r"[^A-Za-z0-9]+", "", text).upper()
    if text in {"PEROXIDASE72", "PEROX"}:
        return "PEROX72"
    return text


def display_target(value):
    key = clean_target_key(value)
    return TARGET_DISPLAY.get(key, key)


def valid_gene_id(value):
    text = clean_text(value)
    if not text:
        return ""
    return text if re.search(r"[A-Za-z]", text) else ""


def is_sample_label(value):
    text = clean_text(value).upper()
    return text == "WT" or text.startswith("TRV2::")


def find_table_title(sheet, header_row):
    for row in range(header_row - 1, 0, -1):
        value = sheet.cell(row, 2).value
        if value is None:
            continue
        text = clean_text(value)
        if not text or is_sample_label(text) or text.lower() == "gels":
            continue
        if text.lower() == "sample type":
            continue
        return text
    return f"table_{header_row}"


def table_ranges(sheet):
    headers = []
    for row in range(1, sheet.max_row + 1):
        if clean_text(sheet.cell(row, 2).value).lower() != "sample type":
            continue
        if clean_text(sheet.cell(row, 4).value).lower() != "band intensity (auc)":
            continue
        headers.append(row)

    for i, header in enumerate(headers):
        end = headers[i + 1] - 1 if i + 1 < len(headers) else sheet.max_row
        yield header, end


def read_rtpcr_workbook(src):
    from openpyxl import load_workbook

    wb = load_workbook(src, data_only=True)
    sheet = wb.active
    rows = []

    for header_row, end_row in table_ranges(sheet):
        table_title = find_table_title(sheet, header_row)
        target_from_samples = ""
        raw_rows = []

        for row in range(header_row + 1, end_row + 1):
            sample_type = clean_text(sheet.cell(row, 2).value)
            auc_value = pd.to_numeric(pd.Series([sheet.cell(row, 4).value]), errors="coerce").iloc[0]
            if not is_sample_label(sample_type) or pd.isna(auc_value):
                continue
            if sample_type.upper().startswith("TRV2::"):
                target_from_samples = target_from_samples or clean_target_key(sample_type)
            raw_rows.append(
                {
                    "source_sheet": sheet.title,
                    "source_row": row,
                    "table_header_row": header_row,
                    "table_title": table_title,
                    "sample_type_raw": sample_type,
                    "gene_id": valid_gene_id(sheet.cell(row, 3).value),
                    "band_intensity_auc": float(auc_value),
                }
            )

        target_key = target_from_samples or clean_target_key(table_title)
        for raw_row in raw_rows:
            is_wt = raw_row["sample_type_raw"].upper() == "WT"
            raw_row["target"] = target_key
            raw_row["target_display"] = display_target(target_key)
            raw_row["group"] = "WT" if is_wt else f"TRV2::{display_target(target_key)}"
            rows.append(raw_row)

    data = pd.DataFrame(rows)
    if data.empty:
        raise ValueError(f"No RT-PCR data rows were parsed from {src}")

    data["gene_key"] = data["gene_id"].replace("", pd.NA)
    data["gene_key"] = data.groupby(["target", "table_header_row"])["gene_key"].transform(
        lambda values: values.ffill().bfill()
    )
    data["gene_key"] = data["gene_key"].fillna("single_gene")
    data["bio_rep"] = data.groupby(["target", "gene_key", "group"]).cumcount() + 1
    data["target_order"] = data["target"].map({target: i for i, target in enumerate(PANEL_ORDER)}).fillna(999)
    first_gene_order = (
        data[["target", "gene_key"]]
        .drop_duplicates()
        .assign(gene_order=lambda df: df.groupby("target").cumcount())
    )
    data = data.merge(first_gene_order, on=["target", "gene_key"], how="left")
    data["group_order"] = np.where(data["group"].eq("WT"), 0, 1)
    return data.sort_values(
        ["target_order", "target", "gene_order", "group_order", "bio_rep"]
    ).reset_index(drop=True)


def read_actin(src):
    from openpyxl import load_workbook
    sheet = load_workbook(src, data_only=True).active
    raw = pd.DataFrame([
        {"actin_label": clean_text(sheet.cell(r, 1).value).upper(), "source_row": r,
         "actin_auc": float(sheet.cell(r, 2).value)}
        for r in range(3, 36)
        if clean_text(sheet.cell(r, 1).value) and is_numeric(sheet.cell(r, 2).value)
    ])
    summary = raw.groupby("actin_label", sort=False).agg(
        actin_n=("actin_auc", "size"), actin_mean_auc=("actin_auc", "mean"),
        actin_sd_auc=("actin_auc", lambda x: x.std(ddof=1))).reset_index()
    if raw.empty or not summary["actin_n"].eq(3).all():
        raise ValueError("The Actin block must contain exactly three numeric readings per line")
    return raw, summary


def read_optimized(src):
    from openpyxl import load_workbook
    sheet = load_workbook(src, data_only=True).active
    rows = []
    for block in REPLACEMENT_BLOCKS:
        for r in range(block["start"], block["end"] + 1):
            label = clean_text(sheet.cell(r, 5).value).upper()
            value = sheet.cell(r, 6).value
            if not label or not is_numeric(value):
                continue
            control = label == "GUS"
            rows.append({
                "source_file": str(src), "source_sheet": sheet.title, "source_row": r,
                "table_header_row": pd.NA, "table_title": "202609 optimized conditions",
                "sample_type_raw": label, "gene_id": block["gene_key"],
                "band_intensity_auc": float(value), "target": block["target"],
                "target_display": display_target(block["target"]),
                "group": "WT" if control else f"TRV2::{display_target(block['target'])}",
                "gene_key": block["gene_key"], "data_status": "optimized replacement",
                "workbook_assay_label": label,
                "actin_label": "GUS" if control else ACTIN_LABEL_FOR_TARGET[block["target"]],
            })
    return pd.DataFrame(rows)


def combine_and_normalize(original, optimized, actin_summary):
    original = original.copy()
    original["source_file"] = str(SRC)
    original["data_status"] = "original retained"
    original["workbook_assay_label"] = original["sample_type_raw"]
    original["actin_label"] = np.where(
        original["group"].eq("WT"), "GUS", original["target"].map(ACTIN_LABEL_FOR_TARGET))
    keys = optimized[["target", "gene_key"]].drop_duplicates()
    marked = original.merge(keys.assign(_replace=True), on=["target", "gene_key"], how="left")
    # Legacy EIN3, JAZ, and PYL8 rows lacked gene IDs; their whole old series is replaced.
    replace_unidentified = marked["target"].isin(["EIN3", "JAZ", "PYL8"]) & marked["gene_key"].eq("single_gene")
    retained = marked.loc[marked["_replace"].isna() & ~replace_unidentified].drop(columns="_replace")
    data = pd.concat([retained, optimized], ignore_index=True, sort=False)
    data = data.merge(actin_summary, on="actin_label", how="left", validate="many_to_one")
    if data["actin_mean_auc"].isna().any():
        raise ValueError("A target lacks an Actin denominator")
    data["actin_normalized_signal"] = data["band_intensity_auc"] / data["actin_mean_auc"]
    data["bio_rep"] = data.groupby(["target", "gene_key", "group"]).cumcount() + 1
    data["target_order"] = data["target"].map({t: i for i, t in enumerate(PANEL_ORDER)}).fillna(999)
    order = data[["target", "gene_key"]].drop_duplicates().assign(
        gene_order=lambda x: x.groupby("target").cumcount())
    data = data.drop(columns="gene_order", errors="ignore").merge(order, on=["target", "gene_key"])
    data["group_order"] = np.where(data["group"].eq("WT"), 0, 1)
    return data.sort_values(["target_order", "target", "gene_order", "group_order", "bio_rep"]).reset_index(drop=True)


def sem(values):
    values = pd.Series(values).dropna()
    if len(values) <= 1:
        return 0.0
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def p_to_stars(p_value):
    if pd.isna(p_value) or p_value >= SIG_ALPHA:
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def run_stats(data):
    rows = []
    for (target, target_display, gene_key), subset in data.groupby(["target", "target_display", "gene_key"], sort=False):
        wt = subset.loc[subset["group"] == "WT", "actin_normalized_signal"].dropna()
        vigs_group = next((group for group in subset["group"].unique() if group != "WT"), "")
        vigs = subset.loc[subset["group"] == vigs_group, "actin_normalized_signal"].dropna()
        if len(wt) >= 2 and len(vigs) >= 2:
            primary_test = stats.ttest_ind(wt, vigs, equal_var=False, alternative="greater")
            welch_two_sided = stats.ttest_ind(wt, vigs, equal_var=False, alternative="two-sided")
            equal_var_two_sided = stats.ttest_ind(wt, vigs, equal_var=True, alternative="two-sided")
            statistic = float(primary_test.statistic)
            p_raw = float(primary_test.pvalue)
            welch_two_sided_statistic = float(welch_two_sided.statistic)
            p_welch_two_sided = float(welch_two_sided.pvalue)
            equal_var_two_sided_statistic = float(equal_var_two_sided.statistic)
            p_equal_var_two_sided = float(equal_var_two_sided.pvalue)
        else:
            statistic = np.nan
            p_raw = np.nan
            welch_two_sided_statistic = np.nan
            p_welch_two_sided = np.nan
            equal_var_two_sided_statistic = np.nan
            p_equal_var_two_sided = np.nan

        rows.append(
            {
                "target": target,
                "target_display": target_display,
                "gene_id": "" if gene_key == "single_gene" else gene_key,
                "comparison": f"WT vs {vigs_group}",
                "test": "one-sided Welch t-test; WT > target VIGS",
                "sensitivity_test_1": "two-sided Welch t-test",
                "sensitivity_test_2": "two-sided equal-variance t-test",
                "wt_n": int(len(wt)),
                "vigs_n": int(len(vigs)),
                "wt_mean": float(wt.mean()) if len(wt) else np.nan,
                "vigs_mean": float(vigs.mean()) if len(vigs) else np.nan,
                "wt_sem": sem(wt),
                "vigs_sem": sem(vigs),
                "estimate_vigs_minus_wt": float(vigs.mean() - wt.mean()) if len(wt) and len(vigs) else np.nan,
                "t_statistic": statistic,
                "p_raw": p_raw,
                "welch_two_sided_t_statistic": welch_two_sided_statistic,
                "p_welch_two_sided": p_welch_two_sided,
                "equal_var_two_sided_t_statistic": equal_var_two_sided_statistic,
                "p_equal_var_two_sided": p_equal_var_two_sided,
            }
        )

    stats_table = pd.DataFrame(rows)
    valid = stats_table["p_raw"].notna()
    stats_table["p_holm_within_target"] = np.nan
    stats_table["p_holm_all_tests"] = np.nan
    for _, index in stats_table.loc[valid].groupby("target").groups.items():
        _, adjusted, _, _ = multipletests(stats_table.loc[index, "p_raw"], method="holm")
        stats_table.loc[index, "p_holm_within_target"] = adjusted
    if valid.any():
        _, adjusted, _, _ = multipletests(stats_table.loc[valid, "p_raw"], method="holm")
        stats_table.loc[valid, "p_holm_all_tests"] = adjusted
    stats_table["stars"] = stats_table["p_holm_within_target"].map(p_to_stars)
    return stats_table


def summary_table(data):
    return (
        data.groupby(["target", "target_display", "gene_key", "group"], sort=False)
        .agg(
            n=("actin_normalized_signal", "size"),
            mean_actin_normalized=("actin_normalized_signal", "mean"),
            sem_actin_normalized=("actin_normalized_signal", sem),
            sd_actin_normalized=("actin_normalized_signal", lambda values: pd.Series(values).std(ddof=1)),
        )
        .reset_index()
        .assign(gene_id=lambda df: df["gene_key"].mask(df["gene_key"].eq("single_gene"), ""))
        .drop(columns=["gene_key"])
    )


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3, direction="out")
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def bar_positions_for_panel(n_genes):
    centers = np.arange(n_genes, dtype=float)
    offset = 0.18
    return centers - offset, centers + offset, centers


def format_gene_tick(gene_id):
    if str(gene_id) == "NbL07g03540.1, NbL17g06220.1":
        return "NbL07g03540.1\nNbL17g06220.1"
    return str(gene_id)


def ordered_gene_keys_for_plot(subset, stats_subset):
    gene_keys = list(dict.fromkeys(subset["gene_key"].tolist()))
    if len(gene_keys) <= 1 or stats_subset.empty:
        return gene_keys

    p_lookup = {}
    for _, row in stats_subset.iterrows():
        gene_key = row["gene_id"] if pd.notna(row["gene_id"]) and str(row["gene_id"]).strip() else "single_gene"
        p_lookup[str(gene_key)] = row.get("p_raw", np.nan)

    return sorted(
        gene_keys,
        key=lambda gene_key: (
            pd.isna(p_lookup.get(str(gene_key), np.nan)),
            p_lookup.get(str(gene_key), np.inf),
            gene_keys.index(gene_key),
        ),
    )


def plot_target_panel(ax, subset, stats_subset):
    target_display = subset["target_display"].iloc[0]
    gene_keys = ordered_gene_keys_for_plot(subset, stats_subset)
    wt_positions, vigs_positions, centers = bar_positions_for_panel(len(gene_keys))
    ymax = 0.0

    for idx, gene_key in enumerate(gene_keys):
        gene_data = subset.loc[subset["gene_key"] == gene_key]
        vigs_group = next((group for group in gene_data["group"].unique() if group != "WT"), f"TRV2::{target_display}")

        for group, x, color in [("WT", wt_positions[idx], WT_COLOR), (vigs_group, vigs_positions[idx], VIGS_COLOR)]:
            values = gene_data.loc[gene_data["group"] == group, "actin_normalized_signal"].dropna()
            if values.empty:
                continue
            mean = values.mean()
            err = sem(values)
            ymax = max(ymax, float(values.max()), float(mean + err))
            ax.bar(
                x,
                mean,
                width=BAR_WIDTH / 2,
                color=color,
                alpha=BAR_ALPHA,
                edgecolor="black",
                linewidth=1.0,
                zorder=2,
            )
            ax.errorbar(
                x,
                mean,
                yerr=err,
                fmt="none",
                ecolor=ERRORBAR_COLOR,
                elinewidth=1.0,
                capsize=2.7,
                capthick=1.0,
                zorder=5,
            )
            jitter = np.linspace(-0.045, 0.045, len(values)) if len(values) > 1 else np.array([0.0])
            ax.scatter(
                x + jitter,
                values,
                s=124,
                facecolor=color,
                edgecolor="#202020",
                linewidth=0.65,
                alpha=POINT_ALPHA,
                zorder=4,
            )

        stat_row = stats_subset.loc[stats_subset["gene_id"].fillna("").eq("" if gene_key == "single_gene" else gene_key)]
        if not stat_row.empty and stat_row["stars"].iloc[0]:
            gene_max = gene_data["actin_normalized_signal"].max()
            y = gene_max + max(ymax * 0.08, 0.04)
            ax.plot([wt_positions[idx], vigs_positions[idx]], [y, y], color=TEXT_COLOR, linewidth=0.8, clip_on=False)
            ax.text(
                centers[idx],
                y + max(ymax * 0.02, 0.01),
                stat_row["stars"].iloc[0],
                ha="center",
                va="bottom",
                fontsize=33.5,
                color=TEXT_COLOR,
                clip_on=False,
            )
            ymax = max(ymax, y + max(ymax * 0.10, 0.05))

    ax.set_title(target_display, pad=7)
    ax.set_xticks(centers)
    labels = [format_gene_tick(gene_key) for gene_key in gene_keys]
    tick_labels = ax.set_xticklabels(labels, rotation=30, ha="right")
    for label in tick_labels:
        label.set_fontsize(RTPCR_RCPARAMS["xtick.labelsize"] - 1.5)

    ax.set_ylim(0, ymax * 1.28 if ymax else 1)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    style_axis(ax)


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{output_stem}.{file_format}", dpi=300, bbox_inches="tight")


def plot_compiled(data, stats_table):
    # Five equal base columns allow two-gene panels to span two columns while
    # single-assay panels occupy half that width.
    panel_slots = {
        "EAH": (0, slice(0, 2)), "EAS": (0, 2), "EDS5": (0, 3), "EIN3": (0, 4),
        "ETR1": (1, slice(0, 2)), "JAZ": (1, 2), "LYM3": (1, slice(3, 5)),
        "PEROX72": (2, 0), "PYL8": (2, slice(1, 3)), "SQS": (2, slice(3, 5)),
    }
    fig = plt.figure(figsize=(24.0, 20.0), constrained_layout=False)
    grid = fig.add_gridspec(3, 5)
    axes_by_target = {}

    for target in PANEL_ORDER:
        ax = fig.add_subplot(grid[panel_slots[target]])
        axes_by_target[target] = ax
        subset = data.loc[data["target"] == target]
        stats_subset = stats_table.loc[stats_table["target"] == target]
        plot_target_panel(ax, subset, stats_subset)

    for target in ["EAH", "ETR1", "PEROX72"]:
        axes_by_target[target].set_ylabel("Actin-normalised\nsignal (AUC)")

    legend_handles = [
        Patch(facecolor=WT_COLOR, edgecolor="black", alpha=BAR_ALPHA, label="TRV2::GUS"),
        Patch(facecolor=VIGS_COLOR, edgecolor="black", alpha=BAR_ALPHA, label="target VIGS line"),
    ]
    legend = fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=2,
        frameon=False,
        handlelength=1.3,
        columnspacing=1.6,
    )
    for text in legend.get_texts():
        text.set_color(TEXT_COLOR)
        text.set_fontfamily(FONT_STACK)
        if text.get_text() == "TRV2::GUS":
            text.set_fontstyle("italic")

    fig.subplots_adjust(left=0.090, right=0.985, top=0.870, bottom=0.080, wspace=0.62, hspace=1.55)
    save(fig, "compiled_barplots")
    plt.close(fig)


def write_outputs(data, stats_table, actin_raw, actin_summary, replacements):
    summary = summary_table(data)
    data.to_excel(LOG_DIR / "rtpcr_cleaned_long.xlsx", index=False)
    summary.to_excel(LOG_DIR / "rtpcr_summary_means.xlsx", index=False)
    stats_table.to_excel(LOG_DIR / "rtpcr_welch_t_tests.xlsx", index=False)
    mapping = replacements.groupby(
        ["target", "gene_key", "data_status", "actin_label"], sort=False
    ).agg(
        source_rows=("source_row", lambda x: ", ".join(map(str, x))),
        assay_labels=("workbook_assay_label", lambda x: ", ".join(dict.fromkeys(x))),
        target_band_n=("band_intensity_auc", "size"),
    ).reset_index()
    with pd.ExcelWriter(LOG_DIR / "rtpcr_202609_interpretation_log.xlsx") as writer:
        actin_raw.to_excel(writer, sheet_name="actin_raw", index=False)
        actin_summary.to_excel(writer, sheet_name="actin_means", index=False)
        mapping.to_excel(writer, sheet_name="replacement_mapping", index=False)
        data.to_excel(writer, sheet_name="combined_normalized_data", index=False)
return summary


# %%
def main():
    original = read_rtpcr_workbook(SRC)
    actin_raw, actin_summary = read_actin(NEW_SRC)
    replacements = read_optimized(NEW_SRC)
    data = combine_and_normalize(original, replacements, actin_summary)
    stats_table = run_stats(data)
    summary = write_outputs(data, stats_table, actin_raw, actin_summary, replacements)
    plot_compiled(data, stats_table)
    print(f"Parsed rows: {len(data)}")
    print(f"Target-gene tests: {len(stats_table)}")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")
    print(
        stats_table[
            ["target_display", "gene_id", "wt_n", "vigs_n", "p_raw", "p_holm_within_target", "p_holm_all_tests", "stars"]
        ]
    )
    return data, summary, stats_table


if __name__ == "__main__":
    main()
