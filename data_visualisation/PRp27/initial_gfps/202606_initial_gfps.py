# %% [markdown]
# # Initial PRp27 GFP endpoint datasets
#
# Endpoint GFP and virB plots for the PRp27 initial GFP project. The script reads
# the original Excel workbooks, writes raw/cleaned analysis tables and outlier
# review logs, runs endpoint statistics, and saves PNG/PDF/SVG figures.

# %%
from itertools import combinations
from pathlib import Path
import os
import sys

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.patches import FancyBboxPatch
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
elif (Path.cwd() / "PRp27" / "initial_gfps").is_dir():
    HERE = Path.cwd() / "PRp27" / "initial_gfps"
else:
    HERE = Path.cwd()

STYLE_DIR = HERE.parent / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))

from figure_styles import (  # noqa: E402
    ANALYSIS_RCPARAMS,
    ERRORBAR_COLOR,
    FONT_STACK,
    INDIVIDUAL_POINT_ALPHA,
    PALETTE_SEQUENCE,
    TEXT_COLOR,
    apply_axis_style,
)

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202606_initial_gfps"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
SAVE_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)
SAVE_FMTS = ["png", "pdf", "svg"]

SOURCE_FILES = {
    "20250505": Path(
        "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/prp27/GFP/20250505_6wk_GFP_slightlyDry/prp27-6wk-p19andnop19.xlsx"
    ),
    "20260330": Path(
        "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/prp27/GFP/20260330/20260330.xlsx"
    ),
    "20260415": Path(
        "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/prp27/GFP/20260415/gfp_prp1and2_check.xlsx"
    ),
    "20250708": Path(
        "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/prp27/GFP/prp_CORE_CSPR/20250728_cspr_core_prp27_nlsAndLMU.xlsx"
    ),
    "202405_4wk": Path(
        "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/InitialMutants_Batch/eds1_npr1_prp27_en_052024/mutantBatch_GFP.xlsx"
    ),
}

P19_ORDER = ["no p19", "p19"]
LEAF_ORDER = ["Lower", "Middle", "Upper"]
ERRBAR = "SEM"
OUTLIER_ALPHA = 0.05
BAR_ALPHA = 0.60
BAR_LINEWIDTH = 1.0
SIG_ALPHA = 0.05
EXCLUDED_LINES = ["core 1-2", "706-4", "CSPR"]
GFP_Y_SCALE = 1000
VIRB_Y_SCALE = 100
BAR_WIDTH = 0.34
BAR_CENTER_SPACING = 0.62
WT_REFERENCE_COLOR = "#C44E52"
PANEL_BAND_COLORS = ["#DCEFF8", "#FDE8C6", "#F2E2C0", "#F4EAF4", "#59A14F"]

LINE_ORDERS = {
    "20250505": ["WT", "PRp27#1"],
    "20260330": ["WT", "PRp27#1"],
    "20260415": ["WT", "PRp27#1", "PRp27#2"],
    "20250708": ["WT", "PRp27#1"],
    "202405_4wk": ["WT", "PRp27#1"],
}

COLORS = {
    "WT": PALETTE_SEQUENCE[0],
    "PRp27#1": PALETTE_SEQUENCE[1],
    "PRp27#2": PALETTE_SEQUENCE[2],
    "npr1": PALETTE_SEQUENCE[3],
    "en": PALETTE_SEQUENCE[4],
    "eds1": PALETTE_SEQUENCE[5],
    "CSPR": PALETTE_SEQUENCE[2],
    "706-4": PALETTE_SEQUENCE[3],
    "core 1-2": PALETTE_SEQUENCE[4],
}

MARKERS = {
    "WT": "o",
    "PRp27#1": "^",
    "PRp27#2": "s",
    "npr1": "D",
    "en": "P",
    "eds1": "X",
    "CSPR": "s",
    "706-4": "D",
    "core 1-2": "P",
}


INITIAL_RCPARAMS = {
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.title_fontsize": 8,
}

plt.rcParams.update(INITIAL_RCPARAMS)


# %%
def standard_line_label(value):
    text = str(value).strip()
    lookup = {
        "WT": "WT",
        "wt": "WT",
        "cas9 WT": "WT",
        "Cas9 WT": "WT",
        "prp27": "PRp27#1",
        "prp27-": "PRp27#1",
        "PRP27": "PRp27#1",
        "prp27 #1": "PRp27#1",
        "prp27 #2": "PRp27#2",
        "npr1": "npr1",
        "en": "en",
        "eds1": "eds1",
        "706-4-": "706-4",
        "706-4": "706-4",
        "cspr1-13-": "CSPR",
        "CSPR": "CSPR",
        "core 1-2": "core 1-2",
    }
    return lookup.get(text, text)


def load_20250505():
    df = pd.read_excel(SOURCE_FILES["20250505"], sheet_name="Sheet1", header=0)
    df["Plant Type"] = df["Plant Type"].ffill()
    return pd.DataFrame(
        {
            "dataset": "20250505",
            "assay": "GFP",
            "source_sheet": "Sheet1",
            "batch": "batch1",
            "condition_type": "p19",
            "condition": df["p19"].astype(str).str.strip(),
            "line_raw": df["Plant Type"],
            "line": df["Plant Type"].map(standard_line_label),
            "background": pd.to_numeric(df["background"], errors="coerce"),
            "signal": pd.to_numeric(df["gfp"], errors="coerce"),
            "response": pd.to_numeric(df["normalised gfp"], errors="coerce"),
        }
    )


def load_20260330():
    frames = []
    for sheet in ["batch1", "batch2"]:
        df = pd.read_excel(SOURCE_FILES["20260330"], sheet_name=sheet, header=1)
        frames.append(
            pd.DataFrame(
                {
                    "dataset": "20260330",
                    "assay": "GFP",
                    "source_sheet": sheet,
                    "batch": sheet,
                    "condition_type": "p19",
                    "condition": df["p19?"].astype(str).str.strip(),
                    "line_raw": df["Sample Type"],
                    "line": df["Sample Type"].map(standard_line_label),
                    "background": pd.to_numeric(df["Background"], errors="coerce"),
                    "signal": pd.to_numeric(df["GFP"], errors="coerce"),
                    "response": pd.to_numeric(df["Normalised GFP"], errors="coerce"),
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def load_20260415():
    df = pd.read_excel(SOURCE_FILES["20260415"], sheet_name="Sheet1", header=1)
    return pd.DataFrame(
        {
            "dataset": "20260415",
            "assay": "GFP",
            "source_sheet": "Sheet1",
            "batch": "batch1",
            "condition_type": "p19",
            "condition": df["Infiltration Type"].astype(str).str.strip(),
            "line_raw": df["Sample Type"],
            "line": df["Sample Type"].map(standard_line_label),
            "background": pd.to_numeric(df["Background"], errors="coerce"),
            "signal": pd.to_numeric(df["Signal"], errors="coerce"),
            "response": pd.to_numeric(df["Normalised Signal"], errors="coerce"),
        }
    )


def load_20250708_lmu():
    df = pd.read_excel(SOURCE_FILES["20250708"], sheet_name="LMU", header=2)
    df = df.dropna(subset=["Plant Type", "Leaf position", "Normalised GFP"])
    leaf_lookup = {"L": "Lower", "M": "Middle", "U": "Upper"}
    return pd.DataFrame(
        {
            "dataset": "20250708",
            "assay": "GFP",
            "source_sheet": "LMU",
            "batch": "batch1",
            "condition_type": "leaf_position",
            "condition": df["Leaf position"].map(leaf_lookup),
            "line_raw": df["Plant Type"],
            "line": df["Plant Type"].map(standard_line_label),
            "background": pd.to_numeric(df["Background"], errors="coerce"),
            "signal": pd.to_numeric(df["GFP"], errors="coerce"),
            "response": pd.to_numeric(df["Normalised GFP"], errors="coerce"),
        }
    )


def load_20250708_virb():
    df = pd.read_excel(SOURCE_FILES["20250708"], sheet_name="virB", header=1)
    df = df.dropna(subset=["Plant Type", "Normalised VirB"])
    return pd.DataFrame(
        {
            "dataset": "20250708",
            "assay": "virB",
            "source_sheet": "virB",
            "batch": "batch1",
            "condition_type": "virB",
            "condition": "virB",
            "line_raw": df["Plant Type"],
            "line": df["Plant Type"].map(standard_line_label),
            "background": pd.to_numeric(df["Background"], errors="coerce"),
            "signal": pd.to_numeric(df["virB Signal"], errors="coerce"),
            "response": pd.to_numeric(df["Normalised VirB"], errors="coerce"),
        }
    )


def load_202405_4wk():
    df = pd.read_excel(SOURCE_FILES["202405_4wk"], sheet_name="RawPlain", header=1)
    df = df.dropna(subset=["plant_type", "mean_gfp_intensity_background"])
    return pd.DataFrame(
        {
            "dataset": "202405_4wk",
            "assay": "GFP",
            "source_sheet": "RawPlain",
            "batch": "batch1",
            "condition_type": "plant_age",
            "condition": "4-week",
            "line_raw": df["plant_type"],
            "line": df["plant_type"].map(standard_line_label),
            "background": pd.to_numeric(df["Background"], errors="coerce"),
            "signal": pd.to_numeric(df["RawGFP"], errors="coerce"),
            "response": pd.to_numeric(
                df["mean_gfp_intensity_background"],
                errors="coerce",
            ),
        }
    )


def load_all():
    data = pd.concat(
        [
            load_20250505(),
            load_20260330(),
            load_20260415(),
            load_20250708_lmu(),
            load_20250708_virb(),
            load_202405_4wk(),
        ],
        ignore_index=True,
    )
    data = data.dropna(subset=["dataset", "assay", "condition", "line", "response"])
    data["bio_rep"] = (
        data.groupby(["dataset", "assay", "condition", "line", "batch"], observed=True)
        .cumcount()
        .add(1)
    )
    return data


def make_analysis_input(raw):
    analysis = raw[~raw["line"].isin(EXCLUDED_LINES)].copy()
    excluded = raw[raw["line"].isin(EXCLUDED_LINES)].copy()
    return analysis, excluded


def grubbs_two_sided_flag(values):
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
    screened = raw.copy()
    for column in [
        "flag_outlier",
        "flag_method",
        "flag_group",
        "flag_group_n",
        "flag_group_mean",
        "flag_group_sd",
        "flag_statistic",
        "flag_p_value",
        "flag_alpha",
        "manual_exclude",
        "manual_exclusion_reason",
    ]:
        screened[column] = False if column in ["flag_outlier", "manual_exclude"] else np.nan

    screened["flag_method"] = ""
    screened["flag_group"] = ""
    screened["manual_exclusion_reason"] = ""
    screened["flag_alpha"] = OUTLIER_ALPHA

    group_cols = ["dataset", "assay", "condition", "line", "batch"]
    for group_key, group in screened.groupby(group_cols, observed=True):
        result = grubbs_two_sided_flag(group["response"])
        if result is None or result["p_value"] >= OUTLIER_ALPHA:
            continue
        idx = result["index"]
        screened.loc[idx, "flag_outlier"] = True
        screened.loc[idx, "flag_method"] = "two-sided Grubbs test"
        screened.loc[idx, "flag_group"] = ", ".join(
            f"{col}={value}" for col, value in zip(group_cols, group_key)
        )
        screened.loc[idx, "flag_group_n"] = result["n"]
        screened.loc[idx, "flag_group_mean"] = result["group_mean"]
        screened.loc[idx, "flag_group_sd"] = result["group_sd"]
        screened.loc[idx, "flag_statistic"] = result["g_statistic"]
        screened.loc[idx, "flag_p_value"] = result["p_value"]

    cleaned = screened[~screened["manual_exclude"]].copy()
    review_log = screened[screened["flag_outlier"] | screened["manual_exclude"]].copy()
    return cleaned, screened, review_log


def sem(values):
    values = pd.Series(values).dropna()
    if len(values) <= 1:
        return 0.0
    return values.std(ddof=1) / np.sqrt(len(values))


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{output_stem}.{file_format}", dpi=300, bbox_inches="tight")


def clear_previous_outputs():
    for folder in [SAVE_DIR, LOG_DIR]:
        for path in folder.glob(f"{SCRIPT_STEM}_*"):
            if path.is_file():
                path.unlink()
    for path in LOG_DIR.glob("initial_gfps_*.xlsx"):
        if path.is_file():
            path.unlink()




def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3, direction="out")
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def style_age_comparison(ax):
    for bar in ax.patches:
        if isinstance(bar, Rectangle):
            bar.set_alpha(BAR_ALPHA)

    for collection in ax.collections:
        if not isinstance(collection, PathCollection):
            continue
        facecolors = collection.get_facecolors().copy()
        if facecolors.size:
            facecolors[:, :3] *= 0.72
            facecolors[:, 3] = 0.90
            collection.set_facecolors(facecolors)
        collection.set_alpha(0.90)
        collection.set_edgecolor("black")
        collection.set_linewidth(0.45)


def axis_ylabel(text):
    return text.replace("Background-normalised ", "Background-normalised\n")


def scaled_ylabel(text, scale):
    label = axis_ylabel(text)
    return f"{label} (x{scale})" if scale != 1 else label


def display_line_label(line):
    if line == "PRp27#1":
        return "prp27-1"
    if line == "PRp27#2":
        return "prp27-2"
    return line


def p_to_stars(p_value):
    if pd.isna(p_value):
        return ""
    if p_value >= SIG_ALPHA:
        return "ns"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def add_star_bracket(ax, x1, x2, y, h, label):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color=ERRORBAR_COLOR, linewidth=0.9)
    label_size = 12.5 if "*" in str(label) else 9.5
    ax.text(
        (x1 + x2) / 2,
        y + h * 1.35,
        label,
        ha="center",
        va="bottom",
        fontsize=label_size,
        color=TEXT_COLOR,
    )


def n_label_from_counts(counts):
    values = [count for count in counts.values() if count > 0]
    if not values:
        return ""
    low = min(values)
    high = max(values)
    return f"n={low}" if low == high else f"n={low}-{high}"


def add_axis_group_band(
    fig,
    axes,
    color,
    pad_x=0.016,
    pad_y=0.026,
    pad_left=None,
    pad_bottom=None,
    alpha=1.0,
):
    axes = np.atleast_1d(axes).ravel()
    positions = [ax.get_position() for ax in axes if ax.get_visible()]
    if not positions:
        return
    x0 = max(0, min(pos.x0 for pos in positions) - (pad_left if pad_left is not None else pad_x))
    y0 = max(0, min(pos.y0 for pos in positions) - (pad_bottom if pad_bottom is not None else pad_y))
    x1 = min(1, max(pos.x1 for pos in positions) + pad_x)
    y1 = min(1, max(pos.y1 for pos in positions) + pad_y)
    fig.add_artist(
        FancyBboxPatch(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            boxstyle="round,pad=0.004,rounding_size=0.012",
            transform=fig.transFigure,
            facecolor=color,
            alpha=alpha,
            edgecolor="none",
            zorder=-10,
        )
    )


def run_family_stats(data, family, line_order, include_batch=False):
    ready = data[data["line"].isin(line_order)].copy()
    ready["line"] = pd.Categorical(ready["line"], categories=line_order)
    ready["batch"] = pd.Categorical(ready["batch"])
    ready = ready.dropna(subset=["response", "line"])
    if ready["line"].nunique() < 2:
        return pd.DataFrame(), pd.DataFrame()

    formula = "response ~ C(line, Sum)"
    if include_batch and ready["batch"].nunique() > 1:
        formula += " + C(batch, Sum)"

    model = ols(formula, data=ready).fit()
    anova = anova_lm(model, typ=3).reset_index(names="term")
    anova = anova[anova["term"] != "Intercept"].copy()
    anova["term"] = anova["term"].replace(
        {"C(line, Sum)": "line", "C(batch, Sum)": "batch"}
    )
    for key, value in family.items():
        anova[key] = value

    rows = []
    p_values = []
    for left, right in combinations(line_order, 2):
        if "WT" not in (left, right):
            continue
        left_values = ready.loc[ready["line"] == left, "response"]
        right_values = ready.loc[ready["line"] == right, "response"]
        if len(left_values) < 2 or len(right_values) < 2:
            continue
        test = stats.ttest_ind(left_values, right_values, equal_var=False, nan_policy="omit")
        rows.append(
            {
                **family,
                "contrast": f"{left} - {right}",
                "left": left,
                "right": right,
                "left_mean": left_values.mean(),
                "right_mean": right_values.mean(),
                "estimate_left_minus_right": left_values.mean() - right_values.mean(),
                "t": float(test.statistic),
                "p_raw": float(test.pvalue),
            }
        )
        p_values.append(float(test.pvalue))

    contrasts = pd.DataFrame(rows)
    if not contrasts.empty:
        _, adjusted, _, _ = multipletests(p_values, method="holm")
        contrasts["p_holm_within_family"] = adjusted
        contrasts["stars"] = [p_to_stars(p) for p in adjusted]

    return anova, contrasts


def compact_letter_display(contrast_rows, line_order):
    if contrast_rows is None or contrast_rows.empty:
        return {line: "a" for line in line_order}

    significant = {}
    for _, row in contrast_rows.iterrows():
        significant[frozenset([row["left"], row["right"]])] = (
            row["p_holm_within_family"] < SIG_ALPHA
        )
    if not any(significant.values()):
        return {line: "a" for line in line_order}

    letters = "abcd"
    masks = range(1, 2 ** len(letters))
    best = None
    for mask_values in __import__("itertools").product(masks, repeat=len(line_order)):
        assignment = {
            line: {letters[i] for i in range(len(letters)) if mask & (1 << i)}
            for line, mask in zip(line_order, mask_values)
        }
        valid = True
        for left, right in combinations(line_order, 2):
            shares_letter = bool(assignment[left] & assignment[right])
            is_sig = significant.get(frozenset([left, right]), False)
            if is_sig and shares_letter:
                valid = False
                break
            if not is_sig and not shares_letter:
                valid = False
                break
        if not valid:
            continue
        used = sorted({letter for values in assignment.values() for letter in values})
        score = (len(used), sum(len(values) for values in assignment.values()))
        if best is None or score < best[0]:
            best = (score, assignment)
    if best is None:
        return {line: "" for line in line_order}
    return {line: "".join(sorted(best[1][line])) for line in line_order}


def run_stats(cleaned):
    anova_tables = []
    contrast_tables = []

    for dataset in ["20250505", "20260415"]:
        data = cleaned[(cleaned["dataset"] == dataset) & (cleaned["assay"] == "GFP")]
        for condition in P19_ORDER:
            family_data = data[data["condition"] == condition]
            family = {
                "dataset": dataset,
                "assay": "GFP",
                "condition_type": "p19",
                "condition": condition,
                "batch": "batch1",
            }
            anova, contrasts = run_family_stats(
                family_data,
                family,
                LINE_ORDERS[dataset],
            )
            anova_tables.append(anova)
            contrast_tables.append(contrasts)

    data_0330 = cleaned[(cleaned["dataset"] == "20260330") & (cleaned["assay"] == "GFP")]
    for batch in ["batch1", "batch2"]:
        for condition in P19_ORDER:
            family_data = data_0330[
                (data_0330["batch"] == batch) & (data_0330["condition"] == condition)
            ]
            family = {
                "dataset": "20260330",
                "assay": "GFP",
                "condition_type": "p19",
                "condition": condition,
                "batch": batch,
            }
            anova, contrasts = run_family_stats(
                family_data,
                family,
                LINE_ORDERS["20260330"],
            )
            anova_tables.append(anova)
            contrast_tables.append(contrasts)

    lmu = cleaned[(cleaned["dataset"] == "20250708") & (cleaned["assay"] == "GFP")]
    for condition in LEAF_ORDER:
        family_data = lmu[lmu["condition"] == condition]
        family = {
            "dataset": "20250708",
            "assay": "GFP",
            "condition_type": "leaf_position",
            "condition": condition,
            "batch": "batch1",
        }
        anova, contrasts = run_family_stats(family_data, family, LINE_ORDERS["20250708"])
        anova_tables.append(anova)
        contrast_tables.append(contrasts)

    virb = cleaned[(cleaned["dataset"] == "20250708") & (cleaned["assay"] == "virB")]
    family = {
        "dataset": "20250708",
        "assay": "virB",
        "condition_type": "virB",
        "condition": "virB",
        "batch": "batch1",
    }
    anova, contrasts = run_family_stats(virb, family, LINE_ORDERS["20250708"])
    anova_tables.append(anova)
    contrast_tables.append(contrasts)

    data_4wk = cleaned[
        (cleaned["dataset"] == "202405_4wk") & (cleaned["assay"] == "GFP")
    ]
    family = {
        "dataset": "202405_4wk",
        "assay": "GFP",
        "condition_type": "plant_age",
        "condition": "4-week",
        "batch": "batch1",
    }
    anova, contrasts = run_family_stats(data_4wk, family, LINE_ORDERS["202405_4wk"])
    anova_tables.append(anova)
    contrast_tables.append(contrasts)

    anova = pd.concat([x for x in anova_tables if not x.empty], ignore_index=True)
    contrasts = pd.concat([x for x in contrast_tables if not x.empty], ignore_index=True)
    if not anova.empty:
        anova["p_bh_fdr"] = np.nan
        for term, idx in anova.groupby("term").groups.items():
            p_values = anova.loc[idx, "PR(>F)"]
            valid = p_values.notna()
            if valid.any():
                _, adjusted, _, _ = multipletests(p_values[valid], method="fdr_bh")
                anova.loc[p_values[valid].index, "p_bh_fdr"] = adjusted

    anova.to_excel(LOG_DIR / "initial_gfps_endpoint_anova.xlsx", index=False)
    contrasts.to_excel(LOG_DIR / "initial_gfps_endpoint_pairwise_contrasts.xlsx", index=False)
    return anova, contrasts


def run_combination_checks(cleaned):
    rows = []

    data_0330 = cleaned[
        (cleaned["dataset"] == "20260330")
        & (cleaned["assay"] == "GFP")
        & (cleaned["line"].isin(LINE_ORDERS["20260330"]))
    ].copy()
    for condition, subset in data_0330.groupby("condition", observed=True):
        ready = subset.dropna(subset=["response", "line", "batch"]).copy()
        if ready["line"].nunique() < 2 or ready["batch"].nunique() < 2:
            continue
        ready["line"] = pd.Categorical(ready["line"], categories=LINE_ORDERS["20260330"])
        model = ols("response ~ C(line, Sum) * C(batch, Sum)", data=ready).fit()
        table = anova_lm(model, typ=3).reset_index(names="term")
        table = table[table["term"] != "Intercept"].copy()
        table["term"] = table["term"].replace(
            {
                "C(line, Sum)": "line",
                "C(batch, Sum)": "batch",
                "C(line, Sum):C(batch, Sum)": "line:batch",
            }
        )
        table["check"] = "20260330_batch_pooling"
        table["condition"] = condition
        rows.append(table)

    date_data = cleaned[
        (cleaned["assay"] == "GFP")
        & (cleaned["condition_type"] == "p19")
        & (cleaned["line"].isin(["WT", "PRp27#1"]))
        & (cleaned["dataset"].isin(["20250505", "20260330", "20260415"]))
    ].copy()
    for condition, subset in date_data.groupby("condition", observed=True):
        ready = subset.dropna(subset=["response", "line", "dataset"]).copy()
        if ready["line"].nunique() < 2 or ready["dataset"].nunique() < 2:
            continue
        ready["line"] = pd.Categorical(ready["line"], categories=["WT", "PRp27#1"])
        ready["dataset"] = pd.Categorical(
            ready["dataset"], categories=["20250505", "20260330", "20260415"]
        )
        model = ols("response ~ C(line, Sum) * C(dataset, Sum)", data=ready).fit()
        table = anova_lm(model, typ=3).reset_index(names="term")
        table = table[table["term"] != "Intercept"].copy()
        table["term"] = table["term"].replace(
            {
                "C(line, Sum)": "line",
                "C(dataset, Sum)": "date",
                "C(line, Sum):C(dataset, Sum)": "line:date",
            }
        )
        table["check"] = "cross_date_pooling"
        table["condition"] = condition
        rows.append(table)

    checks = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    checks.to_excel(LOG_DIR / "initial_gfps_pooling_feasibility_checks.xlsx", index=False)
    return checks


def contrast_subset(contrasts, dataset, assay, condition, batch):
    if contrasts.empty:
        return contrasts
    mask = (
        (contrasts["dataset"] == dataset)
        & (contrasts["assay"] == assay)
        & (contrasts["condition"] == condition)
        & (contrasts["batch"] == batch)
    )
    return contrasts[mask].copy()


def draw_grouped_bars(
    ax,
    data,
    x_col,
    x_order,
    line_order,
    ylabel,
    contrasts,
    title=None,
    dataset=None,
    assay="GFP",
    batch="batch1",
    y_scale=1,
    x_tick_labels=None,
):
    width = min(0.13, 0.56 / max(1, len(line_order)))
    offsets = (np.arange(len(line_order)) - (len(line_order) - 1) / 2) * width * 1.75
    ymax = 0
    bar_positions = {}
    group_tops = {}

    for x_i, condition in enumerate(x_order):
        condition_data = data[data[x_col] == condition]
        family_contrasts = contrast_subset(contrasts, dataset, assay, condition, batch)
        for line_i, line in enumerate(line_order):
            values = condition_data.loc[condition_data["line"] == line, "response"].dropna()
            if values.empty:
                continue
            scaled_values = values / y_scale
            x = x_i + offsets[line_i]
            mean = scaled_values.mean()
            error = sem(scaled_values) if ERRBAR == "SEM" else scaled_values.std(ddof=1)
            ymax = max(ymax, float(scaled_values.max()), float(mean + error))
            group_tops[condition] = max(
                group_tops.get(condition, 0),
                float(scaled_values.max()),
                float(mean + error),
            )
            bar_positions[(condition, line)] = x
            ax.bar(
                x,
                mean,
                width=width,
                color=COLORS.get(line, PALETTE_SEQUENCE[line_i]),
                alpha=BAR_ALPHA,
                edgecolor="black",
                linewidth=BAR_LINEWIDTH,
                zorder=2,
            )
            ax.errorbar(
                x,
                mean,
                yerr=error,
                color=ERRORBAR_COLOR,
                capsize=3,
                linewidth=BAR_LINEWIDTH,
                zorder=3,
            )
            jitter = np.linspace(-width * 0.16, width * 0.16, len(values))
            ax.scatter(
                x + jitter,
                scaled_values,
                s=22,
                marker=MARKERS.get(line, "o"),
                facecolor=COLORS.get(line, PALETTE_SEQUENCE[line_i]),
                edgecolor="none",
                alpha=INDIVIDUAL_POINT_ALPHA,
                zorder=4,
            )
        if family_contrasts is not None and not family_contrasts.empty:
            step = max(ymax, group_tops.get(condition, 1), 1) * 0.09
            bracket_y = group_tops.get(condition, 0) + step * 0.85
            bracket_h = step * 0.48
            for _, row in family_contrasts.iterrows():
                left = bar_positions.get((condition, row["left"]))
                right = bar_positions.get((condition, row["right"]))
                if left is None or right is None:
                    continue
                label = row.get("stars", "") or "ns"
                add_star_bracket(ax, left, right, bracket_y, bracket_h, label)
                ymax = max(ymax, bracket_y + bracket_h + step * 0.3)
                bracket_y += step

    ax.set_xticks(range(len(x_order)))
    ax.set_xticklabels(x_tick_labels or x_order)
    ax.set_ylabel(scaled_ylabel(ylabel, y_scale))
    if title:
        ax.set_title(title, pad=8)
    ax.set_ylim(0, ymax * 1.24 if ymax else 1)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax)


def draw_single_condition_bars(
    ax,
    data,
    line_order,
    ylabel,
    contrasts,
    title,
    dataset,
    assay,
    batch="batch1",
    condition="virB",
    y_scale=1,
):
    width = BAR_WIDTH
    ymax = 0
    family_contrasts = contrast_subset(contrasts, dataset, assay, condition, batch)
    bar_positions = {}
    positions = (np.arange(len(line_order)) - (len(line_order) - 1) / 2) * BAR_CENTER_SPACING
    scaled_means = {}
    counts = {}
    for i, (line, x) in enumerate(zip(line_order, positions)):
        values = data.loc[data["line"] == line, "response"].dropna()
        if values.empty:
            continue
        scaled_values = values / y_scale
        mean = scaled_values.mean()
        error = sem(scaled_values) if ERRBAR == "SEM" else scaled_values.std(ddof=1)
        ymax = max(ymax, float(scaled_values.max()), float(mean + error))
        bar_positions[line] = x
        scaled_means[line] = float(mean)
        counts[line] = len(values)
        ax.bar(
            x,
            mean,
            width=width,
            color=COLORS.get(line, PALETTE_SEQUENCE[i]),
            alpha=BAR_ALPHA,
            edgecolor="black",
            linewidth=BAR_LINEWIDTH,
            zorder=2,
        )
        ax.errorbar(x, mean, yerr=error, color=ERRORBAR_COLOR, capsize=3, linewidth=BAR_LINEWIDTH, zorder=3)
        jitter = np.linspace(-width * 0.14, width * 0.14, len(values))
        ax.scatter(
            x + jitter,
            scaled_values,
            s=22,
            marker=MARKERS.get(line, "o"),
            facecolor=COLORS.get(line, PALETTE_SEQUENCE[i]),
            edgecolor="none",
            alpha=INDIVIDUAL_POINT_ALPHA,
            zorder=4,
        )
    if "WT" in scaled_means:
        prp_positions = [
            bar_positions[line] + width / 2
            for line in line_order
            if line != "WT" and line in bar_positions
        ]
        if prp_positions:
            wt_x = bar_positions["WT"]
            ax.hlines(
                scaled_means["WT"],
                wt_x,
                max(prp_positions),
                colors=WT_REFERENCE_COLOR,
                linestyles=(0, (1.2, 2.0)),
                linewidth=0.8,
                zorder=5,
            )
    if family_contrasts is not None and not family_contrasts.empty:
        family_contrasts = family_contrasts.assign(
            span=family_contrasts.apply(
                lambda row: abs(
                    line_order.index(row["left"]) - line_order.index(row["right"])
                )
                if row["left"] in line_order and row["right"] in line_order
                else 0,
                axis=1,
            )
        ).sort_values("span")
        step = max(ymax, 1) * (0.15 if len(family_contrasts) > 1 else 0.115)
        bracket_y = ymax + step * 0.85
        bracket_h = step * 0.55
        for _, row in family_contrasts.iterrows():
            left = bar_positions.get(row["left"])
            right = bar_positions.get(row["right"])
            if left is None or right is None:
                continue
            label = row.get("stars", "") or "ns"
            add_star_bracket(ax, left, right, bracket_y, bracket_h, label)
            ymax = max(ymax, bracket_y + bracket_h + step * 0.3)
            bracket_y += step
    if counts:
        n_text = n_label_from_counts(counts)
        ax.text(
            0.98,
            0.97,
            n_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=6.8,
            color=TEXT_COLOR,
        )
    ax.set_xticks([bar_positions[line] for line in line_order if line in bar_positions])
    ax.set_xticklabels(
        [display_line_label(line) for line in line_order if line in bar_positions],
        rotation=30,
        ha="right",
    )
    for tick_label in ax.get_xticklabels():
        if tick_label.get_text() in {"prp27-1", "prp27-2"}:
            tick_label.set_fontstyle("italic")
    ax.set_ylabel(scaled_ylabel(ylabel, y_scale))
    ax.set_title(title, pad=8)
    ax.set_ylim(0, ymax * 1.24 if ymax else 1)
    if len(bar_positions) > 0:
        xs = list(bar_positions.values())
        ax.set_xlim(
            min(xs) - BAR_CENTER_SPACING * 0.95,
            max(xs) + BAR_CENTER_SPACING * 0.95,
        )
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax)


def add_legend(
    fig,
    line_order,
    anchor=(0.5, 0.98),
    ncol=None,
    loc="upper center",
    frameon=False,
    title=None,
    fontsize=9,
    title_fontsize=9,
    frame_linewidth=1.5,
    borderpad=0.4,
    labelspacing=0.4,
    columnspacing=1.0,
    handletextpad=0.5,
):
    handles = [
        plt.Line2D(
            [0],
            [0],
            color=COLORS.get(line, PALETTE_SEQUENCE[i]),
            marker=MARKERS.get(line, "o"),
            markerfacecolor=COLORS.get(line, PALETTE_SEQUENCE[i]),
            markeredgecolor="none",
            linewidth=0,
            markersize=5.5,
            label=display_line_label(line),
        )
        for i, line in enumerate(line_order)
    ]
    legend = fig.legend(
        handles=handles,
        loc=loc,
        bbox_to_anchor=anchor,
        ncol=ncol or len(line_order),
        frameon=frameon,
        fancybox=True,
        facecolor="none",
        edgecolor=TEXT_COLOR,
        fontsize=fontsize,
        title=title,
        title_fontsize=title_fontsize,
        borderpad=borderpad,
        labelspacing=labelspacing,
        columnspacing=columnspacing,
        handletextpad=handletextpad,
    )
    for legend_text in legend.get_texts():
        if legend_text.get_text() in {"prp27-1", "prp27-2"}:
            legend_text.set_fontstyle("italic")
    if frameon:
        frame = legend.get_frame()
        frame.set_linewidth(frame_linewidth)
        frame.set_alpha(1)
        frame.set_facecolor("none")
        frame.set_fill(False)
    return legend


def make_figures(cleaned, contrasts):
    endpoint = cleaned[(cleaned["assay"] == "GFP") & (cleaned["condition_type"] == "p19")]

    def plot_gfp_condition(ax, dataset, batch, condition, title, line_order=None):
        if dataset == "20250708":
            leaf_condition = condition if condition in LEAF_ORDER else "Middle"
            subset = cleaned[
                (cleaned["dataset"] == "20250708")
                & (cleaned["assay"] == "GFP")
                & (cleaned["condition"] == leaf_condition)
            ]
            stats_condition = leaf_condition
        else:
            subset = endpoint[
                (endpoint["dataset"] == dataset)
                & (endpoint["batch"] == batch)
                & (endpoint["condition"] == condition)
            ]
            stats_condition = condition
        draw_single_condition_bars(
            ax,
            subset,
            line_order or LINE_ORDERS[dataset],
            "Background-normalised GFP signal",
            contrasts,
            title,
            dataset,
            "GFP",
            batch=batch,
            condition=stats_condition,
            y_scale=GFP_Y_SCALE,
        )

    fig, axes = plt.subplots(5, 2, figsize=(7.4, 15.2), constrained_layout=False)
    comparison_specs = [
        (axes[0, 0], "20250505", "batch1", "no p19", "20250505 no p19"),
        (axes[0, 1], "20250505", "batch1", "p19", "20250505 p19"),
        (axes[1, 0], "20260330", "batch1", "no p19", "20260330 batch1 no p19"),
        (axes[1, 1], "20260330", "batch1", "p19", "20260330 batch1 p19"),
        (axes[2, 0], "20260330", "batch2", "no p19", "20260330 batch2 no p19"),
        (axes[2, 1], "20260330", "batch2", "p19", "20260330 batch2 p19"),
        (axes[3, 0], "20260415", "batch1", "no p19", "20260415 no p19"),
        (axes[3, 1], "20260415", "batch1", "p19", "20260415 p19"),
        (axes[4, 0], "20250708", "batch1", "Middle", "20250708 no p19"),
    ]
    for ax, dataset, batch, condition, title in comparison_specs:
        plot_gfp_condition(ax, dataset, batch, condition, title)
    axes[4, 1].axis("off")
    fig.subplots_adjust(left=0.105, right=0.975, top=0.945, bottom=0.175, wspace=0.34, hspace=1.28)
    add_legend(
        fig,
        ["WT", "PRp27#1", "PRp27#2"],
        anchor=(0.91, 0.255),
        ncol=1,
        loc="upper right",
        frameon=True,
        title="Plant line",
        fontsize=10,
        title_fontsize=10,
        frame_linewidth=0.8,
        borderpad=0.85,
        labelspacing=0.8,
        columnspacing=1.2,
        handletextpad=0.9,
    )
    for i, color in enumerate(PANEL_BAND_COLORS):
        row_axes = axes[i, :] if i < 4 else [axes[4, 0]]
        add_axis_group_band(
            fig,
            row_axes,
            color,
            pad_left=0.09,
            pad_bottom=0.052,
            alpha=0.45 if i == 4 else 1.0,
        )
    save(fig, "gfp_endpoint_all_dates_unpooled")
    plt.close(fig)

    # Alternate presentation of the same unpooled endpoints. The coloured row
    # bands retain the dataset grouping, while the panel headings show only the
    # experimental p19 condition.
    fig, axes = plt.subplots(5, 2, figsize=(7.4, 15.2), constrained_layout=False)
    condition_title_specs = [
        (axes[0, 0], "20250505", "batch1", "no p19", "no p19"),
        (axes[0, 1], "20250505", "batch1", "p19", "p19"),
        (axes[1, 0], "20260330", "batch1", "no p19", "no p19"),
        (axes[1, 1], "20260330", "batch1", "p19", "p19"),
        (axes[2, 0], "20260330", "batch2", "no p19", "no p19"),
        (axes[2, 1], "20260330", "batch2", "p19", "p19"),
        (axes[3, 0], "20260415", "batch1", "no p19", "no p19"),
        (axes[3, 1], "20260415", "batch1", "p19", "p19"),
        (axes[4, 0], "20250708", "batch1", "Middle", "no p19"),
    ]
    for ax, dataset, batch, condition, title in condition_title_specs:
        plot_gfp_condition(ax, dataset, batch, condition, title)
    axes[4, 1].axis("off")
    fig.subplots_adjust(left=0.105, right=0.975, top=0.945, bottom=0.175, wspace=0.34, hspace=1.28)
    add_legend(
        fig,
        ["WT", "PRp27#1", "PRp27#2"],
        anchor=(0.91, 0.255),
        ncol=1,
        loc="upper right",
        frameon=True,
        title="Plant line",
        fontsize=10,
        title_fontsize=10,
        frame_linewidth=0.8,
        borderpad=0.85,
        labelspacing=0.8,
        columnspacing=1.2,
        handletextpad=0.9,
    )
    for i, color in enumerate(PANEL_BAND_COLORS):
        row_axes = axes[i, :] if i < 4 else [axes[4, 0]]
        add_axis_group_band(
            fig,
            row_axes,
            color,
            pad_left=0.09,
            pad_bottom=0.052,
            alpha=0.45 if i == 4 else 1.0,
        )
    save(fig, "gfp_endpoint_all_dates_condition_titles")
    plt.close(fig)

    data_4wk = cleaned[
        (cleaned["dataset"] == "202405_4wk") & (cleaned["assay"] == "GFP")
    ]
    fig, ax = plt.subplots(figsize=(3.8, 4.0), constrained_layout=False)
    draw_single_condition_bars(
        ax,
        data_4wk,
        LINE_ORDERS["202405_4wk"],
        "Background-normalised GFP signal",
        contrasts,
        "",
        "202405_4wk",
        "GFP",
        condition="4-week",
        y_scale=GFP_Y_SCALE,
    )
    add_legend(fig, LINE_ORDERS["202405_4wk"], anchor=(0.5, 0.93), ncol=2)
    fig.subplots_adjust(left=0.18, right=0.97, top=0.74, bottom=0.32)
    save(fig, "gfp_202405_4wk_endpoint")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(5.9, 3.8), constrained_layout=False)
    plot_gfp_condition(axes[0], "20250505", "batch1", "no p19", "20250505 no p19")
    plot_gfp_condition(axes[1], "20250505", "batch1", "p19", "20250505 p19")
    add_legend(fig, LINE_ORDERS["20250505"], anchor=(0.5, 0.985), ncol=2)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.78, bottom=0.31, wspace=0.36)
    add_axis_group_band(fig, axes, PANEL_BAND_COLORS[0])
    save(fig, "gfp_20250505_endpoint")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.8), constrained_layout=False)
    for ax, batch, condition, title in [
        (axes[0, 0], "batch1", "no p19", "batch1 no p19"),
        (axes[0, 1], "batch1", "p19", "batch1 p19"),
        (axes[1, 0], "batch2", "no p19", "batch2 no p19"),
        (axes[1, 1], "batch2", "p19", "batch2 p19"),
    ]:
        plot_gfp_condition(ax, "20260330", batch, condition, f"20260330 {title}")
    add_legend(fig, LINE_ORDERS["20260330"], anchor=(0.5, 0.99), ncol=2)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.85, bottom=0.23, wspace=0.34, hspace=0.62)
    add_axis_group_band(fig, axes[0, :], PANEL_BAND_COLORS[1])
    add_axis_group_band(fig, axes[1, :], PANEL_BAND_COLORS[2])
    save(fig, "gfp_20260330_endpoint_by_batch")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.6), constrained_layout=False)
    plot_gfp_condition(axes[0], "20260415", "batch1", "no p19", "20260415 no p19")
    plot_gfp_condition(axes[1], "20260415", "batch1", "p19", "20260415 p19")
    add_legend(fig, LINE_ORDERS["20260415"], anchor=(0.5, 0.985), ncol=3)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.78, bottom=0.31, wspace=0.38)
    add_axis_group_band(fig, axes, PANEL_BAND_COLORS[3])
    save(fig, "gfp_20260415_endpoint")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.4, 3.8), constrained_layout=False)
    plot_gfp_condition(ax, "20250708", "batch1", "Middle", "20250708 middle")
    add_legend(fig, LINE_ORDERS["20250708"], anchor=(0.5, 0.985), ncol=2)
    fig.subplots_adjust(left=0.18, right=0.97, top=0.80, bottom=0.31)
    add_axis_group_band(fig, [ax], PANEL_BAND_COLORS[4])
    save(fig, "gfp_20250708_middle_endpoint")
    plt.close(fig)

    # Direct age comparison using the final 4-week endpoint and the middle-leaf
    # endpoint from the 20250708 6-week experiment.
    fig, axes = plt.subplots(1, 2, figsize=(6.2, 3.8), constrained_layout=False)
    draw_single_condition_bars(
        axes[0],
        data_4wk,
        LINE_ORDERS["202405_4wk"],
        "Background-normalised GFP signal",
        contrasts,
        "4-week",
        "202405_4wk",
        "GFP",
        condition="4-week",
        y_scale=GFP_Y_SCALE,
    )
    plot_gfp_condition(axes[1], "20250708", "batch1", "Middle", "6-week")
    for ax in axes:
        ax.set_facecolor("#E9EFF6")
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="white", linewidth=1.2)
        ax.set_ylabel("GFP signal intensity (x1000)")
        style_age_comparison(ax)
    fig.subplots_adjust(left=0.115, right=0.975, top=0.89, bottom=0.24, wspace=0.42)
    save(fig, "gfp_4wk_6wk_endpoint")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.8), constrained_layout=False)
    for ax, leaf in zip(axes, LEAF_ORDER):
        plot_gfp_condition(ax, "20250708", "batch1", leaf, f"20250708 {leaf.lower()}")
    add_legend(fig, LINE_ORDERS["20250708"], anchor=(0.5, 0.985), ncol=2)
    fig.subplots_adjust(left=0.105, right=0.98, top=0.78, bottom=0.31, wspace=0.38)
    add_axis_group_band(fig, axes, PANEL_BAND_COLORS[4])
    save(fig, "gfp_20250708_leaf_positions_endpoint")
    plt.close(fig)

    virb = cleaned[(cleaned["dataset"] == "20250708") & (cleaned["assay"] == "virB")]
    fig, ax = plt.subplots(figsize=(3.4, 3.8), constrained_layout=False)
    draw_single_condition_bars(
        ax,
        virb,
        LINE_ORDERS["20250708"],
        "Background-normalised virB signal",
        contrasts,
        "20250708 virB",
        "20250708",
        "virB",
        y_scale=VIRB_Y_SCALE,
    )
    fig.subplots_adjust(left=0.21, right=0.97, top=0.84, bottom=0.32)
    add_axis_group_band(fig, [ax], PANEL_BAND_COLORS[4])
    save(fig, "virb_20250708_endpoint")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.9), constrained_layout=False)
    plot_gfp_condition(axes[0], "20250708", "batch1", "Middle", "GFP middle")
    draw_single_condition_bars(
        axes[1],
        virb,
        LINE_ORDERS["20250708"],
        "Background-normalised virB signal",
        contrasts,
        "virB",
        "20250708",
        "virB",
        y_scale=VIRB_Y_SCALE,
    )
    add_legend(fig, LINE_ORDERS["20250708"], anchor=(0.5, 0.99), ncol=2)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.79, bottom=0.30, wspace=0.38)
    add_axis_group_band(fig, axes, PANEL_BAND_COLORS[4])
    save(fig, "gfp_virb_20250708_middle_endpoint")
    plt.close(fig)


def make_cleaning_reference(screened, cleaned):
    before = screened.assign(cleaning_stage="pre-clean")
    after = cleaned.assign(cleaning_stage="post-clean")
    plot_data = pd.concat([before, after], ignore_index=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.8), constrained_layout=False)
    for ax, assay, ylabel in zip(
        axes,
        ["GFP", "virB"],
        ["Background-normalised GFP signal", "Background-normalised virB signal"],
    ):
        subset = plot_data[plot_data["assay"] == assay]
        positions = {"pre-clean": 0, "post-clean": 1}
        for i, stage in enumerate(["pre-clean", "post-clean"]):
            values = subset.loc[subset["cleaning_stage"] == stage, "response"].dropna()
            mean = values.mean()
            error = sem(values)
            ax.bar(
                i,
                mean,
                width=0.58,
                color=PALETTE_SEQUENCE[i],
                alpha=BAR_ALPHA,
                edgecolor="black",
                linewidth=BAR_LINEWIDTH,
                zorder=2,
            )
            ax.errorbar(i, mean, yerr=error, color=ERRORBAR_COLOR, capsize=3, linewidth=BAR_LINEWIDTH, zorder=3)
            jitter = np.linspace(-0.12, 0.12, len(values))
            ax.scatter(i + jitter, values, s=18, color=PALETTE_SEQUENCE[i], alpha=0.35, edgecolor="none", zorder=4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(list(positions))
        ax.set_ylabel(axis_ylabel(ylabel))
        ax.set_title(assay)
        style_axis(ax)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.86, bottom=0.28, wspace=0.35)
    save(fig, "cleaning_reference_pre_vs_post")
    plt.close(fig)


# %%
clear_previous_outputs()
raw = load_all()
analysis_input, excluded_rows = make_analysis_input(raw)
cleaned, screened, review_log = clean_data(analysis_input)

raw.to_excel(LOG_DIR / "initial_gfps_raw_long.xlsx", index=False)
analysis_input.to_excel(LOG_DIR / "initial_gfps_analysis_input_filtered.xlsx", index=False)
excluded_rows.to_excel(LOG_DIR / "initial_gfps_excluded_sample_types.xlsx", index=False)
screened.to_excel(LOG_DIR / "initial_gfps_outlier_screened_all_rows.xlsx", index=False)
review_log.to_excel(LOG_DIR / "initial_gfps_outlier_log.xlsx", index=False)
cleaned.to_excel(LOG_DIR / "initial_gfps_cleaned_long.xlsx", index=False)

anova, contrasts = run_stats(cleaned)
pooling_checks = run_combination_checks(cleaned)
make_figures(cleaned, contrasts)
make_cleaning_reference(screened, cleaned)

print(f"Loaded {len(raw)} rows from {len(SOURCE_FILES)} source workbooks.")
print(f"Excluded {len(excluded_rows)} non-target rows before analysis: {', '.join(EXCLUDED_LINES)}.")
print(f"Flagged {int(screened['flag_outlier'].sum())} potential outlier rows; removed {int(screened['manual_exclude'].sum())}.")
if not pooling_checks.empty:
    print(f"Pooling feasibility checks written with {len(pooling_checks)} rows.")
print(f"Figures saved to: {SAVE_DIR}")
print(f"Logs saved to: {LOG_DIR}")
