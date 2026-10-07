# %% [markdown]
# # VIGS screen batch GFP analysis
#
# First-pass workflow for per-batch VIGS screen plots. The script inventories
# batch folders, extracts usable Excel GFP tables, compiles batch metadata, runs
# conservative within-batch control comparisons, and saves per-batch bar plots.

# %%
from pathlib import Path
import os
import re
import sys
import zipfile

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

from PIL import Image, ImageOps
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
else:
    HERE = Path.cwd() / "VIGS" / "VIGS Screen"

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))

from figure_styles import (  # noqa: E402
    ANALYSIS_RCPARAMS,
    ERRORBAR_COLOR,
    FONT_STACK,
    INDIVIDUAL_POINT_ALPHA,
    TEXT_COLOR,
    apply_axis_style,
)

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_vigs_screen"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
IMAGE_DIR = HERE / "Images"
PANEL_IMAGE_DIR = HERE / "BatchPhotos"
EXTRACTED_IMAGE_DIR = IMAGE_DIR / "extracted_pptx_media"
for folder in [SAVE_DIR, LOG_DIR, IMAGE_DIR, PANEL_IMAGE_DIR, EXTRACTED_IMAGE_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

SAVE_FMTS = ["png", "pdf", "svg"]
SOURCE_ROOT = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/4. Immunity Silencing/VIGS/Batches"
)
VIGS_LIBRARY = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/4. Immunity Silencing/VIGS/VIGSLibrary.xlsx"
)
CONSOLIDATED_BATCHES = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/4. Immunity Silencing/VIGS/Consolidated_VIGS_Batches.xlsx"
)

# The standalone workbook filed under Batch4 contains the Batch5 measurements
# (plus JAZ). The genuine 20240616 Batch4 table is retained in the consolidated
# workbook, so use that named sheet explicitly and preserve the override in the
# provenance columns/logs.
BATCH_SOURCE_OVERRIDES = {
    "Batch3": {
        "path": SOURCE_ROOT / "Batch3" / "VIGSdata_20240527_and0604.xlsx",
        "sheet": "Raw",
    },
    "Batch4": {
        "path": SOURCE_ROOT / "Batch4" / "20240616" / "VIGSdata_20240616.xlsx",
        "sheet": "All",
    },
    "Batch9": {
        "path": SOURCE_ROOT / "Batch9" / "20241204" / "20241204and17.xlsx",
        "sheet": "Sheet1",
    },
}

CONTROL_PRIORITY = ["GUS"]
BASE_COLOR = "#3E6E86"
CONTROL_COLOR = "#6D6D6D"
BAR_ALPHA = 0.68
BAR_LINEWIDTH = 1.2
SIG_ALPHA = 0.05
MAX_CONTACT_IMAGES = 6

LINE_ALIASES = {
    "ACIF": "ACIF1",
    "ACIF1": "ACIF1",
    "CRT3": "CRT3A",
    "CRT3A": "CRT3A",
    "MAPK": "MPK3/6",
    "MPK": "MPK3/6",
    "MPK3/6": "MPK3/6",
    "CPDK": "CDPK",
    "PEROX": "PEROX",
    "PEROX7": "PEROX",
    "PEROX72": "PEROX",
    "PEROX72-LIKE": "PEROX",
    "PEROXIDASE": "PEROX",
    "PEROXIDASE 72-LIKE": "PEROX",
    "PEROXIDASE 72-LIKE/49-LIKE": "PEROX",
    "PR17": "PRP27",
    "PRP27": "PRP27",
    "RBOHA": "RBOHA",
    "RBOHB": "RBOHB",
    "HEVA": "HEVA",
    "HEVAMINE": "HEVA",
}

TARGET_DISPLAY_ALIASES = {
    "ACIF1": "ACIF1",
    "CRT3A": "CRT3a",
    "MPK3/6": "MPK3/6",
    "PEROX": "PEROX72",
    "PRP27": "PRp27",
    "RBOHA": "RbohA",
    "RBOHB": "RbohB",
    "HEVA": "HevA",
    "NBD030152": "NbD030152",
}

TARGETS_EXCLUDED_FROM_SCREEN = {
    "GOX4": "construct issues",
    "PYL1": "construct issues",
    "STP13": "construct/source issue in target list",
    "LYK5": "universal exclusion",
    "LYM1": "excluded in target list",
    "NBD030152": "universal exclusion",
    "BIK1": "universal exclusion",
}

BIOLOGICAL_UNIVERSAL_EXCLUSIONS = {
    "CDPK", "EAS", "EDS5", "ETR1", "HEVA", "JAZ",
    "PYL8", "RBOHA", "RBOHB", "SOBIR1",
}

PLOT_EXCLUDED_LINES = set(TARGETS_EXCLUDED_FROM_SCREEN) | BIOLOGICAL_UNIVERSAL_EXCLUSIONS

BATCH_MANIFEST = {
    "Batch1": {
        "date": "20250529",
        "included": [
            "ACIF1",
            "CERK1",
            "COI1",
            "CORE",
            "CRT3A",
            "CSPR",
            "EAH",
            "EAS",
            "EDS1",
            "EDS5",
            "EIN2",
            "EIN3",
            "GUS",
            "ICS",
            "MPK3/6",
            "NOA1",
            "NPR1",
            "NRC234",
            "PEROX",
            "PR3",
            "PRP27",
            "RBOHB",
            "SAG101",
            "WAK1",
            "WRKY",
        ],
        "excluded": ["BAK1", "ETR1", "GOX4", "PYL1", "RBOHA", "PYL8", "STP13"],
    },
    "Batch2": {
        "date": "20240318",
        "included": ["GUS", "WAK1", "PR3", "EIN3", "EIN2", "CSPR"],
        "excluded": ["PYL1"],
    },
    "Batch3": {
        "date": "20240527_and0604",
        "included": ["GUS", "ACIF1", "CORE", "CRT3A", "EDS1", "PAL", "PRP27", "SAG101", "WRKY"],
        "excluded": ["STP13", "ETR1"],
    },
    "Batch4": {
        "date": "20240616",
        "included": ["PR3", "NOA1", "GUS", "EIN2"],
        "excluded": ["GOX4"],
    },
    "Batch5": {
        "date": "20240629",
        "included": ["GUS", "MPK3/6", "COI1", "LYM1", "LYM3", "RBOHA"],
        "excluded": ["JAZ", "PYL8", "SOBIR1"],
    },
    "Batch6": {"date": "20240730", "included": ["GUS", "CERK1"], "excluded": []},
    "Batch7": {
        "date": "20240909",
        "included": ["GUS", "CNGC2"],
        "excluded": ["CDPK", "LIK1", "BIK1"],
    },
    "Batch8": {
        "date": "20241008",
        "included": ["GUS", "NPR1", "LYK3"],
        "excluded": ["EDS5", "LYK5", "HevA"],
    },
    "Batch9": {
        "date": "20241115",
        "included": ["GUS", "LYK3", "NRC234", "EAH", "PEROX"],
        "excluded": [],
    },
    "Batch10": {
        "date": "20250623",
        "included": ["GUS", "CSPR", "ICS", "EDS1", "BAK1", "EAH", "CERK1", "MPK3/6"],
        "excluded": ["EIN3", "EAS", "COI1", "JAZ", "CORE", "EDS5", "ACIF1", "GOX4", "CRT3A"],
    },
    "Batch11": {"date": "20251031", "included": ["GUS"], "excluded": ["BIK1", "SOBIR1"]},
    "Batch12": {"date": "20251229", "included": ["GUS"], "excluded": ["ICS", "EAS", "ETR1"]},
    "Batch13": {"date": "20260105", "included": ["GUS"], "excluded": ["PYL1", "PYL8"]},
    "Batch14": {
        "date": "20260211",
        "included": ["GUS", "PYL8", "NRC234", "WAK1", "SAG101"],
        "excluded": ["PRP27", "BIK1", "PYL1", "RBOHB", "HevA", "LYK3", "WRKY", "RBOHA", "NPR1", "SOBIR1", "CDPK"],
    },
    "Batch15": {
        "date": "20260212",
        "included": [
            "GUS", "NRC234", "SQS", "PEROX", "LYM3", "MPK3/6",
        ],
        "excluded": ["RBOHB", "JAZ", "PAL", "CNGC2", "NOA1", "PYL8", "LYM1"],
    },
    "Batch16": {
        "date": "20260530",
        "included": ["GUS", "PEROX", "LYK3", "LIK1", "EAS", "NRC234", "EAH"],
        "excluded": ["LYM3", "JAZ", "SQS"],
    },
}

EXCESS_BATCH_EXCLUSIONS = [
    ("Batch1", "EAH"),
    ("Batch1", "EAS"),
    ("Batch1", "ETR1"),
    ("Batch8", "LYK3"),
    ("Batch1", "MPK3/6"),
    ("Batch9", "NRC234"),
    ("Batch1", "NRC234"),
    ("Batch1", "PEROX"),
    ("Batch1", "PYL1"),
    ("Batch5", "PYL8"),
    ("Batch1", "PYL8"),
    ("Batch1", "SQS"),
    ("Batch5", "JAZ"),
]

MANUAL_BATCH_EXCLUSIONS = [
    {
        "batch": "Batch3",
        "plant_line": "PAL",
        "exclusion_source": "quality exclusion: second-date set lacks same-day GUS control",
    },
] + [
    {
        "batch": batch,
        "plant_line": target,
        "exclusion_source": "quality exclusion: excess batch; only three highest-numbered batch attempts retained",
    }
    for batch, target in EXCESS_BATCH_EXCLUSIONS
]


VIGS_RCPARAMS = {
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 11.5,
    "axes.titlesize": 12.5,
    "axes.labelsize": 11.5,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 10.5,
}
plt.rcParams.update(VIGS_RCPARAMS)


# %%
def sem(values):
    values = pd.Series(values).dropna()
    if len(values) <= 1:
        return 0.0
    return values.std(ddof=1) / np.sqrt(len(values))


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


def clean_label(value):
    if pd.isna(value):
        return ""
    label = str(value).strip()
    label = re.sub(r"^\d+\.", "", label)
    return label.strip()


def canonical_line(value):
    label = clean_label(value)
    upper = label.upper()
    if upper in {"1.GUS", "GUS", "GUS "}:
        return "GUS"
    if upper == "WT":
        return "GUS"
    upper = upper.replace("PRP27", "PRP27")
    return LINE_ALIASES.get(upper, upper)


def display_line(value):
    line = canonical_line(value)
    return TARGET_DISPLAY_ALIASES.get(line, line)


def sort_key_label(label):
    upper = str(label).upper()
    if upper == "GUS":
        return (0, upper)
    return (1, upper)


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3, direction="out")
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{output_stem}.{file_format}", dpi=300, bbox_inches="tight")




# %%
def infer_dates_from_path(path):
    return sorted(set(re.findall(r"20\d{6}", str(path))))


def is_second_part_file(path):
    text = str(path).lower()
    return any(token in text for token in ["part2", "part_2", "and17", "0604"])


def batch_sort_key(path):
    match = re.search(r"Batch(\d+)", path.name, flags=re.IGNORECASE)
    return int(match.group(1)) if match else 999


def inventory_batches():
    rows = []
    flags = []
    for batch_dir in sorted(
        [p for p in SOURCE_ROOT.iterdir() if p.is_dir() and p.name.lower() != "scrap"],
        key=batch_sort_key,
    ):
        excel_files = sorted(batch_dir.rglob("*.xlsx"))
        pptx_files = sorted(batch_dir.rglob("*.pptx"))
        prism_files = sorted(batch_dir.rglob("*.prism"))
        dates = sorted(set(date for f in excel_files + pptx_files + prism_files for date in infer_dates_from_path(f)))
        data_excels = [f for f in excel_files if not is_second_part_file(f)]
        second_part_excels = [f for f in excel_files if is_second_part_file(f)]

        primary_excel = None
        status = "ready"
        reason = ""
        if not data_excels:
            status = "flagged"
            reason = "No non-second-part Excel file found."
        elif len(data_excels) > 1:
            status = "flagged"
            reason = "Multiple candidate primary Excel files found."
        else:
            primary_excel = data_excels[0]
        if second_part_excels:
            flags.append(
                {
                    "batch": batch_dir.name,
                    "flag_type": "ignored_second_part_excel",
                    "details": "; ".join(str(p) for p in second_part_excels),
                }
            )
        if any("and" in f.name.lower() and len(infer_dates_from_path(f)) >= 2 for f in excel_files):
            flags.append(
                {
                    "batch": batch_dir.name,
                    "flag_type": "combined_date_filename",
                    "details": "; ".join(str(f) for f in excel_files if "and" in f.name.lower()),
                }
            )
            status = "flagged"
            reason = "Excel filename suggests combined date file; needs manual split decision."

        rows.append(
            {
                "batch": batch_dir.name,
                "batch_path": str(batch_dir),
                "status": status,
                "flag_reason": reason,
                "primary_excel": str(primary_excel) if primary_excel else "",
                "excel_files": "; ".join(str(p) for p in excel_files),
                "pptx_files": "; ".join(str(p) for p in pptx_files),
                "prism_files": "; ".join(str(p) for p in prism_files),
                "dates_detected": "; ".join(dates),
                "n_excel": len(excel_files),
                "n_pptx": len(pptx_files),
                "n_prism": len(prism_files),
            }
        )
        if status == "flagged":
            flags.append({"batch": batch_dir.name, "flag_type": "batch_not_processed", "details": reason})
    return pd.DataFrame(rows), pd.DataFrame(flags)


# %%
TARGET_COLUMN_CANDIDATES = [
    "sampletype",
    "sample type",
    "plant type",
    "plant line",
    "vigs target",
    "vigs target ",
    "vigs target",
    "planttype",
]

VALUE_COLUMN_PREFERENCE = [
    "gus-relative",
    "(signal-background)/gus",
    "normalised gfp",
    "normalized gfp",
    "normalised signal",
    "normalized signal",
    "signal - background",
    "mean gfp intensity background",
    "gfpsignal",
    "gfp signal",
]


def normalized_columns(df):
    return {str(column).strip().lower(): column for column in df.columns}


def read_candidate_sheet(path, forced_sheet=None):
    xl = pd.ExcelFile(path)
    if forced_sheet is not None:
        if forced_sheet not in xl.sheet_names:
            raise ValueError(f"Requested sheet {forced_sheet!r} is absent from {path}.")
        return forced_sheet, read_sheet_with_detected_header(path, forced_sheet)
    preferred = ["raw", "sheet1", "batch10", "all", "processed", "pivot"]
    for wanted in preferred:
        for sheet in xl.sheet_names:
            if sheet.strip().lower() == wanted:
                df = read_sheet_with_detected_header(path, sheet)
                if len(df) > 0:
                    return sheet, df
    sheet = xl.sheet_names[0]
    return sheet, read_sheet_with_detected_header(path, sheet)


def read_sheet_with_detected_header(path, sheet):
    preview = pd.read_excel(path, sheet_name=sheet, header=None, nrows=8)
    header_tokens = {
        "sampletype",
        "sample type",
        "plant type",
        "plant line",
        "vigs target",
        "replicate",
        "background",
        "rawgfpsignal",
        "gfpsignal",
        "normalised gfp",
        "normalized gfp",
        "normalised signal",
    }
    best_row = 0
    best_score = -1
    for idx, row in preview.iterrows():
        values = {str(value).strip().lower() for value in row.dropna().tolist()}
        score = len(values & header_tokens)
        if score > best_score:
            best_row = int(idx)
            best_score = score
    if best_score <= 0:
        best_row = 0
    return pd.read_excel(path, sheet_name=sheet, header=best_row)


def melt_wide_table(df):
    numeric_columns = []
    for column in df.columns:
        values = pd.to_numeric(df[column], errors="coerce")
        if values.notna().sum() >= 2:
            numeric_columns.append(column)
    if len(numeric_columns) < 2:
        return pd.DataFrame()
    long = df[numeric_columns].melt(var_name="plant_line", value_name="normalized_gfp")
    long["plant_line"] = long["plant_line"].map(canonical_line)
    long["normalized_gfp"] = pd.to_numeric(long["normalized_gfp"], errors="coerce")
    long = long.dropna(subset=["normalized_gfp"])
    long["replicate"] = long.groupby("plant_line").cumcount() + 1
    long["source_value_column"] = "wide_numeric_columns"
    return long[["plant_line", "replicate", "normalized_gfp", "source_value_column"]]


def parse_excel_table(path, forced_sheet=None):
    sheet, df = read_candidate_sheet(path, forced_sheet=forced_sheet)
    df = df.dropna(how="all").copy()
    columns = normalized_columns(df)
    target_col = next((columns[name] for name in TARGET_COLUMN_CANDIDATES if name in columns), None)

    if target_col is None:
        long = melt_wide_table(df)
        if long.empty:
            raise ValueError("No recognizable target column or wide numeric sample columns.")
        long["source_sheet"] = sheet
        long["parse_mode"] = "wide_numeric"
        return long

    target = df[target_col].map(canonical_line)
    replicate_col = next((columns[name] for name in ["replicate", "rep"] if name in columns), None)
    background_col = next((columns[name] for name in ["background", "mean background"] if name in columns), None)
    raw_col = next((columns[name] for name in ["gfp", "gfp signal", "rawgfpsignal", "signal"] if name in columns), None)
    value_col = next((columns[name] for name in VALUE_COLUMN_PREFERENCE if name in columns), None)

    source_value_column = value_col
    if value_col is None and raw_col is not None and background_col is not None:
        values = pd.to_numeric(df[raw_col], errors="coerce") - pd.to_numeric(df[background_col], errors="coerce")
        source_value_column = f"{raw_col} - {background_col}"
    elif value_col is not None:
        values = pd.to_numeric(df[value_col], errors="coerce")
    else:
        raise ValueError("No recognizable GFP value column.")

    long = pd.DataFrame(
        {
            "plant_line": target,
            "replicate": pd.to_numeric(df[replicate_col], errors="coerce") if replicate_col else np.nan,
            "normalized_gfp": values,
            "source_value_column": str(source_value_column),
        }
    )
    long = long.dropna(subset=["plant_line", "normalized_gfp"])
    long = long[long["plant_line"].astype(str).str.len() > 0].copy()
    if long["replicate"].isna().all():
        long["replicate"] = long.groupby("plant_line").cumcount() + 1
    long["source_sheet"] = sheet
    long["parse_mode"] = "long_table"
    return long


def choose_control(lines):
    upper_to_label = {str(line).upper(): line for line in lines}
    for control in CONTROL_PRIORITY:
        if control in upper_to_label:
            return upper_to_label[control]
    return None


def add_relative_values(data):
    data = data.copy()
    control = choose_control(data["plant_line"].unique())
    if control is None:
        data["control_line"] = ""
        data["relative_to_control"] = np.nan
        return data, None
    control_mean = data.loc[data["plant_line"] == control, "normalized_gfp"].mean()
    data["control_line"] = control
    data["relative_to_control"] = data["normalized_gfp"] / control_mean if control_mean else np.nan
    return data, control


# %%
def run_batch_stats(data, batch, control):
    values_by_line = [
        group["normalized_gfp"].dropna()
        for _, group in data.groupby("plant_line", sort=False)
        if len(group["normalized_gfp"].dropna()) >= 2
    ]
    if len(values_by_line) >= 2:
        omnibus = stats.kruskal(*values_by_line)
        omnibus_row = {
            "batch": batch,
            "test": "Kruskal-Wallis",
            "statistic": float(omnibus.statistic),
            "p_value": float(omnibus.pvalue),
        }
    else:
        omnibus_row = {"batch": batch, "test": "Kruskal-Wallis", "statistic": np.nan, "p_value": np.nan}

    rows = []
    p_values = []
    control_values = data.loc[data["plant_line"] == control, "normalized_gfp"].dropna()
    if control and len(control_values) >= 2:
        for line, subset in data.groupby("plant_line", sort=False):
            if line == control:
                continue
            values = subset["normalized_gfp"].dropna()
            if len(values) < 2:
                continue
            test = stats.mannwhitneyu(values, control_values, alternative="two-sided")
            rows.append(
                {
                    "batch": batch,
                    "control_line": control,
                    "plant_line": line,
                    "line_n": len(values),
                    "control_n": len(control_values),
                    "line_mean": values.mean(),
                    "control_mean": control_values.mean(),
                    "estimate_line_minus_control": values.mean() - control_values.mean(),
                    "test": "Mann-Whitney U vs control",
                    "statistic": float(test.statistic),
                    "p_raw": float(test.pvalue),
                }
            )
            p_values.append(float(test.pvalue))
    pairwise = pd.DataFrame(rows)
    if not pairwise.empty:
        _, adjusted, _, _ = multipletests(p_values, method="holm")
        pairwise["p_holm_within_batch"] = adjusted
        pairwise["stars"] = [p_to_stars(p) for p in adjusted]
    return pd.DataFrame([omnibus_row]), pairwise


def n_label(counts):
    values = [value for value in counts.values() if value > 0]
    if not values:
        return ""
    low, high = min(values), max(values)
    return f"n={low}" if low == high else f"n={low}-{high}"


def y_scale_for(values):
    max_value = pd.Series(values).dropna().max()
    if pd.isna(max_value) or max_value < 10:
        return 1
    if max_value >= 1000:
        return 1000
    return 1


def plot_batch_bar(ax, data, batch, control, pairwise, master_style=False):
    means = data.groupby("plant_line")["relative_to_control"].mean()
    ordered_lines = [line for line in means.sort_values(ascending=True).index.tolist() if line != control]
    order = ([control] if control in means.index else []) + ordered_lines
    positions = np.arange(len(order), dtype=float)
    y_scale = y_scale_for(data["relative_to_control"])
    ymax = 0
    counts = {}
    colors = [CONTROL_COLOR if line == control else BASE_COLOR for line in order]
    for i, line in enumerate(order):
        values = data.loc[data["plant_line"] == line, "relative_to_control"].dropna() / y_scale
        if values.empty:
            continue
        mean = values.mean()
        err = sem(values)
        ymax = max(ymax, float(values.max()), float(mean + err))
        counts[line] = len(values)
        ax.bar(
            positions[i],
            mean,
            width=0.68,
            color=colors[i],
            alpha=BAR_ALPHA,
            edgecolor="black",
            linewidth=BAR_LINEWIDTH,
            zorder=2,
        )
        ax.errorbar(
            positions[i],
            mean,
            yerr=err,
            fmt="none",
            ecolor=ERRORBAR_COLOR,
            elinewidth=1.1,
            capsize=3.0,
            capthick=1.1,
            zorder=5,
        )
        jitter = np.linspace(-0.11, 0.11, len(values))
        if master_style:
            # Set fill and edge alpha independently: the fill remains translucent
            # while the thinner black outline is completely opaque.
            ax.scatter(
                positions[i] + jitter, values, s=220, marker="o",
                facecolor=to_rgba(colors[i], 0.72),
                edgecolor=(0.0, 0.0, 0.0, 1.0), linewidth=2.5,
                alpha=None, zorder=8,
            )
        else:
            ax.scatter(
                positions[i] + jitter, values, s=24, marker="o",
                facecolor=colors[i], edgecolor="none", linewidth=0,
                alpha=INDIVIDUAL_POINT_ALPHA, zorder=4,
            )

    if control in order:
        ax.axhline(
            1 / y_scale,
            color="#C43B32",
            linestyle=(0, (6, 4)) if master_style else (0, (1.5, 2.4)),
            linewidth=3.2 if master_style else 1.0,
            alpha=1.0,
            zorder=3 if master_style else 1,
        )
    if not pairwise.empty:
        star_map = dict(zip(pairwise["plant_line"], pairwise["stars"]))
        for i, line in enumerate(order):
            if line == control or line not in star_map:
                continue
            if star_map[line] == "ns":
                continue
            values = data.loc[data["plant_line"] == line, "relative_to_control"].dropna() / y_scale
            if values.empty:
                continue
            ax.text(
                positions[i],
                values.max() + max(ymax * 0.08, 0.08),
                star_map[line],
                ha="center",
                va="bottom",
                fontsize=45.0 if master_style else 11.5,
                color=TEXT_COLOR,
            )
            ymax = max(ymax, float(values.max() + max(ymax * 0.15, 0.12)))

    text = n_label(counts)
    if text:
        ax.text(
            0.99, 1.075 if master_style else 0.98,
            text, transform=ax.transAxes, ha="right",
            va="bottom" if master_style else "top",
            fontsize=35.0 if master_style else 9.0, color=TEXT_COLOR,
            clip_on=False,
        )
    ax.set_xticks(positions)
    label_rotation = 60 if len(order) > 12 else 45
    ax.set_xticklabels([display_line(line) for line in order], rotation=label_rotation, ha="right")
    ax.set_ylabel("" if master_style else "GFP expression\n(relative to control)")
    ax.set_ylim(0, ymax * 1.25 if ymax else 1)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax)
    if master_style:
        ax.tick_params(axis="both", labelsize=36.5)


# %%
def extract_pptx_media(pptx_path, batch):
    output_dir = EXTRACTED_IMAGE_DIR / safe_name(batch) / safe_name(pptx_path.stem)
    output_dir.mkdir(parents=True, exist_ok=True)
    extracted = []
    with zipfile.ZipFile(pptx_path) as archive:
        for member in archive.namelist():
            if not member.startswith("ppt/media/"):
                continue
            suffix = Path(member).suffix.lower()
            if suffix not in {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}:
                continue
            output = output_dir / Path(member).name
            if not output.exists():
                output.write_bytes(archive.read(member))
            extracted.append(output)
    return extracted


def readable_image(path):
    try:
        image = ImageOps.exif_transpose(Image.open(path)).convert("RGB")
    except Exception:
        return None
    return image


def make_contact_sheet(image_paths, batch):
    images = [readable_image(path) for path in image_paths]
    images = [image for image in images if image is not None]
    if not images:
        return None
    images = images[:MAX_CONTACT_IMAGES]
    thumb_w, thumb_h = 300, 220
    thumbs = []
    for image in images:
        image.thumbnail((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (thumb_w, thumb_h), "white")
        canvas.paste(image, ((thumb_w - image.width) // 2, (thumb_h - image.height) // 2))
        thumbs.append(canvas)
    sheet = Image.new("RGB", (thumb_w * len(thumbs), thumb_h), "white")
    for i, thumb in enumerate(thumbs):
        sheet.paste(thumb, (i * thumb_w, 0))
    output = IMAGE_DIR / f"{SCRIPT_STEM}_{safe_name(batch)}_contact_sheet.png"
    sheet.save(output)
    return output


def plot_batch_figure(batch_data, batch_meta, pairwise):
    batch = batch_meta["batch"]
    control = batch_data["control_line"].dropna().iloc[0] if batch_data["control_line"].notna().any() else None
    n_lines = batch_data["plant_line"].nunique()
    fig_width = max(7.2, min(14.5, 3.8 + n_lines * 0.34))
    fig, ax = plt.subplots(figsize=(fig_width, 4.9), constrained_layout=False)
    plot_batch_bar(ax, batch_data, batch, control, pairwise)
    fig.subplots_adjust(left=0.14, right=0.97, top=0.94, bottom=0.34)
    save(fig, f"{safe_name(batch)}_relative_gfp")
    plt.close(fig)
    return {
        "batch": batch,
        "pptx_media_extracted": 0,
        "panel_images": "",
    }


# %%
def plant_line_batch_counts(data):
    if data.empty:
        return pd.DataFrame(columns=["plant_line", "n_batches_present", "batches"])
    rows = []
    for line, subset in data.groupby("plant_line"):
        batches = sorted(subset["batch"].dropna().unique(), key=lambda value: batch_sort_key(Path(str(value))))
        rows.append(
            {
                "plant_line": line,
                "n_batches_present": len(batches),
                "batches": "; ".join(batches),
            }
        )
    return pd.DataFrame(rows).sort_values(["n_batches_present", "plant_line"], ascending=[False, True])


def batch_manifest_long():
    rows = []
    for batch, details in BATCH_MANIFEST.items():
        date = details.get("date", "")
        seen = set()
        for status, targets in [("included", details.get("included", [])), ("excluded", details.get("excluded", []))]:
            for target in targets:
                line = canonical_line(target)
                if (batch, line) in seen:
                    continue
                seen.add((batch, line))
                global_reason = TARGETS_EXCLUDED_FROM_SCREEN.get(line, "")
                is_excluded = status == "excluded" or bool(global_reason)
                reason = "marked excluded in batch table" if status == "excluded" else ""
                if global_reason:
                    reason = "; ".join(part for part in [reason, global_reason] if part)
                rows.append(
                    {
                        "batch": batch,
                        "date": date,
                        "plant_line": line,
                        "display_line": display_line(line),
                        "manifest_status": "excluded" if is_excluded else "included",
                        "exclusion_reason": reason,
                    }
                )
    return pd.DataFrame(rows).sort_values(
        ["batch", "manifest_status", "display_line"],
        key=lambda series: series.map(lambda value: batch_sort_key(Path(str(value))) if str(value).startswith("Batch") else value),
    )


def batch_manifest_counts(manifest):
    if manifest.empty:
        return pd.DataFrame()
    return (
        manifest.groupby(["batch", "date"], as_index=False)
        .agg(
            n_targets_total=("plant_line", "nunique"),
            n_targets_included=("manifest_status", lambda values: int((values == "included").sum())),
            n_targets_excluded=("manifest_status", lambda values: int((values == "excluded").sum())),
            included_targets=(
                "display_line",
                lambda values: "; ".join(
                    manifest.loc[values.index[manifest.loc[values.index, "manifest_status"] == "included"], "display_line"]
                ),
            ),
            excluded_targets=(
                "display_line",
                lambda values: "; ".join(
                    manifest.loc[values.index[manifest.loc[values.index, "manifest_status"] == "excluded"], "display_line"]
                ),
            ),
        )
        .sort_values("batch", key=lambda series: series.map(lambda value: batch_sort_key(Path(str(value)))))
    )


def manifest_plot_check(manifest, plotted_data):
    rows = []
    batches = sorted(
        set(manifest["batch"].unique()) | (set(plotted_data["batch"].unique()) if not plotted_data.empty else set()),
        key=lambda value: batch_sort_key(Path(str(value))),
    )
    for batch in batches:
        manifest_subset = manifest.loc[manifest["batch"] == batch]
        manifest_included = set(manifest_subset.loc[manifest_subset["manifest_status"] == "included", "plant_line"])
        manifest_excluded = set(manifest_subset.loc[manifest_subset["manifest_status"] == "excluded", "plant_line"])
        plotted_lines = (
            set(plotted_data.loc[plotted_data["batch"] == batch, "plant_line"].unique()) if not plotted_data.empty else set()
        )
        rows.append(
            {
                "batch": batch,
                "manifest_included_count": len(manifest_included),
                "manifest_excluded_count": len(manifest_excluded),
                "plotted_line_count": len(plotted_lines),
                "manifest_included_missing_from_plotted": "; ".join(
                    display_line(line) for line in sorted(manifest_included - plotted_lines)
                ),
                "plotted_not_in_manifest_included": "; ".join(
                    display_line(line) for line in sorted(plotted_lines - manifest_included)
                ),
                "excluded_lines_removed_from_plot": "; ".join(
                    display_line(line) for line in sorted(manifest_excluded - plotted_lines)
                ),
                "excluded_lines_still_in_plot": "; ".join(display_line(line) for line in sorted(manifest_excluded & plotted_lines)),
            }
        )
    return pd.DataFrame(rows)


def explicit_batch_exclusions():
    manifest = batch_manifest_long()
    rows = [
        {
            "batch": row["batch"],
            "plant_line": row["plant_line"],
            "exclusion_source": row["exclusion_reason"] or "batch manifest",
        }
        for _, row in manifest.loc[manifest["manifest_status"] == "excluded"].iterrows()
    ]
    rows.extend(
        {
            "batch": entry["batch"],
            "plant_line": canonical_line(entry["plant_line"]),
            "exclusion_source": entry["exclusion_source"],
        }
        for entry in MANUAL_BATCH_EXCLUSIONS
    )
    if CONSOLIDATED_BATCHES.exists():
        xl = pd.ExcelFile(CONSOLIDATED_BATCHES)
        for sheet in xl.sheet_names:
            if not sheet.lower().replace(" ", "").startswith("batch"):
                continue
            batch_number = re.search(r"batch\s*(\d+)", sheet, flags=re.IGNORECASE)
            batch = f"Batch{batch_number.group(1)}" if batch_number else sheet.replace(" ", "")
            table = pd.read_excel(CONSOLIDATED_BATCHES, sheet_name=sheet, header=None, dtype=object)
            for value in table.to_numpy().ravel():
                if not isinstance(value, str) or "exclud" not in value.lower():
                    continue
                target = re.split(r"\bexclud", value, flags=re.IGNORECASE)[0]
                target = re.sub(r"[^A-Za-z0-9/]+$", "", target).strip()
                if not target:
                    continue
                rows.append(
                    {
                        "batch": batch,
                        "plant_line": canonical_line(target),
                        "exclusion_source": f"{CONSOLIDATED_BATCHES.name}:{sheet}:{value}",
                    }
                )
    if not rows:
        return pd.DataFrame(columns=["batch", "plant_line", "exclusion_source"])
    exclusions = pd.DataFrame(rows).drop_duplicates(["batch", "plant_line"])
    return exclusions.sort_values(["batch", "plant_line"]).reset_index(drop=True)


def apply_batch_exclusions(data, exclusions):
    if data.empty:
        return data
    batch = data["batch"].iloc[0]
    excluded_lines = set(PLOT_EXCLUDED_LINES)
    excluded_lines.update(
        exclusions.loc[exclusions["batch"] == batch, "plant_line"] if not exclusions.empty else []
    )
    excluded_lines.update(line for line in data["plant_line"].unique() if str(line).startswith("NBD"))
    if not excluded_lines:
        return data
    return data.loc[~data["plant_line"].isin(excluded_lines)].copy()


def read_target_metadata():
    columns = [
        "target_key",
        "Target (abbrv.)",
        "Target (full)",
        "Putative immunity role",
        "NbL Targets",
        "Information Source",
        "Plasmid Internal ID",
        "Plasmid ID",
        "Plasmid Source",
    ]
    if not VIGS_LIBRARY.exists():
        return pd.DataFrame(columns=columns)
    source = pd.read_excel(VIGS_LIBRARY, sheet_name="VIGSConstructs", header=2)
    source = source.dropna(subset=["VIGS Target"]).copy()
    source = source.loc[~source["VIGS Target"].astype(str).str.contains("non-target", case=False, na=False)].copy()
    rows = []
    for _, row in source.iterrows():
        key = canonical_line(row["VIGS Target"])
        full = clean_label(row.get("Full Names", ""))
        role = clean_label(row.get("Role", "")) or clean_label(row.get("Rough Functional Grouping", ""))
        reference = clean_label(row.get("Reference", ""))
        notes = clean_label(row.get("Notes", ""))
        rows.append(
            {
                "target_key": key,
                "Target (abbrv.)": display_line(key),
                "Target (full)": full,
                "Putative immunity role": role,
                "NbL Targets": clean_label(row.get("NbL (LAB360)", "")),
                "Information Source": reference,
                "Plasmid Internal ID": clean_label(row.get("Plasmid ID", "")),
                "Plasmid ID": f"pTRV2::{display_line(key)}",
                "Plasmid Source": notes,
            }
        )
    metadata = pd.DataFrame(rows)
    return metadata.drop_duplicates("target_key", keep="first")


def target_inventory_table(data, exclusions):
    metadata = read_target_metadata()
    manifest = batch_manifest_long()
    manifest_counts = manifest.rename(columns={"plant_line": "target_key"})
    included = manifest_counts.loc[manifest_counts["manifest_status"] == "included", ["target_key", "batch"]].drop_duplicates()
    excluded = manifest_counts.loc[manifest_counts["manifest_status"] == "excluded", ["target_key", "batch"]].drop_duplicates()
    counted = manifest_counts[["target_key", "batch"]].drop_duplicates()
    all_keys = sorted(set(metadata.get("target_key", [])) | set(counted["target_key"]))
    if metadata.empty:
        metadata = pd.DataFrame({"target_key": all_keys, "Target (abbrv.)": [display_line(key) for key in all_keys]})
    else:
        missing = [key for key in all_keys if key not in set(metadata["target_key"])]
        if missing:
            metadata = pd.concat(
                [
                    metadata,
                    pd.DataFrame(
                        {
                            "target_key": missing,
                            "Target (abbrv.)": [display_line(key) for key in missing],
                        }
                    ),
                ],
                ignore_index=True,
            )

    rows = []
    for _, row in metadata.iterrows():
        key = row["target_key"]
        included_batches = set(included.loc[included["target_key"] == key, "batch"])
        excluded_batches = set(excluded.loc[excluded["target_key"] == key, "batch"])
        counted_batches = sorted(
            counted.loc[counted["target_key"] == key, "batch"].unique(),
            key=lambda value: batch_sort_key(Path(str(value))),
        )
        excluded_batch_list = sorted(excluded_batches, key=lambda value: batch_sort_key(Path(str(value))))
        screen_exclusion_note = TARGETS_EXCLUDED_FROM_SCREEN.get(key, "")
        screen_status = "excluded from target list" if screen_exclusion_note else "included"
        output = row.drop(labels=["target_key"]).to_dict()
        output.update(
            {
                "Screen status": screen_status,
                "Screen exclusion note": screen_exclusion_note,
                "Batches counted so far": "; ".join(counted_batches),
                "VIGS batch instances so far": len(counted_batches),
                "Excluded batch instances": len(excluded_batch_list),
                "Excluded batches": "; ".join(excluded_batch_list),
                "_sort_excluded": 1 if screen_exclusion_note else 0,
            }
        )
        rows.append(output)
    table = pd.DataFrame(rows)
    table = table.sort_values(["_sort_excluded", "Target (abbrv.)"], ascending=[True, True])
    return table.drop(columns=["_sort_excluded"]).reset_index(drop=True)


def compiled_gus_data(data):
    if data.empty:
        return pd.DataFrame(), pd.DataFrame()
    excluded = pd.DataFrame(columns=["batch", "control_line", "exclusion_reason"])
    ready = data.loc[(data["control_line"] == "GUS") & (data["relative_to_control"] > 0)].copy()
    return ready, excluded


def compiled_batch_means(data):
    if data.empty:
        return pd.DataFrame()
    return (
        data.groupby(["batch", "plant_line"], as_index=False)
        .agg(
            batch_mean_relative_gfp=("relative_to_control", "mean"),
            batch_sem_relative_gfp=("relative_to_control", sem),
            n_replicates=("relative_to_control", "count"),
        )
        .sort_values(["plant_line", "batch"])
    )


def run_compiled_batch_adjusted_stats(data):
    ready, excluded = compiled_gus_data(data)
    if ready.empty or ready["batch"].nunique() < 2 or ready["plant_line"].nunique() < 2:
        empty_terms = pd.DataFrame()
        empty_contrasts = pd.DataFrame()
        return empty_terms, empty_contrasts, excluded

    ready = ready.copy()
    ready["log2_relative_gfp"] = np.log2(ready["relative_to_control"])
    model = ols(
        'log2_relative_gfp ~ C(plant_line, Treatment(reference="GUS")) + C(batch)',
        data=ready,
    ).fit()
    robust = model.get_robustcov_results(cov_type="HC3")
    params = pd.Series(robust.params, index=model.params.index)
    p_values = pd.Series(robust.pvalues, index=model.params.index)
    conf_int = pd.DataFrame(robust.conf_int(), index=model.params.index, columns=["ci_low", "ci_high"])
    model_terms = pd.DataFrame(
        {
            "term": params.index,
            "estimate_log2": params.values,
            "p_value_hc3": p_values.values,
            "ci_low_log2": conf_int["ci_low"].values,
            "ci_high_log2": conf_int["ci_high"].values,
            "model": "OLS log2(relative GFP) ~ plant line + batch, HC3 robust SE",
        }
    )

    contrast_rows = []
    for term, estimate in params.items():
        if "C(plant_line" not in term or "[T." not in term:
            continue
        match = re.search(r"\[T\.(.*)\]", term)
        if not match:
            continue
        line = match.group(1)
        contrast_rows.append(
            {
                "plant_line": line,
                "estimate_log2_vs_GUS": float(estimate),
                "fold_change_vs_GUS": float(2**estimate),
                "ci_low_fold_change": float(2 ** conf_int.loc[term, "ci_low"]),
                "ci_high_fold_change": float(2 ** conf_int.loc[term, "ci_high"]),
                "p_raw_hc3": float(p_values.loc[term]),
                "test": "Batch-adjusted OLS term vs GUS on log2(relative GFP), HC3 robust SE",
            }
        )
    contrasts = pd.DataFrame(contrast_rows)
    if not contrasts.empty:
        _, adjusted, _, _ = multipletests(contrasts["p_raw_hc3"], method="holm")
        contrasts["p_holm"] = adjusted
        contrasts["stars"] = [p_to_stars(p) for p in adjusted]
    return model_terms, contrasts, excluded


def plot_compiled_across_batches(data, contrasts):
    ready, excluded = compiled_gus_data(data)
    batch_means = compiled_batch_means(ready)
    if batch_means.empty:
        return batch_means, excluded

    summary = (
        batch_means.groupby("plant_line", as_index=False)
        .agg(
            mean_relative_gfp=("batch_mean_relative_gfp", "mean"),
            sem_relative_gfp=("batch_mean_relative_gfp", sem),
            n_batches=("batch", "nunique"),
        )
    )
    order = ["GUS"] + [
        line
        for line in summary.sort_values(["mean_relative_gfp", "plant_line"], ascending=[True, True])["plant_line"]
        if line != "GUS"
    ]
    positions = np.arange(len(order), dtype=float)
    position_map = dict(zip(order, positions))
    fig_width = max(18.0, min(30.0, 5.0 + len(order) * 0.58))
    fig, ax = plt.subplots(
        figsize=(fig_width, (fig_width / np.sqrt(2)) * (2 / 3)),
        constrained_layout=False,
    )
    ymax = 0
    annotation_top = 0
    star_map = dict(zip(contrasts.get("plant_line", []), contrasts.get("stars", []))) if not contrasts.empty else {}

    for line_i, line in enumerate(order):
        row = summary.loc[summary["plant_line"] == line].iloc[0]
        x = position_map[line]
        color = CONTROL_COLOR if line == "GUS" else BASE_COLOR
        ax.bar(
            x,
            row["mean_relative_gfp"],
            width=0.70,
            color=color,
            alpha=BAR_ALPHA,
            edgecolor="black",
            linewidth=BAR_LINEWIDTH,
            zorder=2,
        )
        ax.errorbar(
            x,
            row["mean_relative_gfp"],
            yerr=row["sem_relative_gfp"],
            fmt="none",
            ecolor=ERRORBAR_COLOR,
            elinewidth=1.1,
            capsize=3.0,
            capthick=1.1,
            zorder=5,
        )
        # Show every measured plant, while retaining batch means as the unit
        # summarized by the bars and error bars.
        points = ready.loc[ready["plant_line"] == line, "relative_to_control"].dropna()
        rng = np.random.default_rng(202607 + line_i)
        jitter = rng.uniform(-0.20, 0.20, len(points)) if len(points) > 1 else np.array([0.0])
        ax.scatter(
            x + jitter,
            points,
            # Previous diameter was sqrt(34) pt; s=62 adds approximately 2 pt.
            s=62,
            facecolor="#3E3E3E" if line == "GUS" else "#315F78",
            edgecolor="none",
            alpha=1.0,
            zorder=4,
        )
        ymax = max(ymax, float(points.max()), float(row["mean_relative_gfp"] + row["sem_relative_gfp"]))
        if line in star_map and star_map[line] != "ns":
            star_y = max(
                points.max(), row["mean_relative_gfp"] + row["sem_relative_gfp"]
            ) + max(ymax * 0.07, 0.08)
            ax.text(
                x,
                star_y,
                star_map[line],
                ha="center",
                va="bottom",
                fontsize=33.5,
                color=TEXT_COLOR,
            )
            annotation_top = max(annotation_top, float(star_y))

    ax.axhline(
        1, color="#C43B32", linestyle=(0, (7, 4)), linewidth=3.2,
        alpha=1.0, zorder=3,
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [display_line(line) for line in order], rotation=60, ha="right",
        fontsize=26.5,
    )
    ax.set_ylabel("GFP Expression (GUS-relative)", fontsize=27.5)
    upper_limit = max(ymax * 1.08, annotation_top * 1.06) if ymax else 1
    ax.set_ylim(0, upper_limit)
    ax.set_xlim(-0.35, len(order) - 0.65)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax)
    ax.tick_params(axis="y", labelsize=26.5)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.94, bottom=0.33)
    save(fig, "compiled_across_batches_relative_gfp")
    plt.close(fig)
    return batch_means, excluded


def plot_all_batches_panel(data, pairwise, crossbatch_image=None, output_name="all_batches_summary_panel"):
    if data.empty:
        return
    candidate_batches = []
    for batch, subset in data.groupby("batch", sort=False):
        non_control_lines = [line for line in subset["plant_line"].unique() if line != "GUS"]
        if len(non_control_lines) == 0:
            continue
        candidate_batches.append(batch)
    candidate_batches = sorted(candidate_batches, key=lambda value: batch_sort_key(Path(str(value))))
    if not candidate_batches:
        return

    n_cols = 3
    slot_count = len(candidate_batches) + (1 if "Batch1" in candidate_batches else 0)
    n_rows = int(np.ceil(slot_count / n_cols))
    fig_width = 38.0
    image = readable_image(crossbatch_image) if crossbatch_image else None
    image_rows = 1 if image is not None else 0
    fig_height = max(9.0, 7.8 * n_rows + (13.0 if image is not None else 0))
    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=False)
    grid = fig.add_gridspec(
        n_rows + image_rows,
        n_cols,
        height_ratios=[1.0] * n_rows + ([2.25] if image is not None else []),
        hspace=0.78,
        wspace=0.33,
    )
    axes = []
    slot = 0
    for batch in candidate_batches:
        span = 2 if batch == "Batch1" else 1
        row = slot // n_cols
        col = slot % n_cols
        if col + span > n_cols:
            slot += n_cols - col
            row = slot // n_cols
            col = 0
        axes.append(fig.add_subplot(grid[row, col : col + span]))
        slot += span
    axes = np.asarray(axes)

    if image is not None:
        image_ax = fig.add_subplot(grid[n_rows, :])
        image_ax.imshow(np.asarray(image))
        image_ax.set_axis_off()

    for ax, batch in zip(axes, candidate_batches):
        subset = data.loc[data["batch"] == batch].copy()
        batch_pairwise = pairwise.loc[pairwise["batch"] == batch].copy() if not pairwise.empty else pd.DataFrame()
        control = subset["control_line"].dropna().iloc[0] if subset["control_line"].notna().any() else "GUS"
        plot_batch_bar(ax, subset, batch, control, batch_pairwise, master_style=True)
    if image is not None:
        fig.subplots_adjust(left=0.06, right=0.99, top=0.985, bottom=0.035)
    else:
        fig.subplots_adjust(left=0.085, right=0.99, top=0.975, bottom=0.15, wspace=0.38, hspace=1.05)
    shared_label_size = 40.5
    fig.supylabel(
        "GFP expression (relative to GUS control)", x=0.020,
        fontsize=shared_label_size, fontweight="bold", color=TEXT_COLOR,
    )
    fig.supxlabel(
        "VIGS target", y=0.012, fontsize=shared_label_size,
        fontweight="bold", color=TEXT_COLOR,
    )
    guide_color = "#9A9A9A"
    fig.add_artist(Line2D(
        [0.052, 0.052], [0.145, 0.975], transform=fig.transFigure,
        color=guide_color, linewidth=1.2, solid_capstyle="butt",
    ))
    fig.add_artist(Line2D(
        [0.085, 0.99], [0.060, 0.060], transform=fig.transFigure,
        color=guide_color, linewidth=1.2, solid_capstyle="butt",
    ))
    save(fig, output_name)
    plt.close(fig)


# %%
def main():
    inventory, flags = inventory_batches()
    manifest = batch_manifest_long()
    manifest_counts = batch_manifest_counts(manifest)
    batch_exclusions = explicit_batch_exclusions()
    inventory.to_excel(LOG_DIR / "vigs_screen_batch_inventory.xlsx", index=False)
    manifest.to_excel(LOG_DIR / "vigs_screen_batch_manifest.xlsx", index=False)
    manifest_counts.to_excel(LOG_DIR / "vigs_screen_batch_manifest_counts.xlsx", index=False)
    if not flags.empty:
        flags.to_excel(LOG_DIR / "vigs_screen_flagged_cases.xlsx", index=False)

    raw_tables = []
    metadata_rows = []
    parse_flags = [] if flags.empty else flags.to_dict("records")
    plot_rows = []
    omnibus_tables = []
    pairwise_tables = []

    for _, meta in inventory.iterrows():
        if meta["status"] != "ready" and meta["batch"] not in BATCH_SOURCE_OVERRIDES:
            continue
        batch = meta["batch"]
        override = BATCH_SOURCE_OVERRIDES.get(batch)
        path = Path(override["path"]) if override else Path(meta["primary_excel"])
        forced_sheet = override["sheet"] if override else None
        if override:
            override_details = {
                "Batch3": (
                    f"Used {path}, sheet {forced_sheet}, as the single source workbook "
                    "for Batch3; its combined-date filename is retained as a provenance warning."
                ),
                "Batch4": (
                    f"Used {path}, sheet {forced_sheet}. The workbook's All sheet contains "
                    "the genuine Batch4 measurements; its Raw sheet contains the Batch5 "
                    "measurements and must not be selected."
                ),
                "Batch9": (
                    f"Used {path}, sheet {forced_sheet}, retaining only the 20241204 "
                    "measurements with same-day GUS. The 20241115 workbook and appended "
                    "20241217 LIK1 rows lack same-day GUS controls."
                ),
            }[batch]
            parse_flags.append(
                {
                    "batch": batch,
                    "flag_type": "source_override",
                    "details": override_details,
                }
            )
        try:
            data = parse_excel_table(path, forced_sheet=forced_sheet)
            data.insert(0, "batch", batch)
            if batch == "Batch9":
                data = data.loc[
                    ~(
                        ((data["plant_line"] == "LIK1") & (data["replicate"] > 6))
                        | (data["plant_line"] == "LIK1 (2)")
                    )
                ].copy()
            data = apply_batch_exclusions(data, batch_exclusions)
            data, control = add_relative_values(data)
            if control is None:
                parse_flags.append(
                    {"batch": batch, "flag_type": "batch_not_processed", "details": "No GUS or WT control found."}
                )
                continue
            data["source_excel"] = str(path)
            data["dates_detected"] = meta["dates_detected"]
            raw_tables.append(data)

            omnibus, pairwise = run_batch_stats(data, batch, control)
            omnibus_tables.append(omnibus)
            if not pairwise.empty:
                pairwise_tables.append(pairwise)

            plot_rows.append(plot_batch_figure(data, meta, pairwise))
            metadata_rows.append(
                {
                    "batch": batch,
                    "source_excel": str(path),
                    "source_sheet": data["source_sheet"].iloc[0],
                    "parse_mode": data["parse_mode"].iloc[0],
                    "source_value_column": data["source_value_column"].iloc[0],
                    "control_line": control,
                    "dates_detected": meta["dates_detected"],
                    "plant_lines": "; ".join(sorted(data["plant_line"].unique(), key=sort_key_label)),
                    "n_replicates_total": len(data),
                    "n_plant_lines": data["plant_line"].nunique(),
                    "pptx_files": meta["pptx_files"],
                    "prism_files": meta["prism_files"],
                }
            )
        except Exception as exc:
            parse_flags.append(
                {
                    "batch": batch,
                    "flag_type": "parse_failed",
                    "details": f"{type(exc).__name__}: {exc}",
                }
            )

    combined = pd.concat(raw_tables, ignore_index=True) if raw_tables else pd.DataFrame()
    metadata = pd.DataFrame(metadata_rows)
    plot_log = pd.DataFrame(plot_rows)
    omnibus = pd.concat(omnibus_tables, ignore_index=True) if omnibus_tables else pd.DataFrame()
    pairwise = pd.concat(pairwise_tables, ignore_index=True) if pairwise_tables else pd.DataFrame()
    flagged = pd.DataFrame(parse_flags)
    repeat_counts = plant_line_batch_counts(combined)
    manifest_check = manifest_plot_check(manifest, combined)
    target_table = target_inventory_table(combined, batch_exclusions)
    compiled_model_terms, compiled_contrasts, compiled_excluded = run_compiled_batch_adjusted_stats(combined)
    compiled_batch_mean_table, compiled_excluded_from_plot = plot_compiled_across_batches(
        combined,
        compiled_contrasts,
    )
    plot_all_batches_panel(combined, pairwise)
    crossbatch_image = PANEL_IMAGE_DIR / "CrossBatch.png"
    if crossbatch_image.exists():
        plot_all_batches_panel(
            combined,
            pairwise,
            crossbatch_image=crossbatch_image,
            output_name="all_batches_summary_panel_with_crossbatch_images",
        )
    if not compiled_excluded_from_plot.empty:
        compiled_excluded = pd.concat(
            [compiled_excluded, compiled_excluded_from_plot],
            ignore_index=True,
        ).drop_duplicates()

    combined.to_excel(LOG_DIR / "vigs_screen_normalized_long.xlsx", index=False)
    metadata.to_excel(LOG_DIR / "vigs_screen_metadata_compiled.xlsx", index=False)
    manifest_check.to_excel(LOG_DIR / "vigs_screen_batch_manifest_plot_check.xlsx", index=False)
    plot_log.to_excel(LOG_DIR / "vigs_screen_image_extraction_log.xlsx", index=False)
    omnibus.to_excel(LOG_DIR / "vigs_screen_kruskal_wallis.xlsx", index=False)
    pairwise.to_excel(LOG_DIR / "vigs_screen_control_pairwise_mannwhitney.xlsx", index=False)
    target_table.to_excel(LOG_DIR / "vigs_screen_target_table.xlsx", index=False)
    target_table.to_excel(LOG_DIR / "vigs_screen_plant_line_batch_counts.xlsx", index=False)
    batch_exclusions.to_excel(LOG_DIR / "vigs_screen_batch_exclusions.xlsx", index=False)
    compiled_batch_mean_table.to_excel(LOG_DIR / "vigs_screen_compiled_batch_means.xlsx", index=False)
    compiled_model_terms.to_excel(LOG_DIR / "vigs_screen_compiled_batch_adjusted_model_terms.xlsx", index=False)
    compiled_contrasts.to_excel(LOG_DIR / "vigs_screen_compiled_batch_adjusted_line_contrasts.xlsx", index=False)
    compiled_excluded.to_excel(LOG_DIR / "vigs_screen_compiled_excluded_batches.xlsx", index=False)
    flagged.to_excel(LOG_DIR / "vigs_screen_flagged_cases.xlsx", index=False)

    print(f"Inventoried {len(inventory)} non-scrap batch folders.")
    print(f"Processed {combined['batch'].nunique() if not combined.empty else 0} batches into {len(combined)} rows.")
    print(f"Flagged {len(flagged)} cases for review.")
    print(f"Figures saved to: {SAVE_DIR}")
    print(f"Logs saved to: {LOG_DIR}")


if __name__ == "__main__":
    main()
