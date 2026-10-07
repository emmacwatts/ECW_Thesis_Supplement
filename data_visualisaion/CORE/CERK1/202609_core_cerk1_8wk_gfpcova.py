from pathlib import Path
import io
import json
import os
import sys
import zipfile

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

SOURCE_FILE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "6. PosHits/CORE/8wk_CORE/8wkCORE_wCERK/20250521_8wk_CORE/"
    "8wk_CORE_20250521.prism"
)

STEM = Path(__file__).stem
OUT = HERE / f"{STEM}_figures"
LOG = HERE / f"{STEM}_logs"
OUT.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

ORDER = ["WT", "cerk1-1"]
CONDITIONS = ["No p19", "+ p19"]
SOURCE_GENOTYPES = {"cas9": "WT", "cerk136": "cerk1-1"}
BAR_COLOR = "#86A5C1"
FORMATS = ("png", "pdf", "svg")

plt.rcParams.update({
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 18,
    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})


def _read_json(archive, member):
    return json.loads(archive.read(member))


def load_prism_data():
    """Read the two Prism source tables and retain only Cas9 and cerk1-3-6."""
    rows = []
    with zipfile.ZipFile(SOURCE_FILE) as archive:
        document = _read_json(archive, "document.json")
        attributes = document["sheetAttributesMap"]
        for sheet_uid in document["sheets"]["data"]:
            sheet = _read_json(archive, f"data/sheets/{sheet_uid}/sheet.json")
            source_condition = attributes[sheet_uid]["title"]
            condition = "+ p19" if source_condition == "w p19" else "No p19"
            table = sheet["table"]
            values = pd.read_csv(
                io.BytesIO(archive.read(f"data/tables/{table['uid']}/data.csv")),
                header=None,
            )
            for column, dataset_uid in enumerate(table["dataSets"]):
                dataset = _read_json(archive, f"data/sets/{dataset_uid}.json")
                source_genotype = dataset["title"]
                if source_genotype not in SOURCE_GENOTYPES:
                    continue
                for replicate, value in enumerate(values.iloc[:, column].dropna(), start=1):
                    rows.append({
                        "condition": condition,
                        "source_condition": source_condition,
                        "genotype": SOURCE_GENOTYPES[source_genotype],
                        "source_genotype": source_genotype,
                        "replicate": replicate,
                        "normalised_gfp": float(value),
                        "source_file": str(SOURCE_FILE),
                        "source_sheet": source_condition,
                    })
    return pd.DataFrame(rows)


def analyse(data):
    summary = (
        data.groupby(["condition", "genotype"], sort=False)["normalised_gfp"]
        .agg(n="size", mean="mean", sd="std", sem="sem")
        .reset_index()
    )
    tests = []
    for condition in CONDITIONS:
        subset = data.loc[data.condition.eq(condition)]
        wt = subset.loc[subset.genotype.eq("WT"), "normalised_gfp"]
        cerk = subset.loc[subset.genotype.eq("cerk1-1"), "normalised_gfp"]
        test = stats.ttest_ind(wt, cerk, equal_var=False, alternative="two-sided")
        tests.append({
            "condition": condition,
            "comparison": "WT vs cerk1-1",
            "test": "two-sided Welch t-test",
            "wt_n": len(wt),
            "cerk_n": len(cerk),
            "wt_mean": wt.mean(),
            "cerk_mean": cerk.mean(),
            "t": test.statistic,
            "df": getattr(test, "df", np.nan),
            "p": test.pvalue,
        })
    return summary, pd.DataFrame(tests)


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    ax.spines[["top", "right"]].set_visible(False)


def add_panel(ax, data, summary, test_row, condition, show_ylabel):
    panel_data = data.loc[data.condition.eq(condition)]
    panel_summary = (
        summary.loc[summary.condition.eq(condition)]
        .set_index("genotype").reindex(ORDER).reset_index()
    )
    x = np.arange(2)
    means = panel_summary["mean"].to_numpy() / 1000
    sems = panel_summary["sem"].to_numpy() / 1000
    ax.bar(
        x, means, yerr=sems, width=0.6, color=BAR_COLOR, alpha=0.9,
        edgecolor="#66727D", linewidth=1.25, capsize=5,
        error_kw={"elinewidth": 1.4, "ecolor": TEXT_COLOR, "capthick": 1.4},
    )
    rng = np.random.default_rng(260928)
    for i, genotype in enumerate(ORDER):
        values = panel_data.loc[panel_data.genotype.eq(genotype), "normalised_gfp"].to_numpy() / 1000
        jitter = rng.uniform(-0.105, 0.105, len(values))
        ax.scatter(
            np.full(len(values), i) + jitter, values, s=32,
            facecolor="#587C9C", edgecolor="white", linewidth=0.55,
            alpha=0.8, zorder=4,
        )
    ax.set_xticks(x, ["WT", "cerk1-1"])
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_title(condition, fontsize=20, pad=14)
    if show_ylabel:
        ax.set_ylabel("Relative GFP intensity\n(a.u. × 10³)")
    ax.set_xlim(-0.62, 1.62)
    style_axis(ax)

    data_max = max((panel_data.normalised_gfp / 1000).max(), (means + sems).max())
    bracket_y = data_max + max(means) * 0.19
    tick = max(means) * 0.045
    ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
            color=TEXT_COLOR, lw=1.2, clip_on=False)
    p = float(test_row["p"])
    label = "ns" if p >= 0.05 else f"P = {p:.3g}"
    ax.text(0.5, bracket_y + tick * 0.25, label, ha="center", va="bottom", fontsize=16)
    ax.set_ylim(0, bracket_y + max(means) * 0.17)


def save_logs(data, summary, tests):
    with pd.ExcelWriter(LOG / "core_cerk1_8wk_processed_data_and_statistics.xlsx") as writer:
        data.to_excel(writer, sheet_name="source_values", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        tests.to_excel(writer, sheet_name="welch_tests", index=False)
    (LOG / "core_cerk1_8wk_methods_note.txt").write_text(
        "Eight-week GFP data were read directly from the Prism source tables 'w p19' "
        "and 'no p19'. Only source columns cas9 and cerk136 (the cerk1-3-6 line) were "
        "retained; 704-6, 705-1, and cerk146 were excluded. Cas9 is displayed as WT "
        "and cerk1-3-6 as italic cerk1-1. Bars show means +/- SEM and all five source "
        "values. Two-sided Welch t-tests compare WT with cerk1-1 separately within "
        "each p19 condition.\n"
    )
    (LOG / "core_cerk1_8wk_figure_caption.txt").write_text(
        "GFP accumulation in 8-week-old WT and cerk1-1 plants in the absence "
        "(No p19) or presence (+ p19) of the silencing suppressor p19. Relative "
        "GFP intensity is shown for WT (Cas9) and cerk1-1 (source line "
        "cerk1-3-6). Bars represent mean +/- SEM; points represent individual "
        "biological replicates (n = 5). Comparisons between genotypes within "
        "each treatment were performed using two-sided Welch's t-tests. ns, "
        "not significant.\n"
    )


def make_figure(data, summary, tests):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 6.0), sharey=True, facecolor="white")
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.15, top=0.86, wspace=0.28)
    for i, (ax, condition) in enumerate(zip(axes, CONDITIONS)):
        test_row = tests.loc[tests.condition.eq(condition)].iloc[0]
        add_panel(ax, data, summary, test_row, condition, show_ylabel=(i == 0))
    for fmt in FORMATS:
        fig.savefig(OUT / f"{STEM}_wt_vs_cerk1_1_p19.{fmt}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    data = load_prism_data()
    summary, tests = analyse(data)
    save_logs(data, summary, tests)
    make_figure(data, summary, tests)
    print(summary.to_string(index=False))
    print(tests[["condition", "t", "df", "p"]].to_string(index=False))
