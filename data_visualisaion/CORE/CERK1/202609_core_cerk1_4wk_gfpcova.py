from pathlib import Path
import io
import os
import sys
import zipfile

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats

STYLE_DIR = HERE.parents[1] / "PRp27" / "PRp27_timecourse_202606"
if str(STYLE_DIR) not in sys.path:
    sys.path.insert(0, str(STYLE_DIR))
from figure_styles import ANALYSIS_RCPARAMS, FONT_STACK, TEXT_COLOR, apply_axis_style

SOURCE_DIR = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/"
    "1.General/ReportsandPapers/IssyVIGS/NPR1_FollowOn_Exps/4wk_GFPCOVA"
)
GFP_FILE = SOURCE_DIR / "GFP" / "WTvNPR1.xlsx"
ELISA_FILE = SOURCE_DIR / "COVA_ELISA" / "cova_elisa_npr1cerkwt.xlsx"
SOURCE_PPTX = SOURCE_DIR / "IssyVIGS_4wk_GFPCOVA_finalfig.pptx"

STEM = Path(__file__).stem
OUT = HERE / f"{STEM}_figures"
LOG = HERE / f"{STEM}_logs"
OUT.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

ORDER = ["WT", "cerk1-1"]
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


def load_gfp():
    data = pd.read_excel(GFP_FILE, sheet_name="Raw", header=2)
    data = data.rename(columns={
        "Plant Type": "source_genotype", "Repeat": "repeat",
        "Replicate": "replicate", "Background": "background",
        "GFP": "gfp", "Normalised GFP": "normalised_gfp",
    })
    data = data.loc[
        ((data.source_genotype == "WT") & (data.repeat == 1))
        | (data.source_genotype == "cerk1-1382")
    ].copy()
    data["genotype"] = data.source_genotype.replace({"cerk1-1382": "cerk1-1"})
    data["source_file"] = str(GFP_FILE)
    data["source_sheet"] = "Raw"
    return data.reset_index(drop=True)


def load_elisa():
    raw = pd.read_excel(ELISA_FILE, sheet_name="Processed", header=None)
    data = raw.iloc[41:113, :4].copy()
    data.columns = ["source_genotype", "infiltration", "input_ug_ml", "a450"]
    data["input_ug_ml"] = pd.to_numeric(data.input_ug_ml, errors="coerce")
    data["a450"] = pd.to_numeric(data.a450, errors="coerce")
    data = data.loc[
        data.source_genotype.isin(["WT", "cerk1-"])
        & data.infiltration.eq("COVA")
        & data.input_ug_ml.eq(0.16)
    ].copy()
    data["genotype"] = data.source_genotype.replace({"cerk1-": "cerk1-1"})
    data["replicate"] = data.groupby("genotype", sort=False).cumcount() + 1
    data["source_file"] = str(ELISA_FILE)
    data["source_sheet"] = "Processed"
    return data.reset_index(drop=True)


def analyse(data, value, assay):
    summary = (
        data.groupby("genotype")[value]
        .agg(n="size", mean="mean", sd="std", sem="sem")
        .reindex(ORDER).reset_index()
    )
    wt = data.loc[data.genotype.eq("WT"), value]
    cerk = data.loc[data.genotype.eq("cerk1-1"), value]
    test = stats.ttest_ind(wt, cerk, equal_var=False, alternative="two-sided")
    comparison = pd.DataFrame([{
        "assay": assay,
        "comparison": "WT vs cerk1-1",
        "test": "two-sided Welch t-test",
        "wt_n": len(wt), "cerk_n": len(cerk),
        "wt_mean": wt.mean(), "cerk_mean": cerk.mean(),
        "t": test.statistic, "df": getattr(test, "df", np.nan), "p": test.pvalue,
    }])
    return summary, comparison


def style_axis(ax):
    ax.set_facecolor("#E9EFF6")
    ax.grid(axis="y", color="white", linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    apply_axis_style(ax)
    ax.spines[["top", "right"]].set_visible(False)


def add_bar_panel(ax, data, summary, test_row, value, ylabel, scale=1):
    x = np.arange(2)
    means = summary["mean"].to_numpy() / scale
    sems = summary["sem"].to_numpy() / scale
    ax.bar(
        x, means, yerr=sems, width=0.6,
        color=BAR_COLOR, alpha=0.9,
        edgecolor="#66727D", linewidth=1.25, capsize=5,
        error_kw={"elinewidth": 1.4, "ecolor": TEXT_COLOR, "capthick": 1.4},
    )
    rng = np.random.default_rng(260928)
    for i, genotype in enumerate(ORDER):
        values = data.loc[data.genotype.eq(genotype), value].to_numpy() / scale
        jitter = rng.uniform(-0.105, 0.105, len(values))
        ax.scatter(
            np.full(len(values), i) + jitter, values, s=32,
            facecolor="#587C9C", edgecolor="white", linewidth=0.55,
            alpha=0.8, zorder=4,
        )
        top = max(values.max(), means[i] + sems[i])
        ax.text(i, top + max(means) * 0.055, f"{len(values)}", ha="center",
                va="bottom", fontsize=15, color="#666666")
    ax.set_xticks(x, ["WT", "cerk1-1"])
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.62, 1.62)
    style_axis(ax)

    data_max = max((data[value] / scale).max(), (means + sems).max())
    bracket_y = data_max + max(means) * 0.19
    tick = max(means) * 0.045
    ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
            color=TEXT_COLOR, lw=1.2, clip_on=False)
    label = "ns" if float(test_row["p"]) >= 0.05 else f"P = {float(test_row['p']):.3g}"
    ax.text(0.5, bracket_y + tick * 0.25, label, ha="center", va="bottom", fontsize=16)
    ax.set_ylim(0, bracket_y + max(means) * 0.17)


def load_source_western_panel():
    with zipfile.ZipFile(SOURCE_PPTX) as archive:
        return Image.open(io.BytesIO(archive.read("ppt/media/image2.png"))).convert("RGB")


def stitch_lane_blocks(image, y0, y1):
    # Retain only the two GFP-treated WT lanes and the two GFP-treated CERK lanes.
    # The source blocks are joined directly; a dashed line discloses the splice.
    # Stop before the heavy right-hand frame line in the presentation asset.
    blocks = [(983, 1232), (1480, 1718)]
    pieces = [np.asarray(image.crop((x0, y0, x1, y1))) for x0, x1 in blocks]
    return np.concatenate(pieces, axis=1)


def add_western_panel(ax):
    source = load_source_western_panel()
    # A tighter vertical crop preserves the useful molecular-weight region and
    # avoids artificially stretching a large blank part of the source panel.
    blot = stitch_lane_blocks(source, 610, 820)
    ponceau = stitch_lane_blocks(source, 965, 1115)
    ax.set_axis_off()
    # Bounds are matched to the cropped image proportions so aspect='auto'
    # displays the scans without visible distortion or bbox expansion.
    blot_ax = ax.inset_axes([0.195, 0.44, 0.55, 0.32])
    ponceau_ax = ax.inset_axes([0.195, 0.07, 0.55, 0.22])
    for image_ax, panel in [(blot_ax, blot), (ponceau_ax, ponceau)]:
        image_ax.imshow(panel, aspect="auto")
        image_ax.set_xticks([]); image_ax.set_yticks([])
        for spine in image_ax.spines.values():
            spine.set_color(TEXT_COLOR); spine.set_linewidth(0.9)
        image_ax.axvline(0.5 * panel.shape[1], color=TEXT_COLOR, lw=1.0, ls="--")

    centers = [0.3325, 0.6075]
    groups = ["WT", "cerk1-1"]
    for x, group in zip(centers, groups):
        ax.text(x, 1.015, group, ha="center", va="top", fontsize=18,
                fontstyle="italic" if group != "WT" else "normal")
        ax.plot([x - 0.11, x + 0.11], [0.905, 0.905], transform=ax.transAxes,
                color=TEXT_COLOR, lw=1.2, clip_on=False)
        ax.text(x - 0.055, 0.86, "1", ha="center", va="top", fontsize=16, zorder=10)
        ax.text(x + 0.055, 0.86, "2", ha="center", va="top", fontsize=16, zorder=10)
    ax.text(0.745, 0.405, r"$\alpha$-GFP", ha="right", va="top", fontsize=18)
    ax.text(0.745, 0.045, "Ponceau S", ha="right", va="top", fontsize=18,
            fontstyle="italic")
    # Keep the underline artists from changing the host axis data limits, which
    # would otherwise inflate bbox_inches='tight' exports.
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def save_logs(gfp, elisa, gfp_summary, elisa_summary, tests):
    with pd.ExcelWriter(LOG / "core_cerk1_4wk_processed_data_and_statistics.xlsx") as writer:
        gfp.to_excel(writer, sheet_name="gfp_source_values", index=False)
        gfp_summary.to_excel(writer, sheet_name="gfp_summary", index=False)
        elisa.to_excel(writer, sheet_name="elisa_source_values", index=False)
        elisa_summary.to_excel(writer, sheet_name="elisa_summary", index=False)
        tests.to_excel(writer, sheet_name="welch_tests", index=False)
    (LOG / "core_cerk1_4wk_methods_note.txt").write_text(
        "Four-week WT and cerk1-1 data were retained; NPR1 and all other plant "
        "types were excluded. GFP values are source-workbook Normalised GFP "
        "(GFP minus background) from WT repeat 1 and the matching cerk1-1382 block. "
        "COVA ELISA values are A450 readings from COVA-infiltrated samples at "
        "0.16 micrograms/mL input. Bars show means +/- SEM and all source values. "
        "Two-sided Welch t-tests are logged for WT versus cerk1-1. The western "
        "panel is reconstructed from the original presentation-embedded scan; only "
        "panel displays only the GFP-treated WT and CERK lanes, with the removed-lane "
        "join marked by a dashed vertical line. No figure caption is included.\n"
    )


def make_figure(gfp, elisa, gfp_summary, elisa_summary, tests):
    fig = plt.figure(figsize=(13.8, 8.8), facecolor="white")
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1.08, 0.92], height_ratios=[1.0, 0.9],
        left=0.11, right=0.97, bottom=0.10, top=0.91,
        hspace=0.34, wspace=0.34,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_b = fig.add_subplot(gs[:, 1])
    # Match the width and left edge of panel A to the scan area in panel C.
    left_cell = ax_a.get_position()
    aligned_left = left_cell.x0 + 0.195 * left_cell.width
    aligned_width = 0.55 * left_cell.width
    ax_a.set_position([aligned_left, left_cell.y0, aligned_width, left_cell.height])
    # Reclaim the unused part of the original left grid cell so the full-height
    # right panel sits close to A/C without crowding its y-axis label.
    b_cell = ax_b.get_position()
    b_left = aligned_left + aligned_width + 0.11
    ax_b.set_position([b_left, b_cell.y0, aligned_width, b_cell.height])
    add_bar_panel(ax_a, gfp, gfp_summary, tests.iloc[0], "normalised_gfp",
                  "Relative GFP intensity\n(a.u. × 10³)", scale=1000)
    add_bar_panel(ax_b, elisa, elisa_summary, tests.iloc[1], "a450", "COVA ELISA (A450 nm)")
    add_western_panel(ax_c)
    a_box = ax_a.get_position()
    b_box = ax_b.get_position()
    c_box = ax_c.get_position()
    letter_specs = [
        (aligned_left - 0.067, a_box.y1 + 0.038, "A"),
        (aligned_left - 0.067, c_box.y1 + 0.038, "B"),
        (b_box.x0 - 0.067, b_box.y1 + 0.038, "C"),
    ]
    for x, y, letter in letter_specs:
        fig.text(x, y, letter, ha="left", va="top", fontsize=20,
                 fontweight="bold", color=TEXT_COLOR)
    for fmt in FORMATS:
        fig.savefig(OUT / f"{STEM}_wt_vs_cerk1_1.{fmt}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    gfp = load_gfp()
    elisa = load_elisa()
    gfp_summary, gfp_test = analyse(gfp, "normalised_gfp", "GFP")
    elisa_summary, elisa_test = analyse(elisa, "a450", "COVA ELISA")
    tests = pd.concat([gfp_test, elisa_test], ignore_index=True)
    save_logs(gfp, elisa, gfp_summary, elisa_summary, tests)
    make_figure(gfp, elisa, gfp_summary, elisa_summary, tests)
    print(gfp_summary.to_string(index=False))
    print(elisa_summary.to_string(index=False))
    print(tests[["assay", "t", "df", "p"]].to_string(index=False))
