"""Create the NPR1 pre-infiltration GFP plot and its western-blot composite."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps
from scipy import stats


HERE = Path(__file__).resolve().parent
DATA = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/6. PosHits/NPR1/NPR1_Preinfiltration_Check/Pre-infiltration_results.xlsx"
)
PANEL_B = HERE / "preinfil_npr1.png"
PRIMING_GRAPHIC = HERE / "priming npr1.png"
OUTPUT = HERE / "figures"
PANEL_B_RELABELED = OUTPUT / "202607_npr1_preinfiltration_panel_B_relabelled.png"

TEXT_COLOR = "#4A4A4A"
GRID_COLOR = "#FFFFFF"
BAR_OUTLINE = "#6F7780"
COLORS = {"No pre-infil": "#A6B0BA", "pre-infil": "#536FA8"}
POINT_COLORS = {"No pre-infil": "#65727E", "pre-infil": "#315F78"}
GENOTYPES = ["WT", "NPR1-"]
TREATMENTS = ["No pre-infil", "pre-infil"]


def load_data():
    data = pd.read_excel(DATA, sheet_name="Raw", skiprows=4)
    data = data.rename(
        columns={
            "Plant Type": "genotype",
            "Treatment Type": "treatment",
            "Final GFP Intensity": "intensity",
        }
    )
    data = data.loc[:, ["genotype", "treatment", "intensity"]].dropna()
    data["genotype"] = data["genotype"].replace({"NPR1-": "NPR1-"})
    data["intensity"] = pd.to_numeric(data["intensity"], errors="coerce")
    return data.dropna(subset=["intensity"])


def run_stats(data):
    results = []
    for genotype in GENOTYPES:
        no_pre = data.loc[
            (data["genotype"] == genotype) & (data["treatment"] == "No pre-infil"),
            "intensity",
        ].to_numpy()
        pre = data.loc[
            (data["genotype"] == genotype) & (data["treatment"] == "pre-infil"),
            "intensity",
        ].to_numpy()
        result = stats.ttest_ind(no_pre, pre, equal_var=False, alternative="two-sided")
        results.append({
            "genotype": genotype,
            "test": "two-sided Welch t-test",
            "no_pre_infil_n": len(no_pre),
            "pre_infil_n": len(pre),
            "t_statistic": result.statistic,
            "df": result.df,
            "p_raw": result.pvalue,
        })

    table = pd.DataFrame(results)
    order = table["p_raw"].argsort().to_numpy()
    adjusted = np.empty(len(table), dtype=float)
    running_max = 0.0
    for rank, index in enumerate(order):
        running_max = max(running_max, (len(table) - rank) * table.loc[index, "p_raw"])
        adjusted[index] = min(running_max, 1.0)
    table["p_holm"] = adjusted
    table["annotation"] = table["p_holm"].map(p_value_label)
    return table


def p_value_label(p_value):
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def plot_panel_a(data, output_stem):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 19,
        "axes.titlesize": 21,
        "axes.labelsize": 20,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 17,
        "text.color": TEXT_COLOR,
        "axes.labelcolor": TEXT_COLOR,
        "axes.edgecolor": TEXT_COLOR,
        "xtick.color": TEXT_COLOR,
        "ytick.color": TEXT_COLOR,
        "savefig.dpi": 300,
    })

    centers = np.arange(len(GENOTYPES), dtype=float)
    offsets = np.linspace(-0.19, 0.19, len(TREATMENTS))
    width = 0.32
    fig, ax = plt.subplots(figsize=(7.5, 5.4), facecolor="white")
    ax.set_facecolor("#E9EFF6")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID_COLOR, linewidth=1.0)

    for genotype_index, genotype in enumerate(GENOTYPES):
        for treatment_index, treatment in enumerate(TREATMENTS):
            position = centers[genotype_index] + offsets[treatment_index]
            values = data.loc[
                (data["genotype"] == genotype) & (data["treatment"] == treatment), "intensity"
            ].to_numpy()
            mean = values.mean()
            sem = values.std(ddof=1) / np.sqrt(len(values))
            ax.bar(
                position,
                mean / 1000,
                width=width,
                color=COLORS[treatment],
                alpha=0.84,
                edgecolor=BAR_OUTLINE,
                linewidth=0.8,
                zorder=2,
            )
            ax.errorbar(
                position,
                mean / 1000,
                yerr=sem / 1000,
                fmt="none",
                ecolor=TEXT_COLOR,
                elinewidth=1.2,
                capsize=3,
                capthick=1.2,
                zorder=4,
            )
            jitter = np.linspace(-0.06, 0.06, len(values))
            ax.scatter(
                np.full(len(values), position) + jitter,
                values / 1000,
                s=52,
                color=POINT_COLORS[treatment],
                edgecolors="white",
                linewidths=0.6,
                alpha=0.9,
                zorder=5,
            )

    ax.set_xticks(centers, ["WT", "npr1#1"])
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_ylabel("GFP Intensity (x1000)")
    ax.set_ylim(0, 52)
    ax.set_yticks(np.arange(0, 51, 5))
    ax.set_xlim(-0.6, 1.6)
    ax.tick_params(axis="x", length=0, pad=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(TEXT_COLOR)
    ax.spines["bottom"].set_color(TEXT_COLOR)
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLORS[treatment])
        for treatment in reversed(TREATMENTS)
    ]
    ax.legend(
        legend_handles,
        ["pre-infil", "No pre-infil"],
        frameon=False,
        ncol=2,
        mode="expand",
        loc="lower center",
        bbox_to_anchor=(0.5, 1.015),
        borderaxespad=0.0,
        columnspacing=1.0,
        handletextpad=0.55,
    )

    stats_table = run_stats(data)
    for genotype_index, genotype in enumerate(GENOTYPES):
        row = stats_table.loc[stats_table["genotype"] == genotype].iloc[0]
        pair_values = data.loc[data["genotype"] == genotype, "intensity"] / 1000
        pair_errors = []
        for treatment in TREATMENTS:
            values = data.loc[
                (data["genotype"] == genotype) & (data["treatment"] == treatment),
                "intensity",
            ]
            pair_errors.append(values.mean() / 1000 + values.sem() / 1000)
        bracket_y = max(pair_errors) + 2.0
        left = centers[genotype_index] + offsets[0]
        right = centers[genotype_index] + offsets[1]
        bracket_height = 0.55
        ax.plot(
            [left, left, right, right],
            [bracket_y, bracket_y + bracket_height, bracket_y + bracket_height, bracket_y],
            color=TEXT_COLOR,
            linewidth=1.1,
            clip_on=False,
        )
        ax.text(
            (left + right) / 2,
            bracket_y + bracket_height + 0.25,
            row["annotation"],
            ha="center",
            va="bottom",
            fontsize=17,
            color=TEXT_COLOR,
        )
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.15, top=0.82)
    fig.savefig(output_stem.with_suffix(".png"), bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def open_rgb(path):
    image = Image.open(path)
    if image.mode == "RGBA":
        background = Image.new("RGBA", image.size, "white")
        image = Image.alpha_composite(background, image)
    return image.convert("RGB")


def trim_white(image, threshold=249, padding=12):
    mask = Image.new("L", image.size)
    mask.putdata([255 if min(pixel) < threshold else 0 for pixel in image.getdata()])
    bbox = mask.getbbox()
    if bbox is None:
        return image
    left, top, right, bottom = bbox
    return image.crop(
        (max(0, left - padding), max(0, top - padding),
         min(image.width, right + padding), min(image.height, bottom + padding))
    )


def redraw_panel_b_labels(output_path):
    image = open_rgb(PANEL_B)
    draw = ImageDraw.Draw(image)
    font_path = "/System/Library/Fonts/HelveticaNeue.ttc"
    regular = ImageFont.truetype(font_path, 108, index=0)
    italic = ImageFont.truetype(font_path, 108, index=2)

    # Clear only the existing text in the white margins; the blot, borders,
    # lane dividers, and underlines remain from the original panel.
    for box in (
        (250, 85, 570, 205), (790, 85, 1180, 205),
        (220, 260, 1170, 370), (990, 750, 1260, 860),
        (990, 1220, 1240, 1340),
    ):
        draw.rectangle(box, fill="white")

    def centered(text, center_x, baseline_y, font):
        draw.text(
            (center_x, baseline_y),
            text,
            font=font,
            fill=TEXT_COLOR,
            anchor="mm",
        )

    centered("WT", 430, 150, regular)
    centered("npr1#1", 995, 150, italic)
    for text, x in zip(("1", "2", "1", "2"), (300, 550, 860, 1110)):
        centered(text, x, 315, regular)
    centered("Actin", 1140, 805, italic)
    centered("PR1", 1120, 1280, italic)
    image.save(output_path, dpi=(300, 300), optimize=True)


def make_composite(plot_path, output_stem):
    priming = trim_white(open_rgb(PRIMING_GRAPHIC), threshold=252, padding=8)
    plot = trim_white(open_rgb(plot_path))
    blot = trim_white(open_rgb(PANEL_B_RELABELED), threshold=252, padding=8)
    canvas = Image.new("RGB", (2700, 2920), "white")
    draw = ImageDraw.Draw(canvas)
    label_font = ImageFont.truetype("/System/Library/Fonts/HelveticaNeue.ttc", 82, index=0)

    for image, box in (
        (priming, (150, 90, 2550, 1380)),
        (plot, (70, 1510, 1570, 2880)),
        (blot, (1760, 1790, 2630, 2490)),
    ):
        left, top, right, bottom = box
        fitted = ImageOps.contain(image, (right - left, bottom - top), Image.Resampling.LANCZOS)
        x = left + (right - left - fitted.width) // 2
        y = top + (bottom - top - fitted.height) // 2
        canvas.paste(fitted, (x, y))

    draw.text((34, 34), "A", font=label_font, fill=TEXT_COLOR)
    draw.text((34, 1450), "B", font=label_font, fill=TEXT_COLOR)
    draw.text((1780, 1450), "C", font=label_font, fill=TEXT_COLOR)
    canvas.save(output_stem.with_suffix(".png"), dpi=(300, 300), optimize=True)
    canvas.save(output_stem.with_suffix(".pdf"), resolution=300)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = load_data()
    panel_a = OUTPUT / "202607_npr1_preinfiltration_panel_A"
    composite = OUTPUT / "202607_npr1_preinfiltration_composite"
    stats_table = run_stats(data)
    stats_table.to_csv(OUTPUT / "202607_npr1_preinfiltration_stats.csv", index=False)
    plot_panel_a(data, panel_a)
    redraw_panel_b_labels(PANEL_B_RELABELED)
    make_composite(panel_a.with_suffix(".png"), composite)
    print(f"Saved figures to {OUTPUT}")


if __name__ == "__main__":
    main()
