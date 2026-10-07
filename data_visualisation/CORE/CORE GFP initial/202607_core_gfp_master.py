# %% [markdown]
# # CORE GFP master figure
#
# A: combined endpoint GFP with representative signal spots.
# B: no-p19/p19 CORE 1-2 endpoint GFP.
# C: CORE 1-2 western and Ponceau loading control.
# D: representative Cas9 WT and core-1-2 fluorescence images.

# %%
from pathlib import Path
import importlib.util
import os

MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd() / "CORE" / "CORE GFP initial"
SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_core_gfp_master"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
SAVE_DIR.mkdir(exist_ok=True)
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
LOG_DIR.mkdir(exist_ok=True)
SAVE_FMTS = ["png", "pdf", "svg"]

AUC_VALUES = [14584.208, 29444.392, 22358.179, 6481.974, 24634.108, 26614.472]


def load_combined_module():
    path = HERE / "combinedGFP.py"
    spec = importlib.util.spec_from_file_location("combined_gfp", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


combined = load_combined_module()


def load_initial_module():
    path = HERE / "202607_core_gfp_initial.py"
    spec = importlib.util.spec_from_file_location("initial_gfp", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


initial = load_initial_module()
western_ponceau_full = mpimg.imread(HERE / "core-1-2_westernandponceau_RBCL_original_panels.png")
# Approved photographic interiors from the formatted composite.
western_panel = western_ponceau_full[226:535, 118:1151]
ponceau_panel = western_ponceau_full[667:898, 118:1151]
cas9_image = mpimg.imread(HERE / "cas9WT.png")
core_image = mpimg.imread(HERE / "CORE-1-2.png")


def panel_label(fig, ax, label, dx=-0.020, dy=0.012):
    bbox = ax.get_position()
    fig.text(
        bbox.x0 + dx,
        bbox.y1 + dy,
        label,
        ha="left",
        va="top",
        fontsize=24,
        fontweight="bold",
        color="black",
    )


def aligned_panel_label(fig, ax, label, x, dy=0.012):
    bbox = ax.get_position()
    fig.text(
        x,
        bbox.y1 + dy,
        label,
        ha="left",
        va="top",
        fontsize=24,
        fontweight="bold",
        color="black",
    )


def draw_panel_a(fig, subspec, data, stats_table):
    grid = subspec.subgridspec(2, 1, height_ratios=[4.7, 1.12], hspace=0.19)
    plot_ax = fig.add_subplot(grid[0])
    spot_ax = fig.add_subplot(grid[1])

    combined.draw_bars(plot_ax, data, stats_table, show_xlabels=True, font_increase=4)
    plot_ax.yaxis.label.set_fontsize(20)
    plot_ax.tick_params(axis="y", labelsize=20)
    plot_ax.tick_params(axis="x", labelsize=20, pad=7)
    for text in plot_ax.texts:
        text.set_fontsize(20)
    plot_ax.set_xlabel("")

    spots = combined.load_signal_spots(combined.SPOT_DIR)
    spot_ax.set_facecolor("black")
    spot_ax.set_xlim(-0.5, len(combined.SAMPLE_ORDER) - 0.5)
    spot_ax.set_ylim(0, 1)
    spot_ax.set_aspect("equal", adjustable="box")
    spot_ax.set_xticks([])
    spot_ax.set_yticks([])
    for spine in spot_ax.spines.values():
        spine.set_visible(False)

    half_width = 0.40
    for x, sample in enumerate(combined.SAMPLE_ORDER):
        spot_ax.imshow(
            spots[sample],
            extent=(x - half_width, x + half_width, 0.10, 0.90),
            interpolation="lanczos",
            aspect="equal",
        )
    return plot_ax


def draw_panel_b(ax, data, stats_table):
    positions = {
        ("batch1", "no p19", "cas9 WT"): 0.00,
        ("batch1", "no p19", "core 1-2"): 0.58,
        ("batch1", "p19", "cas9 WT"): 1.42,
        ("batch1", "p19", "core 1-2"): 2.00,
    }
    xticks = [0.00, 0.58, 1.42, 2.00]
    xlabels = ["Cas9 WT", "core-1-2", "Cas9 WT", "core-1-2"]
    condition_centers = [0.29, 1.71]
    tops = {}
    ymax = 0.0
    for condition in initial.DISPLAY_COMPARISONS:
        for sample in initial.SAMPLE_ORDER:
            values = data.loc[
                data["batch"].eq(condition["batch"])
                & data["p19"].eq(condition["p19"])
                & data["sample_type"].eq(sample),
                "normalized_gfp",
            ].dropna() / initial.Y_SCALE
            x = positions[(condition["batch"], condition["p19"], sample)]
            mean = values.mean()
            error = initial.sem(values)
            fill = combined.WT_COLOR if sample == "cas9 WT" else combined.BAR_COLOR
            ax.bar(x, mean, width=initial.BAR_WIDTH, color=fill, alpha=initial.BAR_ALPHA,
                   edgecolor="black", linewidth=1.0, zorder=2)
            ax.errorbar(x, mean, yerr=error, fmt="none", ecolor=initial.ERRORBAR_COLOR,
                        elinewidth=1.1, capsize=3, capthick=1.1, zorder=5)
            jitter = np.linspace(-0.045, 0.045, len(values)) if len(values) > 1 else np.array([0.0])
            ax.scatter(x + jitter, values, s=22, facecolor=fill, edgecolor="none",
                       alpha=initial.INDIVIDUAL_POINT_ALPHA, zorder=4)
            label_y = max(mean + error, values.max()) + 0.65
            ax.text(x, label_y, str(len(values)), ha="center", va="bottom",
                    fontsize=14, color="#8A8A8A")
            tops[(condition["batch"], condition["p19"])] = max(
                tops.get((condition["batch"], condition["p19"]), 0), label_y
            )
            ymax = max(ymax, label_y + 0.7)

    for condition in initial.DISPLAY_COMPARISONS:
        row = stats_table.loc[
            stats_table["batch"].eq(condition["batch"])
            & stats_table["p19"].eq(condition["p19"])
        ]
        if row.empty:
            continue
        x1 = positions[(condition["batch"], condition["p19"], "cas9 WT")]
        x2 = positions[(condition["batch"], condition["p19"], "core 1-2")]
        y = tops[(condition["batch"], condition["p19"])] + 4.00
        ax.plot([x1, x2], [y, y], color=initial.TEXT_COLOR, linewidth=0.9, clip_on=False)
        ax.text((x1 + x2) / 2, y + 0.35, row["annotation"].iloc[0], ha="center",
                va="bottom", fontsize=15, color=initial.TEXT_COLOR)
        ymax = max(ymax, y + 1.0)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=20)
    for label, sample_label in zip(ax.get_xticklabels(), xlabels):
        if sample_label.startswith("core"):
            label.set_fontstyle("italic")
    for center, condition in zip(condition_centers, initial.DISPLAY_COMPARISONS):
        ax.text(center, -0.17, condition["label"], transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=20, color=initial.TEXT_COLOR)
    ax.set_ylabel("Relative GFP\n(AU x1000)", fontsize=20, labelpad=3)
    ax.tick_params(axis="y", labelsize=20)
    for text in ax.texts:
        text.set_fontsize(20)
    ax.set_ylim(0, max(30, np.ceil(ymax / 10) * 10))
    initial.style_axis(ax)


def draw_image_panel(ax, image):
    ax.imshow(image, interpolation="lanczos")
    ax.set_axis_off()
    ax.set_anchor("C")


def auc_tables():
    raw = pd.DataFrame(
        {
            "lane": np.arange(1, 7),
            "genotype": ["Cas9 WT"] * 3 + ["core-1-2"] * 3,
            "biological_replicate": [1, 2, 3, 1, 2, 3],
            "auc": AUC_VALUES,
            "measurement": "integrated alpha-GFP band density",
            "source_panel": "core-1-2 western, lanes shown left-to-right",
        }
    )
    summary = (
        raw.groupby("genotype", sort=False)
        .agg(
            n=("auc", "size"),
            mean_auc=("auc", "mean"),
            sem_auc=("auc", "sem"),
            sd_auc=("auc", "std"),
        )
        .reset_index()
    )
    wt = raw.loc[raw["genotype"].eq("Cas9 WT"), "auc"]
    core = raw.loc[raw["genotype"].eq("core-1-2"), "auc"]
    result = stats.ttest_ind(wt, core, equal_var=False, alternative="two-sided")
    test = pd.DataFrame(
        [{
            "comparison": "Cas9 WT vs core-1-2",
            "test": "two-sided Welch t-test",
            "cas9_wt_n": len(wt),
            "core_1_2_n": len(core),
            "cas9_wt_mean_auc": wt.mean(),
            "core_1_2_mean_auc": core.mean(),
            "t_statistic": result.statistic,
            "p_value": result.pvalue,
            "annotation": "ns" if result.pvalue >= 0.05 else "*",
        }]
    )
    return raw, summary, test


def write_auc_workbook(raw, summary, test):
    path = LOG_DIR / "core_1_2_western_auc_quantification.xlsx"
    with pd.ExcelWriter(path) as writer:
        raw.to_excel(writer, sheet_name="lane_auc_values", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
        test.to_excel(writer, sheet_name="welch_test", index=False)
    return path


def draw_auc_plot(ax, raw, summary, test):
    order = ["Cas9 WT", "core-1-2"]
    colors = [combined.WT_COLOR, combined.BAR_COLOR]
    x = np.arange(2)
    means = [summary.loc[summary["genotype"].eq(group), "mean_auc"].iloc[0] / 1000 for group in order]
    errors = [summary.loc[summary["genotype"].eq(group), "sem_auc"].iloc[0] / 1000 for group in order]
    for index, (group, color) in enumerate(zip(order, colors)):
        values = raw.loc[raw["genotype"].eq(group), "auc"].to_numpy() / 1000
        ax.bar(index, means[index], width=0.58, color=color, alpha=0.74,
               edgecolor="black", linewidth=1.0, zorder=2)
        ax.errorbar(index, means[index], yerr=errors[index], fmt="none", color="black",
                    elinewidth=1.1, capsize=3, capthick=1.1, zorder=4)
        ax.scatter(index + np.linspace(-0.055, 0.055, len(values)), values,
                   s=28, color=color, alpha=0.72, edgecolor="none", zorder=3)
    y = max(np.array(means) + np.array(errors)) * 1.22
    ax.plot([0, 1], [y, y], color=combined.TEXT_COLOR, linewidth=0.9)
    ax.text(0.5, y * 1.025, test["annotation"].iloc[0], ha="center", va="bottom", fontsize=20)
    ax.set_ylim(0, y * 1.20)
    ax.set_xticks(x)
    ax.set_xticklabels(["Cas9\nWT", "core-\n1-2"], fontsize=20, rotation=0, va="top")
    ax.get_xticklabels()[1].set_fontstyle("italic")
    ax.set_ylabel("Integrated α-GFP band density\n(AU x1000)", fontsize=20, labelpad=4)
    ax.tick_params(axis="y", labelsize=20)
    initial.style_axis(ax)


def draw_panel_c(fig, subspec, auc_raw, auc_summary, auc_test):
    grid = subspec.subgridspec(1, 2, width_ratios=[1.88, 1.00], wspace=0.34)
    image_ax = fig.add_subplot(grid[0, 0])
    auc_container = fig.add_subplot(grid[0, 1])
    auc_container.set_axis_off()
    auc_ax = auc_container.inset_axes([0.08, 0.16, 0.82, 0.68])
    image_ax.set_axis_off()
    western_ax = image_ax.inset_axes([-0.020, 0.460, 0.72, 0.235])
    ponceau_ax = image_ax.inset_axes([-0.020, 0.105, 0.72, 0.245])
    for inset, image in [(western_ax, western_panel), (ponceau_ax, ponceau_panel)]:
        inset.imshow(image, interpolation="lanczos", aspect="auto")
        inset.set_xticks([])
        inset.set_yticks([])
        for spine in inset.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.0)

    image_ax.text(0.160, 0.860, "Cas9 WT", transform=image_ax.transAxes, ha="center", va="center",
                  fontsize=20, fontfamily=combined.FONT_STACK)
    image_ax.text(0.520, 0.860, "core-1-2", transform=image_ax.transAxes, ha="center", va="center",
                  fontsize=20, fontstyle="italic", fontfamily=combined.FONT_STACK)
    image_ax.plot([0.005, 0.315], [0.810, 0.810], transform=image_ax.transAxes,
                  color="black", linewidth=1.0)
    image_ax.plot([0.355, 0.665], [0.810, 0.810], transform=image_ax.transAxes,
                  color="black", linewidth=1.0)
    for x, label in zip([0.035, 0.160, 0.285, 0.395, 0.520, 0.645], ["1", "2", "3", "1", "2", "3"]):
        image_ax.text(x, 0.755, label, transform=image_ax.transAxes, ha="center", va="center",
                      fontsize=20, fontfamily=combined.FONT_STACK)
    image_ax.text(0.940, 0.580, "~25 kDa", transform=image_ax.transAxes, ha="right", va="center",
                  fontsize=20, fontfamily=combined.FONT_STACK)
    image_ax.text(0.700, 0.405, "α-GFP", transform=image_ax.transAxes, ha="right", va="center",
                  fontsize=20, fontstyle="italic", fontfamily=combined.FONT_STACK)
    image_ax.plot(0.755, 0.227, marker="<", markersize=8, color="black",
                  transform=image_ax.transAxes, clip_on=False)
    image_ax.text(0.940, 0.227, "RBCL", transform=image_ax.transAxes, ha="right", va="center",
                  fontsize=20, fontfamily=combined.FONT_STACK)
    image_ax.text(0.700, 0.050, "Ponceau S", transform=image_ax.transAxes, ha="right", va="center",
                  fontsize=20, fontstyle="italic", fontfamily=combined.FONT_STACK)
    draw_auc_plot(auc_ax, auc_raw, auc_summary, auc_test)
    return image_ax


def draw_panel_d(fig, subspec):
    image_height = min(cas9_image.shape[0], core_image.shape[0])
    images = [cas9_image[:image_height], core_image[:image_height]]
    combined_image = np.concatenate(images, axis=1)
    ax = fig.add_subplot(subspec)
    labels = ["Cas9 WT", "core-1-2"]
    ax.imshow(combined_image, interpolation="lanczos")
    ax.set_axis_off()
    ax.set_anchor("C")
    for x, label in zip([0.25, 0.75], labels):
        ax.text(
            x,
            0.965,
            label,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=20,
            fontweight="normal",
            fontfamily=combined.FONT_STACK,
            fontstyle="italic" if label.startswith("core") else "normal",
            color="white",
        )
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                           facecolor="none", edgecolor="white", linewidth=2.0,
                           zorder=6))
    return ax


def save(fig):
    for file_format in SAVE_FMTS:
        fig.savefig(
            SAVE_DIR / f"{SCRIPT_STEM}_core_gfp_master.{file_format}",
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )


def main():
    raw = combined.load_raw_data(combined.SRC)
    data = combined.analysis_data(raw)
    stats_table = combined.run_stats(data)
    initial_raw = initial.load_raw_data(initial.SRC)
    initial_data = initial.analysis_data(initial_raw)
    initial_stats = initial.run_stats(initial_data)
    auc_raw, auc_summary, auc_test = auc_tables()
    auc_path = write_auc_workbook(auc_raw, auc_summary, auc_test)

    plt.rcParams.update(combined.CORE_RCPARAMS if hasattr(combined, "CORE_RCPARAMS") else {})
    fig = plt.figure(figsize=(20.8, 15.2), facecolor="white")
    outer = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.30, 1.12],
        wspace=0.17,
    )

    panel_a_ax = draw_panel_a(fig, outer[0, 0], data, stats_table)
    right = outer[0, 1].subgridspec(
        3, 1, height_ratios=[0.80, 1.55, 1.00], hspace=0.18
    )
    panel_b_ax = fig.add_subplot(right[0, 0])
    draw_panel_b(panel_b_ax, initial_data, initial_stats)
    panel_c_ax = draw_panel_c(fig, right[1, 0], auc_raw, auc_summary, auc_test)
    panel_d_ax = draw_panel_d(fig, right[2, 0])

    fig.subplots_adjust(left=0.075, right=0.950, top=0.955, bottom=0.055)

    panel_label(fig, panel_a_ax, "A", dx=-0.032, dy=0.020)
    right_label_x = panel_b_ax.get_position().x0 - 0.045
    aligned_panel_label(fig, panel_b_ax, "B", right_label_x)
    aligned_panel_label(fig, panel_c_ax, "C", right_label_x)
    aligned_panel_label(fig, panel_d_ax, "D", right_label_x)
    save(fig)
    plt.close(fig)
    print(f"Master figures saved to: {SAVE_DIR}")
    print(f"AUC workbook saved to: {auc_path}")


if __name__ == "__main__":
    main()
