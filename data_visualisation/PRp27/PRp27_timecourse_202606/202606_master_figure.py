# %% [markdown]
# # Combined GFP/RNA/ELISA/ColonyCount master figure
#
# This script composes selected panels from the analysis scripts without the
# long methods paragraphs. It reuses the plotting helpers from the individual
# scripts so rerunning this file picks up future changes to those plot functions.
# Data are read from the exported `*_logs` workbooks; rerun the individual
# analysis scripts first if source data, cleaning, or statistics change.

# %%
from pathlib import Path
import ast
import os

os.environ.setdefault("MPLBACKEND", "Agg")
MPLCONFIGDIR = Path.cwd() / ".matplotlib-cache"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

LOCAL_CACHE = Path.cwd() / ".cache"
(LOCAL_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(LOCAL_CACHE))

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from figure_styles import (
    COLORS,
    LINESTYLES,
    MARKER_FACES,
    MARKERS,
    MASTER_RCPARAMS,
    TEXT_COLOR,
    apply_axis_style,
    legend_handles,
    make_errorbars_black,
    panel_label,
    recolor_grouped_scatter_by_x,
)


def plant_line_legend(fig, **kwargs):
    """Create a legend while retaining the chapter font for italic line names."""
    legend = fig.legend(**kwargs)
    for legend_text in legend.get_texts():
        if legend_text.get_text() in {"prp27-1", "prp27-2"}:
            legend_text.set_fontstyle("italic")
    return legend

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
elif (Path.cwd() / "PRp27" / "PRp27_timecourse_202606").is_dir():
    HERE = Path.cwd() / "PRp27" / "PRp27_timecourse_202606"
elif (Path.cwd() / "PRp27_timecourse_202606").is_dir():
    HERE = Path.cwd() / "PRp27_timecourse_202606"
else:
    HERE = Path.cwd()

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202606_master_figure"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
SAVE_DIR.mkdir(exist_ok=True)
SPLIT_SAVE_DIR = HERE / "202606_split_figures"
SPLIT_SAVE_DIR.mkdir(exist_ok=True)
SAVE_FMTS = ["png", "pdf", "svg"]


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(
            SAVE_DIR / f"{output_stem}.{file_format}",
            dpi=300,
            bbox_inches="tight",
        )


def save_split(fig, name):
    """Save the smaller thesis-ready figure sets in their dedicated folder."""
    for file_format in SAVE_FMTS:
        fig.savefig(
            SPLIT_SAVE_DIR / f"202606_{name}.{file_format}",
            dpi=300,
            bbox_inches="tight",
        )


SHARED_CONSTANTS = {
    "COLORS",
    "MARKERS",
    "LINESTYLES",
    "MARKER_FACES",
    "ERRBAR",
    "SIG_ALPHA",
    "YLABEL",
    "XLABEL",
    "LINES",
    "DPIS",
    "MAIN_DPIS",
    "DILS",
    "P19_ORDER",
    "FONT_FAMILY",
}


def load_plot_namespace(script_name):
    """Load imports, constants, and functions without running plotting cells."""
    script_path = HERE / script_name
    tree = ast.parse(script_path.read_text(), filename=str(script_path))
    namespace = {"__file__": str(script_path)}

    for node in tree.body:
        should_exec = isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))
        if isinstance(node, ast.Assign):
            target_names = [target.id for target in node.targets if isinstance(target, ast.Name)]
            should_exec = any(name in SHARED_CONSTANTS for name in target_names)
        if should_exec:
            exec(compile(ast.Module([node], type_ignores=[]), str(script_path), "exec"), namespace)

    return namespace


def apply_shared_plot_constants(namespace):
    """Use master style constants in plotting helpers loaded from scripts."""
    namespace["COLORS"] = COLORS
    namespace["MARKERS"] = MARKERS
    if "LINESTYLES" in namespace:
        namespace["LINESTYLES"] = LINESTYLES
    if "MARKER_FACES" in namespace:
        namespace["MARKER_FACES"] = MARKER_FACES


gfp = load_plot_namespace("202606_GFP.py")
rna = load_plot_namespace("202606_RNA.py")
elisa = load_plot_namespace("202606_ELISA.py")
virb = load_plot_namespace("202606_virB.py")
colony = load_plot_namespace("202606_ColonyCount.py")

for plot_namespace in (gfp, rna, elisa, virb, colony):
    apply_shared_plot_constants(plot_namespace)

gfp_clean = pd.read_excel(HERE / "202606_GFP_logs" / "gfp_cleaned_long.xlsx")
gfp_clean["line"] = pd.Categorical(gfp_clean["line"], categories=gfp["LINES"])
gfp_clean["p19"] = pd.Categorical(gfp_clean["p19"], categories=gfp["P19_ORDER"])
gfp_rm = gfp["repeat_means"](gfp_clean)
gfp_point_contrasts = pd.read_excel(
    HERE / "202606_GFP_logs" / "gfp_pointwise_sample_type_pairwise_contrasts.xlsx"
)

rna_clean = pd.read_excel(HERE / "202606_RNA_logs" / "rna_cleaned_long.xlsx")
rna_clean["line"] = pd.Categorical(rna_clean["line"], categories=rna["LINES"])
rna_point_contrasts = pd.read_excel(
    HERE / "202606_RNA_logs" / "rna_pointwise_sample_type_pairwise_contrasts.xlsx"
)

elisa_auc = pd.read_excel(HERE / "202606_ELISA_logs" / "elisa_auc_values.xlsx")
elisa_auc_contrasts = pd.read_excel(
    HERE / "202606_ELISA_logs" / "elisa_auc_pointwise_sample_type_pairwise_contrasts.xlsx"
)

virb_clean = pd.read_excel(HERE / "202606_virB_logs" / "virb_cleaned_long.xlsx")
virb_clean["line"] = pd.Categorical(virb_clean["line"], categories=virb["LINES"])
virb_clean["p19"] = pd.Categorical(virb_clean["p19"], categories=virb["P19_ORDER"])
virb_rm = virb["repeat_means"](virb_clean)
virb_point_contrasts = pd.read_excel(
    HERE / "202606_virB_logs" / "virb_pointwise_sample_type_pairwise_contrasts.xlsx"
)

colony_clean = pd.read_excel(
    HERE / "202606_ColonyCount_logs" / "colony_cleaned_biological_means.xlsx"
)
colony_clean["line"] = pd.Categorical(colony_clean["line"], categories=colony["LINES"])
colony_point_contrasts = pd.read_excel(
    HERE / "202606_ColonyCount_logs" / "colony_pointwise_sample_type_pairwise_contrasts.xlsx"
)
plant_photo = mpimg.imread(HERE / "Images" / "prp27PlantPhotos_relabelled_10cm_scale.png")
rna_formatted = mpimg.imread(HERE / "Images" / "RNA_formatted_cropped.png")
gfp_western = mpimg.imread(HERE / "Images" / "GFP_western_white_cropped.png")
virb_scan = mpimg.imread(HERE / "Images" / "P2first_virB.png")

plt.rcParams.update(MASTER_RCPARAMS)


def annotate_auc_letters(ax, contrasts, auc_data, force_show_dpis=()):
    """Show compact letters for AUC dpi families only when contrasts are non-ns."""
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1
    label_gap_y = 0.10 * y_range
    label_step_y = 0.14 * y_range
    max_label_y = y_max
    x = range(len(elisa["DPIS"]))

    for dpi_index, dpi in enumerate(elisa["DPIS"]):
        rows = rna["contrast_subset"](contrasts, dpi=dpi)
        letters, should_show = rna["compact_letter_display"](rows)
        if not should_show and dpi not in force_show_dpis:
            continue
        if not letters and dpi in force_show_dpis:
            letters = {line: "a" for line in elisa["LINES"]}

        line_positions = []
        for line in elisa["LINES"]:
            values = auc_data[(auc_data["Sample Type"] == line) & (auc_data["dpi"] == dpi)]["auc"]
            if values.empty:
                continue
            mean = values.mean()
            error = values.sem() if elisa["ERRBAR"] == "SEM" else values.std(ddof=1)
            error_extent = mean + (0 if pd.isna(error) else error)
            point_extent = values.max()
            data_extent = error_extent if pd.isna(point_extent) else max(point_extent, error_extent)
            line_positions.append({"line": line, "mean": mean, "data_extent": data_extent})

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
                fontsize=14,
                fontweight="bold",
                ha="center",
                va="center",
                zorder=6,
            )
            max_label_y = max(max_label_y, y + label_step_y)

    if max_label_y > y_max:
        ax.set_ylim(y_min, max_label_y + label_step_y)


def image_panel(ax, image, interpolation="nearest"):
    ax.imshow(image, interpolation=interpolation)
    ax.set_axis_off()
    ax.set_anchor("C")


def labelled_plant_image_panel(ax, image, interpolation="nearest"):
    """Draw the plant photograph with thesis-consistent plant-line labels."""
    image_panel(ax, image, interpolation=interpolation)

    # Mask the labels baked into the source image, then replace them with the
    # same names and typography used by the figure legend.
    image_height, image_width = image.shape[:2]
    ax.add_patch(
        Rectangle(
            (0, 0),
            image_width,
            image_height * 0.15,
            facecolor="black",
            edgecolor="none",
            zorder=2,
        )
    )
    label_y = image_height * 0.09
    for x_fraction, label, fontstyle in (
        (0.155, "WT", "normal"),
        (0.500, "prp27-1", "italic"),
        (0.835, "prp27-2", "italic"),
    ):
        ax.text(
            image_width * x_fraction,
            label_y,
            label,
            color="white",
            fontsize=17,
            fontstyle=fontstyle,
            ha="center",
            va="center",
            zorder=3,
        )


def image_panel_without_white(ax, image, interpolation="nearest"):
    """Draw dark image content while allowing the panel colour to show through."""
    rgb = image[..., :3]
    if rgb.dtype.kind in "ui":
        rgb = rgb.astype(float) / np.iinfo(image.dtype).max
    luminance = rgb.mean(axis=2)
    # Values close to white belong to the exported image canvas rather than the
    # blot itself. Fade that canvas out while retaining grey bands and labels.
    alpha = np.clip((0.92 - luminance) / 0.20, 0, 1)
    rgba = np.dstack([rgb, alpha])
    ax.imshow(rgba, interpolation=interpolation)
    ax.set_axis_off()
    ax.set_anchor("C")


def image_panel_with_opaque_box(
    ax,
    image,
    box=(197, 183, 1326, 585),
    interpolation="nearest",
):
    """Remove the outer white canvas while preserving the blot box verbatim."""
    rgb = image[..., :3]
    if rgb.dtype.kind in "ui":
        rgb = rgb.astype(float) / np.iinfo(image.dtype).max
    luminance = rgb.mean(axis=2)
    left, top, right, bottom = box

    # Layer 1 retains dark text/rules outside the blot but removes its exported
    # white canvas. Keep the blot area out of this thresholded layer entirely.
    outer_alpha = np.clip((0.92 - luminance) / 0.20, 0, 1)
    outer_alpha[top : bottom + 1, left : right + 1] = 0
    ax.imshow(np.dstack([rgb, outer_alpha]), interpolation=interpolation)

    # Layer 2 restores the original blot pixels at full opacity, including the
    # original white-grey membrane background.
    blot_alpha = np.zeros_like(luminance)
    blot_alpha[top : bottom + 1, left : right + 1] = 1
    ax.imshow(np.dstack([rgb, blot_alpha]), interpolation=interpolation)
    ax.set_axis_off()
    ax.set_anchor("C")


def fixed_panel_label(fig, ax, label, x=0.033, dy=0.012):
    bbox = ax.get_position()
    fig.text(
        x,
        bbox.y1 + dy,
        label,
        ha="left",
        va="top",
        fontsize=20,
        fontweight="bold",
        color=TEXT_COLOR,
    )


def fixed_panel_label_at_y(fig, label, y, x=0.033):
    fig.text(
        x,
        y,
        label,
        ha="left",
        va="top",
        fontsize=20,
        fontweight="bold",
        color=TEXT_COLOR,
    )


def style_bars(ax, alpha=0.60, edgecolor="#000000", linewidth=1.5):
    for patch in ax.patches:
        patch.set_alpha(alpha)
        patch.set_edgecolor(edgecolor)
        patch.set_linewidth(linewidth)


def rounded_axis_group_box(
    fig,
    axes,
    facecolor="#EAF1F5",
    alpha=1.0,
    pad=0.018,
    pad_left=None,
    pad_right=None,
    pad_bottom=None,
    pad_top=None,
):
    """Place one lightly shaded rounded box behind a related group of axes."""
    boxes = [ax.get_position() for ax in axes]
    left = min(box.x0 for box in boxes) - (pad if pad_left is None else pad_left)
    bottom = min(box.y0 for box in boxes) - (pad if pad_bottom is None else pad_bottom)
    right = max(box.x1 for box in boxes) + (pad if pad_right is None else pad_right)
    top = max(box.y1 for box in boxes) + (pad if pad_top is None else pad_top)
    patch = FancyBboxPatch(
        (left, bottom),
        right - left,
        top - bottom,
        boxstyle="round,pad=0.008,rounding_size=0.018",
        transform=fig.transFigure,
        facecolor=facecolor,
        alpha=alpha,
        edgecolor="none",
        zorder=-10,
        clip_on=False,
    )
    fig.patches.append(patch)
    return patch


fig = plt.figure(figsize=(20.8, 18.8))
outer = fig.add_gridspec(
    4,
    1,
    height_ratios=[1.82, 1.58, 1.34, 1.02],
    hspace=0.30,
)

top_grid = outer[0, 0].subgridspec(
    1,
    3,
    width_ratios=[2.10, 2.28, 1.16],
    wspace=0.11,
)
photo_ax = fig.add_subplot(top_grid[0, 0])
rna_image_ax = fig.add_subplot(top_grid[0, 1])
rna_ax = fig.add_subplot(top_grid[0, 2])

gfp_grid = outer[1, 0].subgridspec(1, 3, width_ratios=[0.78, 0.78, 2.20], wspace=0.12)
gfp_axes = [fig.add_subplot(gfp_grid[0, idx]) for idx in range(len(gfp["P19_ORDER"]))]
gfp_western_ax = fig.add_subplot(gfp_grid[0, 2])

virb_grid = outer[2, 0].subgridspec(1, 3, width_ratios=[1.28, 1.28, 0.50], wspace=0.10)
virb_axes = [fig.add_subplot(virb_grid[0, idx]) for idx in range(len(virb["P19_ORDER"]))]
virb_scan_ax = fig.add_subplot(virb_grid[0, 2])

bottom_grid = outer[3, 0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.34)
elisa_ax = fig.add_subplot(bottom_grid[0, 0])
colony_ax = fig.add_subplot(bottom_grid[0, 1])

# A. Plant phenotype photo
photo_ax.set_facecolor("black")
image_panel(photo_ax, plant_photo)
photo_ax.set_anchor("C")

# B. RNA formatted panel and RNA grouped plot
image_panel(rna_image_ax, rna_formatted)

rna["draw_grouped_bar_by_dpi"](rna_ax, rna_clean)
rna["annotate_grouped_bar_letters"](
    rna_ax,
    rna_point_contrasts,
    rna_clean,
    force_show_dpis=tuple(rna["DPIS"]),
)
rna_ax.set_ylabel(rna["YLABEL"])
style_bars(rna_ax)

# C. GFP no-p19/p19 timecourse and western blot
for ax, p19 in zip(gfp_axes, gfp["P19_ORDER"]):
    gfp["draw_timecourse"](ax, gfp_rm, p19, "batch_mean", point_mode="batch")
    gfp["annotate_timecourse_letters"](
        ax,
        gfp_point_contrasts,
        gfp_rm,
        p19,
        "batch_mean",
        point_mode="batch",
        force_show_dpis=(0, 1),
    )

gfp_axes[0].set_ylabel(gfp["YLABEL"])
image_panel(gfp_western_ax, gfp_western)

# E. ELISA AUC
elisa["draw_auc_by_dpi"](elisa_ax, elisa_auc)
annotate_auc_letters(
    elisa_ax,
    elisa_auc_contrasts,
    elisa_auc,
    force_show_dpis=tuple(elisa["DPIS"]),
)
recolor_grouped_scatter_by_x(elisa_ax, COLORS)
elisa_ax.set_ylabel(r"AUC of normalised A$_{450}$")
style_bars(elisa_ax)

# D. virB infiltration timecourse and scan
for ax, condition in zip(virb_axes, virb["P19_ORDER"]):
    virb["draw_timecourse"](ax, virb_rm, condition, "batch_mean", point_mode="batch")
    virb["annotate_timecourse_letters"](
        ax,
        virb_point_contrasts,
        virb_rm,
        condition,
        "batch_mean",
        point_mode="batch",
        force_show_dpis=tuple(virb["DPIS"]),
    )

virb_axes[0].set_ylabel(virb["YLABEL"])
image_panel(virb_scan_ax, virb_scan)

# F. Colony count timecourse
colony["draw_timecourse"](colony_ax, colony_clean)
colony["annotate_timecourse_letters"](
    colony_ax,
    colony_point_contrasts,
    colony_clean,
    force_show_dpis=(0, 1),
)
colony_ax.set_ylabel(colony["YLABEL"])
panel_label(colony_ax, "F.")

for ax in [rna_ax, *gfp_axes, elisa_ax, *virb_axes, colony_ax]:
    apply_axis_style(ax)

for ax in [rna_ax, elisa_ax]:
    make_errorbars_black(ax)

plant_line_legend(fig,
    handles=legend_handles(gfp["LINES"], COLORS, MARKERS),
    frameon=False,
    loc="upper right",
    bbox_to_anchor=(0.985, 0.975),
    labelcolor=TEXT_COLOR,
)

fig.subplots_adjust(left=0.070, right=0.925, top=0.965, bottom=0.055)
top_label_y = rna_image_ax.get_position().y1 + 0.012
fixed_panel_label_at_y(fig, "A.", top_label_y)
fixed_panel_label_at_y(fig, "B.", top_label_y, x=rna_image_ax.get_position().x0 - 0.018)
fixed_panel_label(fig, gfp_axes[0], "C.", dy=0.030)
fixed_panel_label(fig, virb_axes[0], "D.")
fixed_panel_label(fig, elisa_ax, "E.")
save(fig, "combined_photo_rna_gfp_elisa_virb_colony_abcdef")

# %% [markdown]
# # Smaller figure sets

# %%
# Figure 1: plant phenotype, GFP protein expression/western, and ELISA.
fig = plt.figure(figsize=(17.2, 10.2))
outer = fig.add_gridspec(
    2,
    2,
    width_ratios=[1.32, 1.08],
    height_ratios=[1.18, 0.82],
    wspace=0.34,
    hspace=0.34,
)
split_photo_ax = fig.add_subplot(outer[0, 0])
split_elisa_ax = fig.add_subplot(outer[1, 0])
protein_grid = outer[:, 1].subgridspec(
    2,
    2,
    height_ratios=[1.0, 0.62],
    wspace=0.18,
    hspace=0.22,
)
split_gfp_axes = [fig.add_subplot(protein_grid[0, index]) for index in range(2)]
split_western_ax = fig.add_subplot(protein_grid[1, :])

labelled_plant_image_panel(split_photo_ax, plant_photo)
for ax, p19 in zip(split_gfp_axes, gfp["P19_ORDER"]):
    gfp["draw_timecourse"](ax, gfp_rm, p19, "batch_mean", point_mode="batch")
    gfp["annotate_timecourse_letters"](
        ax,
        gfp_point_contrasts,
        gfp_rm,
        p19,
        "batch_mean",
        point_mode="batch",
        force_show_dpis=(0, 1),
    )
split_gfp_axes[0].set_ylabel("GFP signal intensity (x1000)")
for ax in split_gfp_axes:
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _position: f"{value / 1000:g}"))
image_panel_with_opaque_box(split_western_ax, gfp_western)
split_western_ax.set_title("p19", pad=12)

elisa["draw_auc_by_dpi"](split_elisa_ax, elisa_auc)
annotate_auc_letters(
    split_elisa_ax,
    elisa_auc_contrasts,
    elisa_auc,
    force_show_dpis=tuple(elisa["DPIS"]),
)
recolor_grouped_scatter_by_x(split_elisa_ax, COLORS)
split_elisa_ax.set_ylabel(r"AUC of normalised A$_{450}$")
style_bars(split_elisa_ax)

for ax in [*split_gfp_axes, split_elisa_ax]:
    apply_axis_style(ax)
for ax in [*split_gfp_axes, split_western_ax]:
    ax.set_facecolor("none")
make_errorbars_black(split_elisa_ax)
fig.subplots_adjust(left=0.085, right=0.955, top=0.89, bottom=0.09)
# The photo occupies less of its nominal axes height because of its landscape
# aspect ratio. Shorten and centre the GFP axes so the visible top-row panels
# have comparable heights.
for ax in split_gfp_axes:
    position = ax.get_position()
    new_height = position.height * 0.64
    ax.set_position(
        [
            position.x0,
            position.y0 + (position.height - new_height) * 0.42,
            position.width,
            new_height,
        ]
    )
western_position = split_western_ax.get_position()
split_western_ax.set_position(
    [
        western_position.x0 - 0.015,
        western_position.y0,
        western_position.width,
        western_position.height,
    ]
)
rounded_axis_group_box(
    fig,
    [*split_gfp_axes, split_western_ax],
    facecolor="#E8D99B",
    alpha=0.25,
    pad=0.038,
    pad_left=0.075,
    pad_top=0.085,
)
plant_line_legend(fig,
    handles=legend_handles(gfp["LINES"], COLORS, MARKERS),
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.975),
    ncol=len(gfp["LINES"]),
    labelcolor=TEXT_COLOR,
    fontsize=17,
    handlelength=2.8,
    columnspacing=2.2,
    handletextpad=0.8,
)
label_x = 0.032
top_label_y = 0.865
fig.text(label_x, top_label_y, "A", ha="left", va="top", fontsize=20, fontweight="bold", color=TEXT_COLOR)
fig.text(
    split_gfp_axes[0].get_position().x0 - 0.055,
    top_label_y,
    "B",
    ha="left",
    va="top",
    fontsize=20,
    fontweight="bold",
    color=TEXT_COLOR,
)
fig.text(label_x, split_elisa_ax.get_position().y1 + 0.035, "C", ha="left", va="top", fontsize=20, fontweight="bold", color=TEXT_COLOR)
save_split(fig, "plant_gfp_western_elisa")
plt.close(fig)

# %%
# Figure 2: the original RNA panel as a standalone figure.
fig = plt.figure(figsize=(11.2, 4.8))
rna_grid = fig.add_gridspec(1, 2, width_ratios=[1.95, 1.0], wspace=0.16)
split_rna_image_ax = fig.add_subplot(rna_grid[0, 0])
split_rna_ax = fig.add_subplot(rna_grid[0, 1])
image_panel(split_rna_image_ax, rna_formatted)
rna["draw_grouped_bar_by_dpi"](split_rna_ax, rna_clean)
rna["annotate_grouped_bar_letters"](
    split_rna_ax,
    rna_point_contrasts,
    rna_clean,
    force_show_dpis=tuple(rna["DPIS"]),
)
split_rna_ax.set_ylabel(rna["YLABEL"])
style_bars(split_rna_ax)
apply_axis_style(split_rna_ax)
make_errorbars_black(split_rna_ax)
plant_line_legend(fig,
    handles=legend_handles(rna["LINES"], COLORS, MARKERS),
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.76, 0.995),
    ncol=len(rna["LINES"]),
    labelcolor=TEXT_COLOR,
)
fig.subplots_adjust(left=0.035, right=0.965, top=0.80, bottom=0.17)
save_split(fig, "rna")
plt.close(fig)

# %%
# Figure 3: pTac and virB timecourses above the full-width colony-count panel.
fig = plt.figure(figsize=(10.8, 8.4))
virb_colony_grid = fig.add_gridspec(
    2,
    2,
    height_ratios=[1.0, 1.0],
    wspace=0.31,
    hspace=0.48,
)
split_virb_axes = [fig.add_subplot(virb_colony_grid[0, index]) for index in range(2)]
split_colony_ax = fig.add_subplot(virb_colony_grid[1, :])

for ax, condition in zip(split_virb_axes, virb["P19_ORDER"]):
    virb["draw_timecourse"](ax, virb_rm, condition, "batch_mean", point_mode="batch")
    virb["annotate_timecourse_letters"](
        ax,
        virb_point_contrasts,
        virb_rm,
        condition,
        "batch_mean",
        point_mode="batch",
        force_show_dpis=tuple(virb["DPIS"]),
    )
split_virb_axes[0].set_ylabel(virb["YLABEL"])

colony["draw_timecourse"](split_colony_ax, colony_clean)
colony["annotate_timecourse_letters"](
    split_colony_ax,
    colony_point_contrasts,
    colony_clean,
    force_show_dpis=(0, 1),
)
split_colony_ax.set_ylabel(colony["YLABEL"])

for ax in [*split_virb_axes, split_colony_ax]:
    apply_axis_style(ax)
plant_line_legend(fig,
    handles=legend_handles(virb["LINES"], COLORS, MARKERS),
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.995),
    ncol=len(virb["LINES"]),
    labelcolor=TEXT_COLOR,
)
fig.subplots_adjust(left=0.105, right=0.965, top=0.88, bottom=0.09)
fixed_panel_label(fig, split_virb_axes[0], "a", x=0.035, dy=0.035)
fixed_panel_label(
    fig,
    split_virb_axes[1],
    "b",
    x=split_virb_axes[1].get_position().x0 - 0.055,
    dy=0.035,
)
fixed_panel_label(fig, split_colony_ax, "c", x=0.035, dy=0.035)
save_split(fig, "ptac_virb_colony")
plt.close(fig)

print(f"Combined master figure saved to: {SAVE_DIR}")
print(f"Smaller figure sets saved to: {SPLIT_SAVE_DIR}")
