# %% [markdown]
# # VIGS ontogenesis GFP expression analysis
#
# Plant age/leaf-position GFP expression plots for the VIGS ontogenesis dataset.
# The script reads the original Excel workbook, exports raw/cleaned tables and
# statistics logs, and saves PNG/PDF/SVG figures.

# %%
from itertools import combinations
from io import BytesIO
from pathlib import Path
import os
import subprocess
import sys
import warnings
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
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
import matplotlib.patheffects as path_effects
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
from scipy import ndimage
from statsmodels.formula.api import mixedlm, ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning, module=r"scipy\.stats\..*")

if "__file__" in globals():
    HERE = Path(__file__).resolve().parent
elif (Path.cwd() / "VIGS" / "ontogenesis").is_dir():
    HERE = Path.cwd() / "VIGS" / "ontogenesis"
else:
    HERE = Path.cwd()

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

SCRIPT_STEM = Path(__file__).stem if "__file__" in globals() else "202607_ontogenesis"
SAVE_DIR = HERE / f"{SCRIPT_STEM}_figures"
LOG_DIR = HERE / f"{SCRIPT_STEM}_logs"
SPLIT_DIR = HERE / "split_figs"
IMAGE_DIR = HERE / "Images"
SAVE_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)
SPLIT_DIR.mkdir(exist_ok=True)
SAVE_FMTS = ["png", "pdf", "svg"]

SOURCE_FILE = Path(
    "/Users/user/Library/CloudStorage/OneDrive-Nexus365/MolecularFarming/4. Immunity Silencing/PlantAgeandMass/GFPExpressionbyPlantAgeandLeaf/GFPExpbyPlantAgeandLeaf.xlsx"
)
SOURCE_SHEET = "GFPexpandMetricsbyPlantAge"

AGE_ORDER = [4, 5, 6, 7, 8]
AGE_LABELS = {age: f"{age}-week" for age in AGE_ORDER}
AGE_BASE_COLOR = "#3E6E86"
AGE_COLORS = {age: AGE_BASE_COLOR for age in AGE_ORDER}
AGE_MARKERS = {4: "o", 5: "^", 6: "s", 7: "D", 8: "P"}
BLACK_BG = "#000000"
BLACK_TEXT = "#E8E8E8"
BLACK_AXIS = "#D8D8D8"
BLACK_BAR = "#B8C4BB"
BLACK_POINT = "#A8C5D3"
ERRBAR = "SEM"
OUTLIER_ALPHA = 0.05
SIG_ALPHA = 0.05
BAR_ALPHA = 0.72
BAR_LINEWIDTH = 1.0
GFP_Y_SCALE = 1000
TOTAL_GFP_Y_SCALE = 1000
LEAF_GFP_METRIC = "gfp_expressed_x_biomass"
LEAF_GFP_LABEL = "GFP x biomass"


ONTOGENESIS_RCPARAMS = {
    **ANALYSIS_RCPARAMS,
    "font.family": "sans-serif",
    "font.sans-serif": FONT_STACK,
    "font.size": 11.5,
    "axes.titlesize": 12.5,
    "axes.labelsize": 11.5,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "legend.title_fontsize": 11,
}

MASTER_ONTOGENESIS_RCPARAMS = {
    **ONTOGENESIS_RCPARAMS,
    "font.size": ONTOGENESIS_RCPARAMS["font.size"] + 2,
    "axes.titlesize": ONTOGENESIS_RCPARAMS["axes.titlesize"] + 2,
    "axes.labelsize": ONTOGENESIS_RCPARAMS["axes.labelsize"] + 2,
    "xtick.labelsize": ONTOGENESIS_RCPARAMS["xtick.labelsize"] + 2,
    "ytick.labelsize": ONTOGENESIS_RCPARAMS["ytick.labelsize"] + 2,
    "legend.fontsize": ONTOGENESIS_RCPARAMS["legend.fontsize"] + 2,
    "legend.title_fontsize": ONTOGENESIS_RCPARAMS["legend.title_fontsize"] + 2,
}

plt.rcParams.update(ONTOGENESIS_RCPARAMS)


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


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=1.0, length=4.5, pad=3, direction="out")
    apply_axis_style(ax)
    for spine in ax.spines.values():
        spine.set_color(TEXT_COLOR)


def scaled_ylabel(text, scale):
    return f"{text} (x{scale})" if scale != 1 else text


def save(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(SAVE_DIR / f"{output_stem}.{file_format}", dpi=300, bbox_inches="tight")


def save_dark(fig, name):
    output_stem = name if name.startswith(SCRIPT_STEM) else f"{SCRIPT_STEM}_{name}"
    for file_format in SAVE_FMTS:
        fig.savefig(
            SAVE_DIR / f"{output_stem}.{file_format}",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.35,
            facecolor=BLACK_BG,
            edgecolor=BLACK_BG,
        )


def save_split_dark(fig, name):
    """Save split black-background figures as PNG and SVG only."""
    for file_format in ["png", "svg"]:
        fig.savefig(
            SPLIT_DIR / f"{name}.{file_format}",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.35,
            facecolor=BLACK_BG,
            edgecolor=BLACK_BG,
        )


def clear_previous_outputs():
    for folder in [SAVE_DIR, LOG_DIR]:
        for path in folder.glob(f"{SCRIPT_STEM}_*"):
            if path.is_file():
                path.unlink()
    for path in LOG_DIR.glob("ontogenesis_*.xlsx"):
        if path.is_file():
            path.unlink()
    for path in LOG_DIR.glob("ontogenesis_*.txt"):
        if path.is_file():
            path.unlink()




def age_legend(fig, anchor=(0.5, 0.98), ncol=5, title="Plant age"):
    handles = [
        plt.Line2D(
            [0],
            [0],
            color=AGE_COLORS[age],
            marker=AGE_MARKERS[age],
            markerfacecolor=AGE_COLORS[age],
            markeredgecolor="none",
            linewidth=0,
            markersize=5.5,
            label=AGE_LABELS[age],
        )
        for age in AGE_ORDER
    ]
    return fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=anchor,
        ncol=ncol,
        title=title,
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )


def add_panel_label(ax, label, x=-0.08, y=1.04, color=TEXT_COLOR):
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=17,
        fontweight="bold",
        color=color,
        clip_on=False,
    )


def master_panel_label(fig, ax, label, x=0.035, y_pad=0.008):
    bbox = ax.get_position()
    fig.text(
        x,
        bbox.y1 + y_pad,
        label,
        ha="left",
        va="bottom",
        fontsize=17,
        fontweight="bold",
        color=TEXT_COLOR,
    )


def crop_image_content(image, white_threshold=246, pad=18):
    image = ImageOps.exif_transpose(image)
    if image.mode in {"RGBA", "LA"}:
        alpha = image.getchannel("A")
        bbox = alpha.getbbox()
    else:
        rgb = image.convert("RGB")
        array = np.asarray(rgb)
        mask = np.any(array < white_threshold, axis=2)
        if mask.any():
            rows, cols = np.where(mask)
            bbox = (
                max(int(cols.min()) - pad, 0),
                max(int(rows.min()) - pad, 0),
                min(int(cols.max()) + pad + 1, image.width),
                min(int(rows.max()) + pad + 1, image.height),
            )
        else:
            bbox = None
    if bbox:
        image = image.crop(bbox)
    return image.convert("RGB")


def image_array(path, crop=True):
    image = Image.open(path)
    if crop:
        image = crop_image_content(image)
    else:
        image = ImageOps.exif_transpose(image).convert("RGB")
    return np.asarray(image)


def image_region_array(path, bbox):
    image = ImageOps.exif_transpose(Image.open(path)).convert("RGB")
    return np.asarray(image.crop(bbox))


def reordered_leaf_number_image_path():
    source = IMAGE_DIR / "LeafNumber_variance.png"
    output = IMAGE_DIR / "LeafNumber_variance_reordered.png"

    image = ImageOps.exif_transpose(Image.open(source)).convert("RGB")
    canvas = Image.new("RGB", (1850, 856), (0, 0, 0))

    plant = image.crop((95, 180, 690, 770))
    plant.thumbnail((690, 590), Image.Resampling.LANCZOS)
    canvas.paste(plant, (55, 150))

    leaf_crops = {
        1: (815, 560, 1015, 785),
        2: (1010, 565, 1215, 820),
        3: (1198, 555, 1450, 820),
        4: (800, 350, 1015, 590),
        5: (1005, 355, 1220, 595),
        6: (1185, 340, 1450, 585),
        7: (795, 175, 1005, 370),
        8: (995, 190, 1205, 370),
    }
    cell_w, cell_h = 235, 270
    x0, gap = 770, 28
    y_rows = {0: 96, 1: 445}
    for row, leaf_numbers in enumerate([[5, 6, 7, 8], [1, 2, 3, 4]]):
        for col, leaf_number in enumerate(leaf_numbers):
            leaf = image.crop(leaf_crops[leaf_number])
            leaf.thumbnail((cell_w, cell_h), Image.Resampling.LANCZOS)
            x = x0 + col * (cell_w + gap) + (cell_w - leaf.width) // 2
            y = y_rows[row] + (cell_h - leaf.height) // 2
            canvas.paste(leaf, (x, y))

    # The original 1 cm scale bar is replaced by a vector 5 cm bar when the
    # panel is assembled, so it remains crisp and readable at A4 size.
    canvas.save(output)
    return output


def isolate_main_leaf(image):
    """Remove disconnected neighbouring-leaf fragments from a fluorescence crop."""
    array = np.asarray(image.convert("RGB"))
    # The plate background is true black; a low threshold retains the dim blue
    # leaf silhouette as well as the bright fluorescence signal.
    foreground = np.max(array, axis=2) > 4
    labels, count = ndimage.label(foreground)
    if count == 0:
        return image.convert("RGB")
    sizes = ndimage.sum(foreground, labels, range(1, count + 1))
    main_label = int(np.argmax(sizes)) + 1
    keep = ndimage.binary_dilation(labels == main_label, iterations=2)
    # Number glyphs can be disconnected from the leaf silhouette; retain bright
    # neutral pixels only when they sit within the main leaf's bounding box.
    rows, cols = np.where(labels == main_label)
    pad = 16
    within = np.zeros_like(keep)
    within[
        max(rows.min() - pad, 0) : min(rows.max() + pad + 1, array.shape[0]),
        max(cols.min() - pad, 0) : min(cols.max() + pad + 1, array.shape[1]),
    ] = True
    neutral_bright = (array.min(axis=2) > 115) & ((array.max(axis=2) - array.min(axis=2)) < 45)
    keep |= neutral_bright & within
    cleaned = np.zeros_like(array)
    cleaned[keep] = array[keep]
    return Image.fromarray(cleaned)


def standalone_leaf_number_image_path():
    """Create a larger, cleaned A-panel plate for the standalone split figure."""
    source = IMAGE_DIR / "LeafNumber_variance.png"
    output = IMAGE_DIR / "LeafNumber_variance_standalone_cleaned.png"
    image = ImageOps.exif_transpose(Image.open(source)).convert("RGB")
    canvas = Image.new("RGB", (4000, 3400), (0, 0, 0))

    def resize_to_fit(item, max_width, max_height):
        scale = min(max_width / item.width, max_height / item.height)
        size = (max(1, round(item.width * scale)), max(1, round(item.height * scale)))
        return item.resize(size, Image.Resampling.LANCZOS)

    plant = image.crop((95, 180, 690, 770))
    plant = resize_to_fit(plant, 1800, 1800)
    canvas.paste(plant, (50 + (1800 - plant.width) // 2, 760 + (1800 - plant.height) // 2))

    leaf_crops = {
        1: (815, 560, 1015, 785), 2: (1010, 565, 1215, 820),
        3: (1198, 555, 1450, 820), 4: (800, 350, 1015, 590),
        5: (1005, 355, 1220, 595), 6: (1185, 340, 1450, 585),
        7: (795, 175, 1005, 370), 8: (995, 190, 1205, 370),
    }
    cell_w, cell_h = 930, 720
    x0, gap = 2020, 70
    y_rows = {0: 40, 1: 870, 2: 1700, 3: 2530}
    for row, leaf_numbers in enumerate([[5, 6], [7, 8], [1, 2], [3, 4]]):
        for col, leaf_number in enumerate(leaf_numbers):
            leaf = isolate_main_leaf(image.crop(leaf_crops[leaf_number]))
            if leaf_number == 2:
                leaf = leaf.crop((0, 0, leaf.width, round(leaf.height * 0.80)))
            elif leaf_number == 6:
                leaf = leaf.crop((0, 0, leaf.width, round(leaf.height * 0.82)))
            leaf = resize_to_fit(leaf, cell_w, cell_h)
            x = x0 + col * (cell_w + gap) + (cell_w - leaf.width) // 2
            y = y_rows[row] + (cell_h - leaf.height) // 2
            canvas.paste(leaf, (x, y))
    canvas.save(output)
    return output


def vertical_leaf_age_pdf_panel_path():
    """Reassemble the improved lower A panel from page 1 of leafAge.pdf."""
    source = IMAGE_DIR / "leafAge.pdf"
    rendered = IMAGE_DIR / "leafAge_page1_300dpi.png"
    output = IMAGE_DIR / "leafAge_vertical_panel.png"
    if not rendered.exists() and not source.exists():
        if output.exists():
            return output
        raise FileNotFoundError(source)
    if not rendered.exists() or (
        source.exists() and rendered.stat().st_mtime < source.stat().st_mtime
    ):
        prefix = rendered.with_suffix("")
        subprocess.run(
            ["pdftoppm", "-f", "1", "-singlefile", "-png", "-r", "300", str(source), str(prefix)],
            check=True,
        )
    page = Image.open(rendered).convert("RGB")
    # Coordinates were measured on the 150 dpi preview and doubled here.
    crop_150 = {
        "plant": (45, 950, 380, 1300),
        1: (425, 1070, 585, 1285), 2: (565, 1060, 725, 1290),
        3: (705, 1055, 880, 1295), 4: (855, 1055, 1035, 1295),
        5: (435, 895, 595, 1065), 6: (575, 890, 740, 1065),
        7: (720, 890, 885, 1065), 8: (865, 890, 1035, 1065),
    }
    crops = {
        key: page.crop(tuple(2 * value for value in bbox))
        for key, bbox in crop_150.items()
    }
    canvas = Image.new("RGB", (3600, 6100), (0, 0, 0))

    def fit(item, width, height):
        scale = min(width / item.width, height / item.height)
        return item.resize((round(item.width * scale), round(item.height * scale)), Image.Resampling.LANCZOS)

    plant = fit(crops["plant"], 3200, 2450)
    canvas.paste(plant, ((3600 - plant.width) // 2, 40))
    cell_w, cell_h = 1660, 880
    for column, numbers in enumerate([[1, 2, 3, 4], [5, 6, 7, 8]]):
        for row, number in enumerate(numbers):
            leaf = fit(crops[number], cell_w, cell_h)
            x = 80 + column * 1760 + (cell_w - leaf.width) // 2
            y = 2470 + row * 900 + (cell_h - leaf.height) // 2
            canvas.paste(leaf, (x, y))
    canvas.save(output)
    return output


def side_plant_image_path():
    source = IMAGE_DIR / "PlantImagesbyAge_Compiled.png"
    output = IMAGE_DIR / "PlantImagesbyAge_SideOnly.png"
    if output.exists() and output.stat().st_mtime >= source.stat().st_mtime:
        return output
    image = ImageOps.exif_transpose(Image.open(source)).convert("RGB")
    width, height = image.size
    side = image.crop((0, 0, width, int(height * 0.56)))
    side.save(output)
    return output


def side_plant_crop(age):
    image = ImageOps.exif_transpose(Image.open(side_plant_image_path())).convert("RGB")
    plant_regions = {
        4: (120, 575, 625, 952),
        5: (665, 520, 1085, 952),
        6: (1088, 500, 1690, 952),
        7: (1690, 130, 2238, 952),
        8: (2238, 52, 2788, 952),
    }
    return np.asarray(image.crop(plant_regions[age]))


def plot_aligned_side_plants_dark(
    ax, age_positions, xlim, image_zoom=0.25, zoom_by_age=None, y_by_age=None
):
    ax.set_facecolor(BLACK_BG)
    ax.set_xlim(xlim)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    for age in AGE_ORDER:
        image = side_plant_crop(age)
        zoom = image_zoom * ((zoom_by_age or {}).get(age, 1.0))
        imagebox = OffsetImage(image, zoom=zoom, interpolation="lanczos")
        imagebox.set_clip_on(False)
        annotation = AnnotationBbox(
            imagebox,
            (age_positions[age], (y_by_age or {}).get(age, 0.02)),
            xycoords="data",
            box_alignment=(0.5, 0.0),
            frameon=False,
            pad=0,
            annotation_clip=False,
        )
        annotation.set_clip_on(False)
        ax.add_artist(annotation)
    return ax


def add_dark_scale_bar(
    ax, length_axes, label="5 cm", x_right=0.975, y=0.055,
    fontsize=32, linewidth=5.0,
):
    """Add a high-contrast horizontal scale bar in axes coordinates."""
    x_left = x_right - length_axes
    ax.plot(
        [x_left, x_right],
        [y, y],
        transform=ax.transAxes,
        color=BLACK_TEXT,
        linewidth=linewidth,
        solid_capstyle="butt",
        clip_on=False,
        zorder=20,
    )
    ax.text(
        (x_left + x_right) / 2,
        y - 0.035,
        label,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=fontsize,
        color=BLACK_TEXT,
        clip_on=False,
        zorder=20,
    )


def plot_black_age_labels(ax, age_positions, xlim, fontsize=36, suffix_by_age=None):
    ax.set_axis_off()
    ax.set_facecolor(BLACK_BG)
    ax.set_xlim(xlim)
    ax.set_ylim(0, 1)
    for age in AGE_ORDER:
        ax.text(
            age_positions[age],
            0.45,
            f"{age} weeks{(suffix_by_age or {}).get(age, '')}",
            transform=ax.transData,
            ha="center",
            va="center",
            fontsize=fontsize,
            color=BLACK_TEXT,
        )


def age_image_path(age, view):
    matches = sorted(IMAGE_DIR.glob(f"*{age}-week*_{view}.tif*"))
    if not matches:
        raise FileNotFoundError(f"No {view} image found for {age}-week plants in {IMAGE_DIR}")
    return matches[0]


def label_age_image(ax, age):
    text = ax.text(
        0.05,
        0.92,
        f"{age}-week",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        fontfamily="Helvetica Neue",
        fontweight="light",
        color="white",
    )
    text.set_path_effects([path_effects.withStroke(linewidth=1.0, foreground="black", alpha=0.45)])


def plot_age_image_grid(fig, spec, add_row_labels=True):
    grid = GridSpecFromSubplotSpec(
        2,
        len(AGE_ORDER),
        subplot_spec=spec,
        wspace=0.03,
        hspace=0.03,
    )
    first_ax = None
    for row, view in enumerate(["top", "side"]):
        for col, age in enumerate(AGE_ORDER):
            ax = fig.add_subplot(grid[row, col])
            if first_ax is None:
                first_ax = ax
            ax.imshow(image_array(age_image_path(age, view), crop=True))
            label_age_image(ax, age)
            ax.set_axis_off()
            if add_row_labels and col == 0:
                ax.text(
                    -0.05,
                    0.5,
                    view.capitalize(),
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    rotation=90,
                    fontsize=8.5,
                    color=TEXT_COLOR,
                    clip_on=False,
                )
    return first_ax


def plot_compiled_plant_images(ax):
    compiled_path = IMAGE_DIR / "PlantImagesbyAge_Compiled.png"
    if not compiled_path.exists():
        raise FileNotFoundError(
            "Expected VIGS/Ontogenesis/Images/PlantImagesbyAge_Compiled.png. "
            "A PDF source is present, but this script embeds the rasterized PNG copy."
        )
    ax.imshow(image_array(compiled_path, crop=False))
    ax.set_axis_off()
    return ax


# %%
def load_raw_data():
    df = pd.read_excel(SOURCE_FILE, sheet_name=SOURCE_SHEET, header=3)
    if "Sample ID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Sample ID"})
    df = df.dropna(subset=["Sample ID", "Plant Age", "Replicate", "Leaf Number"]).copy()
    for column in [
        "Replicate",
        "Plant Age",
        "Leaf Number",
        "Stem Height",
        "Number of Leaves",
        "Total Plant Weight",
        "Leaf weight",
        "Total weight of true leaves",
        "GFP expression",
        "GFP expressed x biomass",
        "Total plant GFP expressed",
    ]:
        df[column] = pd.to_numeric(df[column].replace("-", np.nan), errors="coerce")

    df["plant_age"] = df["Plant Age"].astype("Int64")
    df["replicate"] = df["Replicate"].astype("Int64")
    df["leaf_number"] = df["Leaf Number"].astype("Int64")
    df["plant_id"] = (
        df["plant_age"].astype(str) + "wk_rep" + df["replicate"].astype(str)
    )
    df["sample_id"] = df["Sample ID"].astype(str)
    return df


def make_plant_metrics(raw):
    metric_columns = [
        "Stem Height",
        "Number of Leaves",
        "Total Plant Weight",
        "Total weight of true leaves",
        "Total plant GFP expressed",
    ]
    plant = (
        raw.sort_values(["plant_age", "replicate", "leaf_number"])
        .groupby(["plant_age", "replicate", "plant_id"], as_index=False, observed=True)[
            metric_columns
        ]
        .first()
    )
    return plant.rename(
        columns={
            "Stem Height": "stem_height",
            "Number of Leaves": "number_of_leaves",
            "Total Plant Weight": "total_plant_weight",
            "Total weight of true leaves": "total_true_leaf_weight",
            "Total plant GFP expressed": "total_plant_gfp_expressed",
        }
    )


def make_leaf_data(raw):
    leaf = raw.rename(
        columns={
            "Leaf weight": "leaf_weight",
            "GFP expression": "gfp_expression",
            "GFP expressed x biomass": "gfp_expressed_x_biomass",
        }
    )[
        [
            "sample_id",
            "plant_id",
            "plant_age",
            "replicate",
            "leaf_number",
            "leaf_weight",
            "gfp_expression",
            "gfp_expressed_x_biomass",
        ]
    ].copy()
    return leaf


def grubbs_two_sided_flag(values):
    values = values.dropna()
    n = len(values)
    if n < 4:
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


def flag_outliers(data, metric_columns, id_columns):
    rows = []
    screened = data.copy()
    for metric in metric_columns:
        screened[f"{metric}_flag_outlier"] = False

    for metric in metric_columns:
        for age, group in screened.groupby("plant_age", observed=True):
            result = grubbs_two_sided_flag(group[metric])
            if result is None or result["p_value"] >= OUTLIER_ALPHA:
                continue
            idx = result["index"]
            screened.loc[idx, f"{metric}_flag_outlier"] = True
            rows.append(
                {
                    **{column: screened.loc[idx, column] for column in id_columns},
                    "metric": metric,
                    "plant_age": age,
                    "value": screened.loc[idx, metric],
                    "method": "two-sided Grubbs test",
                    "group_n": result["n"],
                    "group_mean": result["group_mean"],
                    "group_sd": result["group_sd"],
                    "g_statistic": result["g_statistic"],
                    "p_value": result["p_value"],
                    "alpha": OUTLIER_ALPHA,
                    "manual_exclude": False,
                }
            )
    return screened, pd.DataFrame(rows)


# %%
def welch_anova_from_groups(groups):
    cleaned = [pd.Series(group).dropna().astype(float) for group in groups]
    cleaned = [group for group in cleaned if len(group) >= 2]
    k = len(cleaned)
    if k < 2:
        return {"F": np.nan, "df1": np.nan, "df2": np.nan, "p_value": np.nan}

    n = np.array([len(group) for group in cleaned], dtype=float)
    means = np.array([group.mean() for group in cleaned], dtype=float)
    variances = np.array([group.var(ddof=1) for group in cleaned], dtype=float)
    if np.any(variances <= 0) or np.any(pd.isna(variances)):
        result = stats.f_oneway(*cleaned)
        return {
            "F": float(result.statistic),
            "df1": k - 1,
            "df2": sum(n) - k,
            "p_value": float(result.pvalue),
            "fallback": "classic_one_way_anova_zero_variance",
        }

    weights = n / variances
    weight_sum = weights.sum()
    weighted_mean = np.sum(weights * means) / weight_sum
    df1 = k - 1
    numerator = np.sum(weights * (means - weighted_mean) ** 2) / df1
    lambda_term = np.sum((1 / (n - 1)) * (1 - weights / weight_sum) ** 2)
    denominator = 1 + (2 * (k - 2) / (k**2 - 1)) * lambda_term
    f_stat = numerator / denominator
    df2 = (k**2 - 1) / (3 * lambda_term)
    p_value = stats.f.sf(f_stat, df1, df2)
    return {"F": float(f_stat), "df1": float(df1), "df2": float(df2), "p_value": float(p_value)}


def run_age_stats(data, metric, label):
    ready = data.dropna(subset=["plant_age", metric]).copy()
    groups = [ready.loc[ready["plant_age"] == age, metric] for age in AGE_ORDER]
    anova = welch_anova_from_groups(groups)
    anova_row = {
        "analysis": label,
        "metric": metric,
        "test": "Welch one-way ANOVA",
        **anova,
    }

    assumption_rows = []
    for age in AGE_ORDER:
        values = ready.loc[ready["plant_age"] == age, metric].dropna()
        shapiro_p = np.nan
        if 3 <= len(values) <= 5000:
            shapiro_p = stats.shapiro(values).pvalue
        assumption_rows.append(
            {
                "analysis": label,
                "metric": metric,
                "plant_age": age,
                "n": len(values),
                "mean": values.mean(),
                "sd": values.std(ddof=1) if len(values) > 1 else np.nan,
                "shapiro_p": shapiro_p,
            }
        )
    valid_groups = [group.dropna() for group in groups if len(group.dropna()) >= 2]
    levene_p = stats.levene(*valid_groups, center="median").pvalue if len(valid_groups) >= 2 else np.nan
    for row in assumption_rows:
        row["levene_p_across_age"] = levene_p

    pair_rows = []
    p_values = []
    for left, right in combinations(AGE_ORDER, 2):
        left_values = ready.loc[ready["plant_age"] == left, metric].dropna()
        right_values = ready.loc[ready["plant_age"] == right, metric].dropna()
        if len(left_values) < 2 or len(right_values) < 2:
            continue
        test = stats.ttest_ind(left_values, right_values, equal_var=False, nan_policy="omit")
        pair_rows.append(
            {
                "analysis": label,
                "metric": metric,
                "contrast": f"{left}-week - {right}-week",
                "left_age": left,
                "right_age": right,
                "left_n": len(left_values),
                "right_n": len(right_values),
                "left_mean": left_values.mean(),
                "right_mean": right_values.mean(),
                "estimate_left_minus_right": left_values.mean() - right_values.mean(),
                "t": float(test.statistic),
                "p_raw": float(test.pvalue),
            }
        )
        p_values.append(float(test.pvalue))
    pairwise = pd.DataFrame(pair_rows)
    if not pairwise.empty:
        _, adjusted, _, _ = multipletests(p_values, method="holm")
        pairwise["p_holm_within_metric"] = adjusted
        pairwise["stars"] = [p_to_stars(p) for p in adjusted]

    return pd.DataFrame([anova_row]), pd.DataFrame(assumption_rows), pairwise


def run_leaf_gfp_stats(leaf_data):
    ready = leaf_data.dropna(subset=["plant_age", "leaf_number", LEAF_GFP_METRIC, "plant_id"]).copy()
    ready["plant_age"] = pd.Categorical(ready["plant_age"], categories=AGE_ORDER)
    ready["leaf_number"] = pd.Categorical(ready["leaf_number"], categories=sorted(ready["leaf_number"].dropna().unique()))
    ready["leaf_gfp_value"] = ready[LEAF_GFP_METRIC]

    model_rows = []
    model_note = ""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = mixedlm(
                "leaf_gfp_value ~ C(plant_age) + C(leaf_number)",
                ready,
                groups=ready["plant_id"],
            ).fit(reml=False, method="lbfgs")
        model_rows = pd.DataFrame(
            {
                "term": model.params.index,
                "estimate": model.params.values,
                "p_value": model.pvalues.reindex(model.params.index).values,
            }
        )
        model_note = str(model.summary())
    except Exception as exc:
        fallback = ols("leaf_gfp_value ~ C(plant_age) + C(leaf_number)", data=ready).fit()
        table = anova_lm(fallback, typ=2).reset_index(names="term")
        model_rows = table
        model_note = f"MixedLM failed ({exc}); exported type-II OLS ANOVA fallback.\n\n{fallback.summary()}"

    pair_rows = []
    for leaf_number, subset in leaf_data.dropna(subset=[LEAF_GFP_METRIC]).groupby("leaf_number", observed=True):
        p_values = []
        start = len(pair_rows)
        for left, right in combinations(AGE_ORDER, 2):
            left_values = subset.loc[subset["plant_age"] == left, LEAF_GFP_METRIC].dropna()
            right_values = subset.loc[subset["plant_age"] == right, LEAF_GFP_METRIC].dropna()
            if len(left_values) < 2 or len(right_values) < 2:
                continue
            test = stats.ttest_ind(left_values, right_values, equal_var=False, nan_policy="omit")
            pair_rows.append(
                {
                    "leaf_number": leaf_number,
                    "contrast": f"{left}-week - {right}-week",
                    "left_age": left,
                    "right_age": right,
                    "left_n": len(left_values),
                    "right_n": len(right_values),
                    "left_mean": left_values.mean(),
                    "right_mean": right_values.mean(),
                    "estimate_left_minus_right": left_values.mean() - right_values.mean(),
                    "t": float(test.statistic),
                    "p_raw": float(test.pvalue),
                }
            )
            p_values.append(float(test.pvalue))
        if p_values:
            _, adjusted, _, _ = multipletests(p_values, method="holm")
            for offset, adjusted_p in enumerate(adjusted):
                pair_rows[start + offset]["p_holm_within_leaf"] = adjusted_p
                pair_rows[start + offset]["stars"] = p_to_stars(adjusted_p)

    return model_rows, model_note, pd.DataFrame(pair_rows)


def run_all_stats(plant_data, leaf_data):
    anova_tables = []
    assumption_tables = []
    pairwise_tables = []
    for metric, label in [
        ("stem_height", "height_per_age"),
        ("number_of_leaves", "number_of_leaves_per_age"),
        ("total_plant_gfp_expressed", "total_plant_gfp_expressed_per_age"),
    ]:
        anova, assumptions, pairwise = run_age_stats(plant_data, metric, label)
        anova_tables.append(anova)
        assumption_tables.append(assumptions)
        pairwise_tables.append(pairwise)

    anova = pd.concat(anova_tables, ignore_index=True)
    assumptions = pd.concat(assumption_tables, ignore_index=True)
    pairwise = pd.concat([x for x in pairwise_tables if not x.empty], ignore_index=True)
    leaf_model, leaf_model_note, leaf_pairwise = run_leaf_gfp_stats(leaf_data)

    anova.to_excel(LOG_DIR / "ontogenesis_age_welch_anova.xlsx", index=False)
    assumptions.to_excel(LOG_DIR / "ontogenesis_age_assumption_checks.xlsx", index=False)
    pairwise.to_excel(LOG_DIR / "ontogenesis_age_pairwise_contrasts.xlsx", index=False)
    leaf_model.to_excel(LOG_DIR / "ontogenesis_leaf_gfp_model_terms.xlsx", index=False)
    leaf_pairwise.to_excel(LOG_DIR / "ontogenesis_leaf_gfp_pairwise_age_by_leaf.xlsx", index=False)
    (LOG_DIR / "ontogenesis_leaf_gfp_model_summary.txt").write_text(leaf_model_note)
    return anova, assumptions, pairwise, leaf_model, leaf_pairwise


def n_label_from_counts(counts):
    values = [count for count in counts.values() if count > 0]
    if not values:
        return ""
    low = min(values)
    high = max(values)
    return f"n={low}" if low == high else f"n={low}-{high}"


def annotation_fontsize():
    return max(9.5, plt.rcParams["xtick.labelsize"])


def plot_age_bar(
    ax,
    data,
    metric,
    ylabel,
    title,
    pairwise,
    y_scale=1,
    age_positions=None,
    xlim=None,
    bar_width=0.58,
):
    ymax = 0
    counts = {}
    positions = (
        np.array([age_positions[age] for age in AGE_ORDER], dtype=float)
        if age_positions is not None
        else np.arange(len(AGE_ORDER), dtype=float)
    )
    for i, age in enumerate(AGE_ORDER):
        values = data.loc[data["plant_age"] == age, metric].dropna() / y_scale
        if values.empty:
            continue
        mean = values.mean()
        error = sem(values)
        ymax = max(ymax, float(values.max()), float(mean + error))
        counts[age] = len(values)
        ax.bar(
            positions[i],
            mean,
            width=bar_width,
            color=AGE_COLORS[age],
            alpha=BAR_ALPHA,
            edgecolor="black",
            linewidth=BAR_LINEWIDTH,
            zorder=2,
        )
        ax.errorbar(
            positions[i],
            mean,
            yerr=error,
            fmt="none",
            ecolor=ERRORBAR_COLOR,
            elinewidth=1.0,
            capsize=3.2,
            capthick=1.0,
            zorder=5,
        )
        jitter = np.linspace(-bar_width * 0.17, bar_width * 0.17, len(values))
        ax.scatter(
            positions[i] + jitter,
            values,
            s=22,
            marker=AGE_MARKERS[age],
            facecolor=AGE_COLORS[age],
            edgecolor=BLACK_TEXT,
            linewidth=0.55,
            alpha=INDIVIDUAL_POINT_ALPHA,
            zorder=4,
        )

    n_text = n_label_from_counts(counts)
    if n_text:
        ax.text(
            0.98,
            0.97,
            n_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=annotation_fontsize(),
            color=TEXT_COLOR,
        )
    ax.set_xticks(positions)
    ax.set_xticklabels([str(age) for age in AGE_ORDER])
    ax.set_xlabel("Age (weeks)")
    ax.set_ylabel(scaled_ylabel(ylabel, y_scale))
    ax.set_title("")
    ax.set_ylim(0, ymax * 1.24 if ymax else 1)
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax)


def leaf_group_geometry(leaf_data):
    ready = leaf_data.dropna(subset=[LEAF_GFP_METRIC]).copy()
    width = 0.11
    bar_step = width * 1.55
    group_gap = 0.56
    leaf_numbers_by_age = {
        age: sorted(
            ready.loc[ready["plant_age"] == age, "leaf_number"]
            .dropna()
            .astype(int)
            .unique()
        )
        for age in AGE_ORDER
    }
    group_centers = {}
    cursor = 0.0
    for age in AGE_ORDER:
        leaf_count = len(leaf_numbers_by_age[age])
        if leaf_count == 0:
            continue
        half_span = ((leaf_count - 1) / 2) * bar_step + (width / 2)
        group_centers[age] = cursor + half_span
        cursor += (half_span * 2) + group_gap

    group_bounds = {}
    for age in AGE_ORDER:
        leaf_numbers = leaf_numbers_by_age[age]
        if not leaf_numbers:
            continue
        offsets = (np.arange(len(leaf_numbers)) - (len(leaf_numbers) - 1) / 2) * bar_step
        group_bounds[age] = (
            group_centers[age] + offsets[0],
            group_centers[age] + offsets[-1],
        )
    return ready, width, bar_step, leaf_numbers_by_age, group_centers, group_bounds


def leaf_alignment_limits(group_bounds):
    left = min(left for left, _ in group_bounds.values()) - 1.25
    right = max(right for _, right in group_bounds.values()) + 0.46
    return left, right


def black_leaf_group_geometry(leaf_data):
    ready = leaf_data.dropna(subset=[LEAF_GFP_METRIC]).copy()
    width = 0.085
    bar_step = 0.150
    group_spacing = 1.95
    leaf_numbers_by_age = {
        age: sorted(
            ready.loc[ready["plant_age"] == age, "leaf_number"]
            .dropna()
            .astype(int)
            .unique()
        )
        for age in AGE_ORDER
    }
    group_centers = {age: i * group_spacing for i, age in enumerate(AGE_ORDER)}
    # Give the two widest late-stage groups a little extra visual separation.
    group_centers[8] += 0.38
    group_bounds = {}
    for age in AGE_ORDER:
        leaf_numbers = leaf_numbers_by_age[age]
        if not leaf_numbers:
            continue
        offsets = (np.arange(len(leaf_numbers)) - (len(leaf_numbers) - 1) / 2) * bar_step
        group_bounds[age] = (
            group_centers[age] + offsets[0],
            group_centers[age] + offsets[-1],
        )
    return ready, width, bar_step, leaf_numbers_by_age, group_centers, group_bounds


def black_alignment_limits(group_bounds):
    left = min(left for left, _ in group_bounds.values()) - 0.78
    right = max(right for _, right in group_bounds.values()) + 0.52
    return left, right


def plot_leaf_gfp_grouped(ax, leaf_data):
    ready, width, bar_step, leaf_numbers_by_age, group_centers, group_bounds = leaf_group_geometry(leaf_data)

    ymax = 0
    lower = min(-78.0, ready[LEAF_GFP_METRIC].min() / GFP_Y_SCALE * 1.35)
    leaf_label_y = lower * 0.14
    group_line_y = lower * 0.34
    age_label_y = lower * 0.56
    age_axis_label_y = lower * 0.84
    first_leaf_x = None
    counts = {}

    for age in AGE_ORDER:
        age_subset = ready[ready["plant_age"] == age]
        leaf_numbers = leaf_numbers_by_age[age]
        if not leaf_numbers:
            continue
        offsets = (np.arange(len(leaf_numbers)) - (len(leaf_numbers) - 1) / 2) * bar_step
        for leaf_i, (leaf_number, offset) in enumerate(zip(leaf_numbers, offsets)):
            values = age_subset.loc[
                age_subset["leaf_number"] == leaf_number,
                LEAF_GFP_METRIC,
            ].dropna() / GFP_Y_SCALE
            if values.empty:
                continue
            x = group_centers[age] + offset
            first_leaf_x = x if first_leaf_x is None else min(first_leaf_x, x)
            mean = values.mean()
            error = sem(values)
            ymax = max(ymax, float(values.max()), float(mean + error))
            counts[(age, leaf_number)] = len(values)
            if len(leaf_numbers) == 1:
                leaf_alpha = 0.72
            else:
                leaf_alpha = 0.25 + 0.70 * (leaf_i / (len(leaf_numbers) - 1))
            ax.bar(
                x,
                mean,
                width=width,
                color=AGE_COLORS[age],
                alpha=leaf_alpha,
                edgecolor="black",
                linewidth=0.9,
                zorder=2,
            )
            ax.errorbar(
                x,
                mean,
                yerr=error,
                fmt="none",
                ecolor=ERRORBAR_COLOR,
                elinewidth=0.85,
                capsize=2.2,
                capthick=0.85,
                zorder=5,
            )
            jitter = np.linspace(-width * 0.16, width * 0.16, len(values))
            ax.scatter(
                x + jitter,
                values,
                s=14,
                marker=AGE_MARKERS[age],
                facecolor=AGE_COLORS[age],
                edgecolor=BLACK_TEXT,
                linewidth=0.45,
                alpha=max(0.45, min(0.90, leaf_alpha + 0.12)),
                zorder=4,
            )
            ax.text(
                x,
                leaf_label_y,
                str(leaf_number),
                ha="center",
                va="center",
                fontsize=plt.rcParams["xtick.labelsize"],
                rotation=0,
                color=TEXT_COLOR,
                clip_on=False,
            )

    ax.axhline(0, color=TEXT_COLOR, linewidth=0.8, zorder=1)
    tick_length = abs(lower) * 0.055
    for age in AGE_ORDER:
        if age not in group_bounds:
            continue
        left, right = group_bounds[age]
        ax.plot([left, right], [group_line_y, group_line_y], color=TEXT_COLOR, linewidth=1.0, clip_on=False)
        ax.plot(
            [group_centers[age], group_centers[age]],
            [group_line_y, group_line_y - tick_length],
            color=TEXT_COLOR,
            linewidth=1.0,
            clip_on=False,
        )
        ax.text(
            group_centers[age],
            age_label_y,
            str(age),
            ha="center",
            va="center",
            fontsize=plt.rcParams["xtick.labelsize"],
            color=TEXT_COLOR,
            clip_on=False,
        )
    plot_left = min(left for left, _ in group_bounds.values())
    plot_right = max(right for _, right in group_bounds.values())
    ax.text(
        (plot_left + plot_right) / 2,
        age_axis_label_y,
        "Age (weeks)",
        ha="center",
        va="center",
        fontsize=plt.rcParams["axes.labelsize"],
        color=TEXT_COLOR,
        clip_on=False,
    )
    ax.set_xticks([])
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.set_xlabel("")
    ax.set_ylabel(scaled_ylabel(LEAF_GFP_LABEL, GFP_Y_SCALE))
    ax.set_title("")
    upper = ymax * 1.18 if ymax else 1
    ax.set_ylim(lower, upper)
    ax.set_xlim(leaf_alignment_limits(group_bounds))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    positive_ticks = [tick for tick in ax.get_yticks() if tick >= 0]
    ax.set_yticks(positive_ticks)
    n_text = n_label_from_counts(counts)
    if n_text:
        ax.text(
            1.01,
            0.97,
            n_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=annotation_fontsize(),
            color=TEXT_COLOR,
            clip_on=False,
        )
    ax.text(
        (first_leaf_x if first_leaf_x is not None else min(group_centers.values())) - 0.20,
        leaf_label_y,
        "Leaf number:",
        ha="right",
        va="center",
        fontsize=plt.rcParams["xtick.labelsize"],
        color=TEXT_COLOR,
        clip_on=False,
    )
    style_axis(ax)
    ax.spines["bottom"].set_visible(False)
    if positive_ticks:
        ax.spines["left"].set_bounds(0, max(positive_ticks))


def plot_leaf_number_variance_panel(fig, spec):
    image_ax = fig.add_subplot(spec)
    image_ax.imshow(image_array(reordered_leaf_number_image_path(), crop=False))
    image_ax.set_axis_off()
    return image_ax


def style_axis_dark(ax):
    ax.set_facecolor(BLACK_BG)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_color(BLACK_AXIS)
    ax.spines["bottom"].set_color(BLACK_AXIS)
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(colors=BLACK_TEXT, width=1.1, length=5.0, pad=4, direction="out")
    ax.xaxis.label.set_color(BLACK_TEXT)
    ax.yaxis.label.set_color(BLACK_TEXT)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_color(BLACK_TEXT)
        label.set_fontfamily(FONT_STACK)
    for text in ax.texts:
        text.set_color(BLACK_TEXT)
        text.set_fontfamily(FONT_STACK)


def plot_age_bar_dark(
    ax,
    data,
    metric,
    ylabel,
    y_scale=1,
    age_positions=None,
    xlim=None,
    bar_width=0.42,
    bar_color=BLACK_BAR,
    y_nbins=5,
    exact_y_tick_count=None,
):
    ymax = 0
    counts = {}
    positions = (
        np.array([age_positions[age] for age in AGE_ORDER], dtype=float)
        if age_positions is not None
        else np.arange(len(AGE_ORDER), dtype=float)
    )
    for i, age in enumerate(AGE_ORDER):
        values = data.loc[data["plant_age"] == age, metric].dropna() / y_scale
        if values.empty:
            continue
        mean = values.mean()
        error = sem(values)
        ymax = max(ymax, float(values.max()), float(mean + error))
        counts[age] = len(values)
        ax.bar(
            positions[i],
            mean,
            width=bar_width,
            color=bar_color,
            alpha=0.92,
            edgecolor=bar_color,
            linewidth=0.9,
            zorder=2,
        )
        ax.errorbar(
            positions[i],
            mean,
            yerr=error,
            fmt="none",
            ecolor=BLACK_TEXT,
            elinewidth=1.0,
            capsize=3.0,
            capthick=1.0,
            zorder=5,
        )
        jitter = np.linspace(-bar_width * 0.16, bar_width * 0.16, len(values))
        ax.scatter(
            positions[i] + jitter,
            values,
            s=24,
            marker=AGE_MARKERS[age],
            facecolor=BLACK_POINT,
            edgecolor="none",
            alpha=0.55,
            zorder=4,
        )

    n_text = n_label_from_counts(counts)
    if n_text:
        ax.text(
            1.01,
            0.96,
            n_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=annotation_fontsize(),
            color=BLACK_TEXT,
            clip_on=False,
        )
    ax.set_xticks(positions)
    ax.set_xticklabels([str(age) for age in AGE_ORDER])
    ax.set_xlabel("Age (weeks)", labelpad=10)
    ax.set_ylabel(scaled_ylabel(ylabel, y_scale), labelpad=18)
    ax.set_ylim(0, ymax * 1.24 if ymax else 1)
    if xlim is not None:
        ax.set_xlim(xlim)
    if exact_y_tick_count:
        tick_candidates = MaxNLocator(nbins=exact_y_tick_count - 1).tick_values(0, ymax * 1.18)
        tick_candidates = [tick for tick in tick_candidates if tick >= 0]
        top_tick = tick_candidates[-1]
        ax.set_ylim(0, top_tick)
        ax.set_yticks(np.linspace(0, top_tick, exact_y_tick_count))
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins))
    style_axis_dark(ax)


def plot_leaf_gfp_grouped_dark(
    ax, leaf_data, leaf_label_rotation=0, leaf_label_fontsize=None,
    leaf_labels_odd_only=False, leaf_label_y_override=None,
    leaf_number_heading="Leaf number:", leaf_label_every=None,
):
    ready, width, bar_step, leaf_numbers_by_age, group_centers, group_bounds = black_leaf_group_geometry(leaf_data)

    ymax = 0
    leaf_label_y = -0.080 if leaf_label_y_override is None else leaf_label_y_override
    group_line_y = -0.270
    age_label_y = -0.350
    age_axis_label_y = -0.500
    first_leaf_x = None
    counts = {}
    label_transform = ax.get_xaxis_transform()

    for age in AGE_ORDER:
        age_subset = ready[ready["plant_age"] == age]
        leaf_numbers = leaf_numbers_by_age[age]
        if not leaf_numbers:
            continue
        offsets = (np.arange(len(leaf_numbers)) - (len(leaf_numbers) - 1) / 2) * bar_step
        for leaf_i, (leaf_number, offset) in enumerate(zip(leaf_numbers, offsets)):
            values = age_subset.loc[
                age_subset["leaf_number"] == leaf_number,
                LEAF_GFP_METRIC,
            ].dropna() / GFP_Y_SCALE
            if values.empty:
                continue
            x = group_centers[age] + offset
            first_leaf_x = x if first_leaf_x is None else min(first_leaf_x, x)
            mean = values.mean()
            error = sem(values)
            ymax = max(ymax, float(values.max()), float(mean + error))
            counts[(age, leaf_number)] = len(values)
            leaf_alpha = 0.30 + 0.58 * (leaf_i / max(len(leaf_numbers) - 1, 1))
            ax.bar(
                x,
                mean,
                width=width,
                color=BLACK_BAR,
                alpha=leaf_alpha,
                edgecolor=BLACK_AXIS,
                linewidth=0.75,
                zorder=2,
            )
            ax.errorbar(
                x,
                mean,
                yerr=error,
                fmt="none",
                ecolor=BLACK_TEXT,
                elinewidth=0.8,
                capsize=2.0,
                capthick=0.8,
                zorder=5,
            )
            jitter = np.linspace(-width * 0.16, width * 0.16, len(values))
            ax.scatter(
                x + jitter,
                values,
                s=13,
                marker=AGE_MARKERS[age],
                facecolor=BLACK_POINT,
                edgecolor=BLACK_TEXT,
                linewidth=0.45,
                alpha=0.62,
                zorder=4,
            )
            show_leaf_label = (
                (leaf_number - 1) % leaf_label_every == 0
                if leaf_label_every
                else (not leaf_labels_odd_only or leaf_number % 2 == 1)
            )
            if show_leaf_label:
                label_row = 0 if leaf_label_every else (
                    (leaf_i // 2) % 2 if leaf_labels_odd_only else leaf_i % 2
                )
                ax.text(
                    x,
                    leaf_label_y - 0.095 * label_row,
                    str(leaf_number),
                    transform=label_transform,
                    ha="center",
                    va="center",
                    fontsize=leaf_label_fontsize or plt.rcParams["xtick.labelsize"],
                    rotation=leaf_label_rotation,
                    color=BLACK_TEXT,
                    clip_on=False,
                )

    ax.axhline(0, color=BLACK_AXIS, linewidth=1.0, zorder=1)
    tick_length = 0.035
    for age in AGE_ORDER:
        if age not in group_bounds:
            continue
        left, right = group_bounds[age]
        ax.plot(
            [left, right],
            [group_line_y, group_line_y],
            color=BLACK_AXIS,
            linewidth=1.0,
            transform=label_transform,
            clip_on=False,
        )
        ax.plot(
            [group_centers[age], group_centers[age]],
            [group_line_y, group_line_y - tick_length],
            color=BLACK_AXIS,
            linewidth=1.0,
            transform=label_transform,
            clip_on=False,
        )
        ax.text(
            group_centers[age],
            age_label_y,
            str(age),
            transform=label_transform,
            ha="center",
            va="center",
            fontsize=plt.rcParams["xtick.labelsize"],
            color=BLACK_TEXT,
            clip_on=False,
        )
    plot_left = min(left for left, _ in group_bounds.values())
    plot_right = max(right for _, right in group_bounds.values())
    ax.text(
        (plot_left + plot_right) / 2,
        age_axis_label_y,
        "Age (weeks)",
        transform=label_transform,
        ha="center",
        va="center",
        fontsize=plt.rcParams["axes.labelsize"],
        color=BLACK_TEXT,
        clip_on=False,
    )
    ax.set_xticks([])
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.set_ylabel("GFP x biomass\n(x1000)", labelpad=18)
    ax.set_ylim(0, ymax * 1.18 if ymax else 1)
    ax.set_xlim(black_alignment_limits(group_bounds))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    positive_ticks = [tick for tick in ax.get_yticks() if tick >= 0]
    ax.set_yticks(positive_ticks)
    n_text = n_label_from_counts(counts)
    if n_text:
        ax.text(
            1.01,
            0.97,
            n_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=annotation_fontsize(),
            color=BLACK_TEXT,
            clip_on=False,
        )
    ax.text(
        (first_leaf_x if first_leaf_x is not None else min(group_centers.values())) - 0.20,
        leaf_label_y,
        leaf_number_heading,
        transform=label_transform,
        ha="right",
        va="center",
        fontsize=plt.rcParams["xtick.labelsize"],
        color=BLACK_TEXT,
        clip_on=False,
    )
    style_axis_dark(ax)
    ax.spines["bottom"].set_visible(False)
    if positive_ticks:
        ax.spines["left"].set_bounds(0, max(positive_ticks))


def make_age_image_compilation():
    fig = plt.figure(figsize=(10.8, 4.6), constrained_layout=False)
    plot_age_image_grid(fig, fig.add_gridspec(1, 1)[0, 0])
    fig.subplots_adjust(left=0.05, right=0.99, top=0.98, bottom=0.04)
    save(fig, "age_image_compilation")
    plt.close(fig)


def make_master_figure(plant_data, leaf_data, pairwise):
    previous_rcparams = plt.rcParams.copy()
    plt.rcParams.update(MASTER_ONTOGENESIS_RCPARAMS)
    fig = plt.figure(figsize=(14.0, 15.4), constrained_layout=False)
    outer = fig.add_gridspec(
        4,
        6,
        height_ratios=[2.35, 1.22, 1.86, 1.34],
        hspace=0.78,
        wspace=0.62,
    )

    ax_a = plot_leaf_number_variance_panel(fig, outer[0, :3])

    ax_b = fig.add_subplot(outer[0, 3:])
    plot_compiled_plant_images(ax_b)

    ax_c1 = fig.add_subplot(outer[1, :3])
    plot_age_bar(
        ax_c1,
        plant_data,
        "stem_height",
        "Stem height",
        "Height per plant age",
        pairwise,
    )

    ax_c2 = fig.add_subplot(outer[1, 3:])
    plot_age_bar(
        ax_c2,
        plant_data,
        "number_of_leaves",
        "Number of leaves",
        "Number of leaves per plant age",
        pairwise,
    )

    ax_d = fig.add_subplot(outer[2, :])
    plot_leaf_gfp_grouped(ax_d, leaf_data)

    ax_e = fig.add_subplot(outer[3, 1:5])
    plot_age_bar(
        ax_e,
        plant_data,
        "total_plant_gfp_expressed",
        "Total plant GFP expressed",
        "Total plant GFP expressed per age",
        pairwise,
        y_scale=TOTAL_GFP_Y_SCALE,
    )

    fig.subplots_adjust(left=0.060, right=0.985, top=0.985, bottom=0.115)
    master_panel_label(fig, ax_a, "A", x=0.030)
    master_panel_label(fig, ax_b, "B", x=ax_b.get_position().x0 - 0.028)
    master_panel_label(fig, ax_c1, "C", x=0.030)
    master_panel_label(fig, ax_d, "D", x=0.030)
    master_panel_label(fig, ax_e, "E", x=0.030)
    save(fig, "master_figure")
    plt.close(fig)
    plt.rcParams.update(previous_rcparams)


def make_master_figure_alternate(plant_data, leaf_data, pairwise):
    previous_rcparams = plt.rcParams.copy()
    plt.rcParams.update(MASTER_ONTOGENESIS_RCPARAMS)
    _, _, _, _, age_positions, group_bounds = leaf_group_geometry(leaf_data)
    aligned_xlim = leaf_alignment_limits(group_bounds)

    fig = plt.figure(figsize=(14.0, 17.2), constrained_layout=False)
    outer = fig.add_gridspec(
        5,
        1,
        height_ratios=[2.55, 1.85, 1.15, 1.15, 1.30],
        hspace=0.66,
    )

    ax_a = plot_leaf_number_variance_panel(fig, outer[0, 0])

    ax_b = fig.add_subplot(outer[1, 0])
    plot_leaf_gfp_grouped(ax_b, leaf_data)

    ax_c1 = fig.add_subplot(outer[2, 0])
    plot_age_bar(
        ax_c1,
        plant_data,
        "stem_height",
        "Stem height",
        "Height per plant age",
        pairwise,
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.42,
    )

    ax_c2 = fig.add_subplot(outer[3, 0])
    plot_age_bar(
        ax_c2,
        plant_data,
        "number_of_leaves",
        "Number of leaves",
        "Number of leaves per plant age",
        pairwise,
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.42,
    )

    ax_e = fig.add_subplot(outer[4, 0])
    plot_age_bar(
        ax_e,
        plant_data,
        "total_plant_gfp_expressed",
        "Total plant GFP expressed",
        "Total plant GFP expressed per age",
        pairwise,
        y_scale=TOTAL_GFP_Y_SCALE,
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.42,
    )

    fig.subplots_adjust(left=0.085, right=0.975, top=0.985, bottom=0.095)
    master_panel_label(fig, ax_a, "A", x=0.038)
    master_panel_label(fig, ax_b, "B", x=0.038)
    master_panel_label(fig, ax_c1, "C", x=0.038)
    master_panel_label(fig, ax_e, "E", x=0.038)
    save(fig, "master_figure_alternate")
    plt.close(fig)
    plt.rcParams.update(previous_rcparams)


def make_master_figure_black(plant_data, leaf_data, pairwise):
    previous_rcparams = plt.rcParams.copy()
    dark_rcparams = {
        **MASTER_ONTOGENESIS_RCPARAMS,
        "font.size": 32,
        "axes.labelsize": 32,
        "xtick.labelsize": 31,
        "ytick.labelsize": 31,
        "text.color": BLACK_TEXT,
        "axes.labelcolor": BLACK_TEXT,
        "axes.edgecolor": BLACK_AXIS,
        "xtick.color": BLACK_TEXT,
        "ytick.color": BLACK_TEXT,
        "savefig.facecolor": BLACK_BG,
        "figure.facecolor": BLACK_BG,
    }
    plt.rcParams.update(dark_rcparams)
    _, _, _, _, age_positions, group_bounds = black_leaf_group_geometry(leaf_data)
    aligned_xlim = black_alignment_limits(group_bounds)

    fig = plt.figure(figsize=(17.6, 28.0), facecolor=BLACK_BG, constrained_layout=False)
    outer = fig.add_gridspec(
        7,
        1,
        height_ratios=[6.20, 3.70, 3.00, 0.52, 1.45, 1.45, 1.62],
        hspace=1.05,
    )

    ax_a = fig.add_subplot(outer[0, 0])
    ax_a.imshow(image_region_array(reordered_leaf_number_image_path(), (35, 55, 1845, 846)))
    ax_a.set_axis_off()
    ax_a.set_facecolor(BLACK_BG)
    # Five times the length of the original 1 cm bar in this panel.
    add_dark_scale_bar(ax_a, length_axes=0.085, x_right=1.075, y=0.070)

    ax_b = fig.add_subplot(outer[1, 0])
    plot_leaf_gfp_grouped_dark(ax_b, leaf_data)

    ax_plant = fig.add_subplot(outer[2, 0])
    plot_aligned_side_plants_dark(ax_plant, age_positions, aligned_xlim)
    # The first pot is 9 cm wide; this bar is 5/9 of its rendered width.
    add_dark_scale_bar(ax_plant, length_axes=0.029, x_right=0.975, y=0.070)

    ax_labels = fig.add_subplot(outer[3, 0])
    plot_black_age_labels(ax_labels, age_positions, aligned_xlim)

    ax_c1 = fig.add_subplot(outer[4, 0])
    plot_age_bar_dark(
        ax_c1,
        plant_data,
        "stem_height",
        "Stem\nheight",
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.82,
    )

    ax_c2 = fig.add_subplot(outer[5, 0])
    plot_age_bar_dark(
        ax_c2,
        plant_data,
        "number_of_leaves",
        "Number of\nleaves",
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.82,
    )

    ax_e = fig.add_subplot(outer[6, 0])
    plot_age_bar_dark(
        ax_e,
        plant_data,
        "total_plant_gfp_expressed",
        "Total GFP",
        y_scale=TOTAL_GFP_Y_SCALE,
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.82,
        bar_color="#66706A",
    )

    fig.subplots_adjust(left=0.19, right=0.945, top=0.955, bottom=0.090)
    for label, ax in [("A.", ax_a), ("B.", ax_b), ("C.", ax_plant), ("D.", ax_c1), ("E.", ax_e)]:
        bbox = ax.get_position()
        fig.text(
            0.020,
            bbox.y1,
            label,
            ha="left",
            va="top",
            fontsize=38,
            color=BLACK_TEXT,
        )
    save_dark(fig, "master_figure_black")
    plt.close(fig)
    plt.rcParams.update(previous_rcparams)


def make_master_figure_black_alternate(plant_data, leaf_data, pairwise):
    previous_rcparams = plt.rcParams.copy()
    dark_rcparams = {
        **MASTER_ONTOGENESIS_RCPARAMS,
        "font.size": 22,
        "axes.labelsize": 22,
        "xtick.labelsize": 21,
        "ytick.labelsize": 21,
        "text.color": BLACK_TEXT,
        "axes.labelcolor": BLACK_TEXT,
        "axes.edgecolor": BLACK_AXIS,
        "xtick.color": BLACK_TEXT,
        "ytick.color": BLACK_TEXT,
        "savefig.facecolor": BLACK_BG,
        "figure.facecolor": BLACK_BG,
    }
    plt.rcParams.update(dark_rcparams)
    _, _, _, _, age_positions, group_bounds = black_leaf_group_geometry(leaf_data)
    aligned_xlim = black_alignment_limits(group_bounds)

    fig = plt.figure(figsize=(17.6, 22.3), facecolor=BLACK_BG, constrained_layout=False)
    outer = fig.add_gridspec(
        6,
        1,
        height_ratios=[0.34, 2.35, 1.35, 1.42, 5.20, 2.15],
        hspace=0.78,
    )

    ax_labels = fig.add_subplot(outer[0, 0])
    plot_black_age_labels(ax_labels, age_positions, aligned_xlim)

    ax_a = fig.add_subplot(outer[1, 0])
    plot_aligned_side_plants_dark(ax_a, age_positions, aligned_xlim)

    metrics_grid = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[2, 0], wspace=0.40)
    ax_b = fig.add_subplot(metrics_grid[0, 0])
    plot_age_bar_dark(
        ax_b,
        plant_data,
        "stem_height",
        "Stem height",
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.72,
    )
    ax_b2 = fig.add_subplot(metrics_grid[0, 1])
    plot_age_bar_dark(
        ax_b2,
        plant_data,
        "number_of_leaves",
        "Number of leaves",
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.72,
    )

    ax_c = fig.add_subplot(outer[3, 0])
    plot_age_bar_dark(
        ax_c,
        plant_data,
        "total_plant_gfp_expressed",
        "Total plant GFP\nexpressed",
        y_scale=TOTAL_GFP_Y_SCALE,
        age_positions=age_positions,
        xlim=aligned_xlim,
        bar_width=0.82,
    )

    ax_d = fig.add_subplot(outer[4, 0])
    ax_d.imshow(image_region_array(reordered_leaf_number_image_path(), (118, 72, 1835, 846)))
    ax_d.set_axis_off()
    ax_d.set_facecolor(BLACK_BG)

    ax_e = fig.add_subplot(outer[5, 0])
    plot_leaf_gfp_grouped_dark(ax_e, leaf_data)

    fig.subplots_adjust(left=0.19, right=0.945, top=0.955, bottom=0.090)
    for label, ax in [("A.", ax_a), ("B.", ax_b), ("C.", ax_c), ("D.", ax_d), ("E.", ax_e)]:
        bbox = ax.get_position()
        fig.text(
            0.020,
            bbox.y1,
            label,
            ha="left",
            va="top",
            fontsize=28,
            color=BLACK_TEXT,
        )
    save_dark(fig, "master_figure_black_alternate")
    plt.close(fig)
    plt.rcParams.update(previous_rcparams)


def make_split_figures(plant_data, leaf_data):
    """Export the black master as three roomier, publication-ready figure groups."""
    previous_rcparams = plt.rcParams.copy()
    plt.rcParams.update({
        **MASTER_ONTOGENESIS_RCPARAMS,
        "font.size": 32,
        "axes.labelsize": 32,
        "xtick.labelsize": 31,
        "ytick.labelsize": 31,
        "text.color": BLACK_TEXT,
        "axes.labelcolor": BLACK_TEXT,
        "axes.edgecolor": BLACK_AXIS,
        "xtick.color": BLACK_TEXT,
        "ytick.color": BLACK_TEXT,
        "savefig.facecolor": BLACK_BG,
        "figure.facecolor": BLACK_BG,
    })
    _, _, _, _, age_positions, group_bounds = black_leaf_group_geometry(leaf_data)
    aligned_xlim = black_alignment_limits(group_bounds)

    # Standalone A: leaf-number reference and fluorescence images.
    fig_a = plt.figure(figsize=(12.0, 10.2), facecolor=BLACK_BG)
    ax_a = fig_a.add_subplot(111)
    ax_a.imshow(image_array(standalone_leaf_number_image_path(), crop=False))
    ax_a.set_axis_off()
    ax_a.set_facecolor(BLACK_BG)
    add_dark_scale_bar(ax_a, length_axes=0.072, x_right=1.055, y=0.070)
    fig_a.text(0.018, 0.965, "A.", ha="left", va="top", fontsize=38, color=BLACK_TEXT)
    fig_a.subplots_adjust(left=0.065, right=0.925, top=0.965, bottom=0.055)
    save_split_dark(fig_a, "202607_ontogenesis_split_A_leaf_number_images")
    plt.close(fig_a)

    # C + D: plant images, age labels below them, and the two plant metrics.
    fig_cd = plt.figure(figsize=(17.6, 13.4), facecolor=BLACK_BG)
    grid_cd = fig_cd.add_gridspec(
        4, 1, height_ratios=[4.2, 0.62, 1.85, 1.85], hspace=0.82
    )
    ax_c = fig_cd.add_subplot(grid_cd[0, 0])
    plot_aligned_side_plants_dark(ax_c, age_positions, aligned_xlim)
    add_dark_scale_bar(ax_c, length_axes=0.029, x_right=0.975, y=0.070)
    ax_age = fig_cd.add_subplot(grid_cd[1, 0])
    plot_black_age_labels(ax_age, age_positions, aligned_xlim)
    ax_d1 = fig_cd.add_subplot(grid_cd[2, 0])
    plot_age_bar_dark(
        ax_d1, plant_data, "stem_height", "Stem\nheight",
        age_positions=age_positions, xlim=aligned_xlim, bar_width=0.82,
    )
    ax_d2 = fig_cd.add_subplot(grid_cd[3, 0])
    plot_age_bar_dark(
        ax_d2, plant_data, "number_of_leaves", "Number of\nleaves",
        age_positions=age_positions, xlim=aligned_xlim, bar_width=0.82,
    )
    fig_cd.subplots_adjust(left=0.19, right=0.945, top=0.965, bottom=0.075)
    fig_cd.text(0.020, ax_c.get_position().y1, "C.", ha="left", va="top", fontsize=38, color=BLACK_TEXT)
    fig_cd.text(0.020, ax_d1.get_position().y1, "D.", ha="left", va="top", fontsize=38, color=BLACK_TEXT)
    save_split_dark(fig_cd, "202607_ontogenesis_split_CD_plants_and_growth")
    plt.close(fig_cd)

    # B + E: aligned age groups with extra room for B's two-tier leaf labels.
    fig_be = plt.figure(figsize=(28.0, 12.8), facecolor=BLACK_BG)
    grid_be = fig_be.add_gridspec(2, 1, height_ratios=[5.4, 2.25], hspace=1.05)
    ax_b = fig_be.add_subplot(grid_be[0, 0])
    plot_leaf_gfp_grouped_dark(
        ax_b, leaf_data, leaf_label_rotation=0, leaf_label_fontsize=29,
        leaf_labels_odd_only=True,
    )
    ax_e = fig_be.add_subplot(grid_be[1, 0])
    plot_age_bar_dark(
        ax_e, plant_data, "total_plant_gfp_expressed", "Total GFP",
        y_scale=TOTAL_GFP_Y_SCALE, age_positions=age_positions, xlim=aligned_xlim,
        bar_width=0.82, bar_color="#66706A",
    )
    fig_be.subplots_adjust(left=0.19, right=0.945, top=0.955, bottom=0.105)
    fig_be.text(0.020, ax_b.get_position().y1, "B.", ha="left", va="top", fontsize=38, color=BLACK_TEXT)
    fig_be.text(0.020, ax_e.get_position().y1, "E.", ha="left", va="top", fontsize=38, color=BLACK_TEXT)
    save_split_dark(fig_be, "202607_ontogenesis_split_BE_gfp_by_age")
    plt.close(fig_be)
    plt.rcParams.update(previous_rcparams)


def make_reordered_wide_master(plant_data, leaf_data):
    """Wide C-B-E-D master with the revised PDF-derived A panel at right."""
    previous_rcparams = plt.rcParams.copy()
    plt.rcParams.update({
        **MASTER_ONTOGENESIS_RCPARAMS,
        "font.size": 30, "axes.labelsize": 30,
        "xtick.labelsize": 29, "ytick.labelsize": 29,
        "text.color": BLACK_TEXT, "axes.labelcolor": BLACK_TEXT,
        "axes.edgecolor": BLACK_AXIS, "xtick.color": BLACK_TEXT,
        "ytick.color": BLACK_TEXT, "savefig.facecolor": BLACK_BG,
        "figure.facecolor": BLACK_BG,
    })
    _, _, _, _, age_positions, group_bounds = black_leaf_group_geometry(leaf_data)
    aligned_xlim = black_alignment_limits(group_bounds)

    fig = plt.figure(figsize=(46.0, 34.0), facecolor=BLACK_BG)
    outer = fig.add_gridspec(
        4, 2, width_ratios=[2.35, 1.50],
        height_ratios=[5.50, 4.15, 3.45, 2.55],
        hspace=0.86, wspace=0.18,
    )

    c_grid = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[0, 0], height_ratios=[4.0, 0.62], hspace=0.08
    )
    ax_c = fig.add_subplot(c_grid[0, 0])
    plot_aligned_side_plants_dark(ax_c, age_positions, aligned_xlim, image_zoom=0.50)
    add_dark_scale_bar(ax_c, length_axes=0.064, x_right=1.075, y=0.015)
    ax_c_age = fig.add_subplot(c_grid[1, 0])
    plot_black_age_labels(ax_c_age, age_positions, aligned_xlim)

    d_grid = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1, 0], hspace=0.72)
    ax_d1 = fig.add_subplot(d_grid[0, 0])
    plot_age_bar_dark(
        ax_d1, plant_data, "stem_height", "Stem\nheight",
        age_positions=age_positions, xlim=aligned_xlim, bar_width=0.82,
        exact_y_tick_count=3,
    )
    ax_d2 = fig.add_subplot(d_grid[1, 0])
    plot_age_bar_dark(
        ax_d2, plant_data, "number_of_leaves", "Number of\nleaves",
        age_positions=age_positions, xlim=aligned_xlim, bar_width=0.82,
        exact_y_tick_count=3,
    )

    ax_b = fig.add_subplot(outer[2, 0])
    plot_leaf_gfp_grouped_dark(
        ax_b, leaf_data, leaf_label_fontsize=27,
        leaf_labels_odd_only=True, leaf_label_y_override=-0.125,
        leaf_number_heading="Leaf\nno:",
    )

    ax_e = fig.add_subplot(outer[3, 0])
    plot_age_bar_dark(
        ax_e, plant_data, "total_plant_gfp_expressed", "Total GFP\n",
        y_scale=TOTAL_GFP_Y_SCALE, age_positions=age_positions,
        xlim=aligned_xlim, bar_width=0.82, bar_color="#66706A",
        exact_y_tick_count=3,
    )

    ax_a = fig.add_subplot(outer[:, 1])
    ax_a.imshow(image_array(vertical_leaf_age_pdf_panel_path(), crop=False))
    ax_a.set_axis_off()
    ax_a.set_facecolor(BLACK_BG)
    add_dark_scale_bar(ax_a, length_axes=0.245, x_right=0.965, y=0.020)

    fig.subplots_adjust(left=0.095, right=0.975, top=0.975, bottom=0.055)
    for label, ax in [("C.", ax_c), ("B.", ax_b), ("E.", ax_e)]:
        fig.text(0.012, ax.get_position().y1, label, ha="left", va="top", fontsize=38, color=BLACK_TEXT)
    fig.text(
        ax_a.get_position().x0 - 0.025, ax_a.get_position().y1,
        "A.", ha="left", va="top", fontsize=38, color=BLACK_TEXT,
    )
    save_split_dark(fig, "202607_ontogenesis_master_reordered_vertical_A")
    plt.close(fig)
    plt.rcParams.update(previous_rcparams)


def leafage_pptx_media_image(media_name):
    pptx_path = SPLIT_DIR / "leafAge.pptx"
    with zipfile.ZipFile(pptx_path) as archive:
        payload = archive.read(f"ppt/media/{media_name}")
    return np.asarray(Image.open(BytesIO(payload)).convert("RGBA"))


def make_leafage_slide2_final(plant_data, leaf_data):
    """Final 2x2 composite following slide 2 of leafAge.pptx."""
    previous_rcparams = plt.rcParams.copy()
    font_size = 30
    plt.rcParams.update({
        **MASTER_ONTOGENESIS_RCPARAMS,
        "font.size": font_size, "axes.labelsize": font_size,
        "xtick.labelsize": font_size, "ytick.labelsize": font_size,
        "legend.fontsize": font_size, "text.color": BLACK_TEXT,
        "axes.labelcolor": BLACK_TEXT, "axes.edgecolor": BLACK_AXIS,
        "xtick.color": BLACK_TEXT, "ytick.color": BLACK_TEXT,
        "savefig.facecolor": BLACK_BG, "figure.facecolor": BLACK_BG,
    })
    _, _, _, _, age_positions, group_bounds = black_leaf_group_geometry(leaf_data)
    aligned_xlim = black_alignment_limits(group_bounds)

    fig = plt.figure(figsize=(31.0, 12.0), facecolor=BLACK_BG)
    outer = fig.add_gridspec(
        2, 2, width_ratios=[1.0, 1.50], height_ratios=[1.0, 1.0],
        wspace=0.29, hspace=0.82,
    )

    # A: side-view plants, age labels, and a lower-right 5 cm scale bar.
    a_grid = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[0, 0], height_ratios=[4.55, 0.72], hspace=0.02
    )
    ax_a = fig.add_subplot(a_grid[0, 0])
    # Panel A follows the independent, evenly spaced layout from slide 2 rather
    # than inheriting the variable-width geometry of panel C.
    a_age_positions = {age: i * 1.95 for i, age in enumerate(AGE_ORDER)}
    a_xlim = (-0.68, a_age_positions[8] + 0.68)
    plot_aligned_side_plants_dark(
        ax_a, a_age_positions, a_xlim, image_zoom=0.220,
        y_by_age={4: -0.025},
    )
    ax_a_age = fig.add_subplot(a_grid[1, 0])
    plot_black_age_labels(
        ax_a_age, a_age_positions, a_xlim, fontsize=font_size,
        suffix_by_age={5: "*"},
    )
    # This is twice the calibrated 5 cm bar length used previously.
    # Keep it low in the label band so the enlarged plants remain unobstructed.
    add_dark_scale_bar(
        # The rendered 6-week pot is 9 cm wide; 0.050 axes units is 10/9
        # of that displayed pot width.
        ax_a_age, length_axes=0.050, label="10 cm", x_right=0.955, y=-0.800,
        fontsize=font_size, linewidth=3.2,
    )

    # B: slide-2 top-view plant plus the two four-leaf strips.
    ax_b = fig.add_subplot(outer[1, 0])
    ax_b.set_axis_off(); ax_b.set_facecolor(BLACK_BG)
    plant_ax = ax_b.inset_axes([0.00, 0.10, 0.40, 0.84])
    plant_ax.imshow(leafage_pptx_media_image("image9.png")); plant_ax.set_axis_off()
    plant_ax.text(0.04, 0.94, "*", transform=plant_ax.transAxes,
                  ha="left", va="top", fontsize=font_size,
                  fontweight="bold", color=BLACK_TEXT)
    leaf58_ax = ax_b.inset_axes([0.42, 0.54, 0.56, 0.40])
    leaf58_ax.imshow(leafage_pptx_media_image("image10.png")); leaf58_ax.set_axis_off()
    leaf14_ax = ax_b.inset_axes([0.42, 0.10, 0.56, 0.40])
    leaf14_ax.imshow(leafage_pptx_media_image("image11.png")); leaf14_ax.set_axis_off()

    plant_label_positions = {
        1: (0.70, 0.26), 2: (0.36, 0.85), 3: (0.24, 0.21), 4: (0.83, 0.49),
        5: (0.28, 0.57), 6: (0.54, 0.43), 7: (0.47, 0.57), 8: (0.62, 0.36),
    }
    for number, (x, y) in plant_label_positions.items():
        plant_ax.text(x, y, str(number), transform=plant_ax.transAxes,
                      ha="center", va="center", fontsize=font_size, color=BLACK_TEXT)
    for strip_ax, numbers in [(leaf14_ax, [1, 2, 3, 4]), (leaf58_ax, [5, 6, 7, 8])]:
        for x, number in zip([0.125, 0.375, 0.625, 0.875], numbers):
            strip_ax.text(x, 0.88, str(number), transform=strip_ax.transAxes,
                          ha="center", va="center", fontsize=font_size, color=BLACK_TEXT)
    add_dark_scale_bar(
        ax_b, length_axes=0.040, x_right=0.975, y=0.055,
        fontsize=font_size, linewidth=3.2,
    )

    # C occupies the full height of A.
    ax_c = fig.add_subplot(outer[0, 1])
    plot_leaf_gfp_grouped_dark(
        ax_c, leaf_data, leaf_label_fontsize=font_size,
        leaf_label_every=3, leaf_label_y_override=-0.105,
        leaf_number_heading="Leaf\nno:",
    )

    # D and E stack to exactly the same overall height as B.
    de_grid = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1, 1], hspace=1.18)
    ax_d = fig.add_subplot(de_grid[0, 0])
    plot_age_bar_dark(
        ax_d, plant_data, "total_plant_gfp_expressed", "Total GFP\n",
        y_scale=TOTAL_GFP_Y_SCALE, age_positions=age_positions, xlim=aligned_xlim,
        bar_width=0.82, bar_color="#66706A", exact_y_tick_count=3,
    )
    ax_e = fig.add_subplot(de_grid[1, 0])
    plot_age_bar_dark(
        ax_e, plant_data, "number_of_leaves", "Number of\nleaves",
        age_positions=age_positions, xlim=aligned_xlim, bar_width=0.82,
        exact_y_tick_count=3,
    )

    fig.subplots_adjust(left=0.065, right=0.985, top=0.945, bottom=0.110)
    # Panel letters sit just outside the panels.  Their shared coordinates keep
    # A/C and B/D level, and make A/B and C/D/E true vertical columns.
    fig.canvas.draw()
    left_label_x = ax_a.get_position().x0 - 0.018
    # Clear the complete y-axis title, rather than placing letters in the plots.
    right_label_x = ax_c.get_position().x0 - 0.082
    top_label_y = max(ax_a.get_position().y1, ax_c.get_position().y1) + 0.010
    lower_label_y = max(ax_b.get_position().y1, ax_d.get_position().y1) + 0.010
    panel_label_positions = {
        "A": (left_label_x, top_label_y),
        "B": (left_label_x, lower_label_y),
        "C": (right_label_x, top_label_y),
        "D": (right_label_x, lower_label_y),
        "E": (right_label_x, ax_e.get_position().y1 + 0.010),
    }
    for label, (x, y) in panel_label_positions.items():
        fig.text(x, y, label, ha="left", va="bottom", fontsize=font_size,
                 fontweight="bold", color=BLACK_TEXT, zorder=30)
    save_split_dark(fig, "202607_ontogenesis_leafage_slide2_final")
    plt.close(fig)
    plt.rcParams.update(previous_rcparams)


def make_figures(plant_data, leaf_data, pairwise):
    fig, ax = plt.subplots(figsize=(4.4, 3.9), constrained_layout=False)
    plot_age_bar(
        ax,
        plant_data,
        "stem_height",
        "Stem height",
        "Height per plant age",
        pairwise,
    )
    fig.subplots_adjust(left=0.17, right=0.97, top=0.84, bottom=0.33)
    save(fig, "height_per_age")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 3.9), constrained_layout=False)
    plot_age_bar(
        ax,
        plant_data,
        "number_of_leaves",
        "Number of leaves",
        "Number of leaves per plant age",
        pairwise,
    )
    fig.subplots_adjust(left=0.17, right=0.97, top=0.84, bottom=0.33)
    save(fig, "number_of_leaves_per_age")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.9), constrained_layout=False)
    plot_age_bar(
        axes[0],
        plant_data,
        "stem_height",
        "Stem height",
        "Height per plant age",
        pairwise,
    )
    plot_age_bar(
        axes[1],
        plant_data,
        "number_of_leaves",
        "Number of leaves",
        "Number of leaves per plant age",
        pairwise,
    )
    fig.subplots_adjust(left=0.10, right=0.98, top=0.82, bottom=0.30, wspace=0.34)
    save(fig, "plant_metrics_by_age")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.8, 5.0), constrained_layout=False)
    plot_leaf_gfp_grouped(ax, leaf_data)
    fig.subplots_adjust(left=0.08, right=0.96, top=0.90, bottom=0.42)
    save(fig, "gfp_fluorescence_by_leaf_number_and_age")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 3.9), constrained_layout=False)
    plot_age_bar(
        ax,
        plant_data,
        "total_plant_gfp_expressed",
        "Total plant GFP expressed",
        "Total plant GFP expressed per age",
        pairwise,
        y_scale=TOTAL_GFP_Y_SCALE,
    )
    fig.subplots_adjust(left=0.17, right=0.97, top=0.82, bottom=0.33)
    save(fig, "total_plant_gfp_expressed_by_age")
    plt.close(fig)

    make_master_figure(plant_data, leaf_data, pairwise)
    make_master_figure_alternate(plant_data, leaf_data, pairwise)
    make_master_figure_black(plant_data, leaf_data, pairwise)
    make_master_figure_black_alternate(plant_data, leaf_data, pairwise)
    make_split_figures(plant_data, leaf_data)
    make_reordered_wide_master(plant_data, leaf_data)
    make_leafage_slide2_final(plant_data, leaf_data)


# %%
clear_previous_outputs()
raw = load_raw_data()
plant_metrics = make_plant_metrics(raw)
leaf_data = make_leaf_data(raw)

plant_screened, plant_outliers = flag_outliers(
    plant_metrics,
    ["stem_height", "number_of_leaves", "total_plant_gfp_expressed"],
    ["plant_id", "plant_age", "replicate"],
)
leaf_screened, leaf_outliers = flag_outliers(
    leaf_data,
    ["gfp_expression", "gfp_expressed_x_biomass"],
    ["sample_id", "plant_id", "plant_age", "replicate", "leaf_number"],
)
outlier_log = pd.concat(
    [
        plant_outliers.assign(data_level="plant"),
        leaf_outliers.assign(data_level="leaf"),
    ],
    ignore_index=True,
)

raw.to_excel(LOG_DIR / "ontogenesis_raw_clean_header.xlsx", index=False)
plant_metrics.to_excel(LOG_DIR / "ontogenesis_plant_metrics.xlsx", index=False)
leaf_data.to_excel(LOG_DIR / "ontogenesis_leaf_data.xlsx", index=False)
plant_screened.to_excel(LOG_DIR / "ontogenesis_plant_metrics_outlier_screened.xlsx", index=False)
leaf_screened.to_excel(LOG_DIR / "ontogenesis_leaf_data_outlier_screened.xlsx", index=False)
outlier_log.to_excel(LOG_DIR / "ontogenesis_outlier_log.xlsx", index=False)

anova, assumptions, pairwise, leaf_model, leaf_pairwise = run_all_stats(plant_metrics, leaf_data)
make_figures(plant_metrics, leaf_data, pairwise)

print(f"Loaded {len(raw)} leaf rows from {SOURCE_FILE.name}.")
print(f"Summarised {len(plant_metrics)} plant replicates across ages {AGE_ORDER}.")
print(f"Flagged {len(outlier_log)} potential outlier rows; removed 0.")
print(f"Figures saved to: {SAVE_DIR}")
print(f"Logs saved to: {LOG_DIR}")
