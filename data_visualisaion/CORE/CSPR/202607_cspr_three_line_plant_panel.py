"""Compose Cas9 WT and two CSPR lines as a black-background phenotype panel."""

from pathlib import Path
import os
import subprocess

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))

import matplotlib.pyplot as plt
from PIL import Image

SOURCES = [
    (HERE / "wt.HEIC", "Cas9 WT", False),
    (HERE / "cspr311.HEIC", "cspr-3-1-1", True),
    (HERE / "cspr7064.HEIC", "cspr-706-4-1", True),
]
OUT = HERE / "202607_cspr_three_line_plant_panel_figures"
CONVERTED = OUT / "converted_sources"
OUT.mkdir(parents=True, exist_ok=True)
CONVERTED.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "savefig.dpi": 300,
})


def load_trimmed_heic(path):
    converted = CONVERTED / f"{path.stem}.png"
    if not converted.exists() or converted.stat().st_mtime < path.stat().st_mtime:
        subprocess.run(
            ["heif-convert", str(path), str(converted)],
            check=True, capture_output=True, text=True,
        )
    image = Image.open(converted).convert("RGBA")
    bounds = image.getchannel("A").getbbox()
    return image.crop(bounds) if bounds else image


def make_panel():
    fig, ax = plt.subplots(figsize=(12.4, 6.4), facecolor="black")
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_facecolor("black")
    ax.set_axis_off()

    centres = [0.18, 0.50, 0.82]
    width = 0.29
    # Equal-height boxes and a shared lower baseline preserve relative plant stature.
    for centre, (path, _, _) in zip(centres, SOURCES):
        image_ax = ax.inset_axes([centre - width / 2, 0.035, width, 0.78])
        image_ax.set_facecolor("black")
        image_ax.imshow(load_trimmed_heic(path))
        image_ax.set_anchor("S")
        image_ax.set_axis_off()

    label_style = {
        "transform": ax.transAxes,
        "ha": "center",
        "va": "center",
        "color": "white",
        "fontsize": 20,
        "fontfamily": "Helvetica Neue",
        "fontweight": 300,
    }
    for centre, (_, label, italic) in zip(centres, SOURCES):
        ax.text(
            centre, 0.91, label,
            fontstyle="italic" if italic else "normal",
            **label_style,
        )

    stem = "202607_cspr_cas9_wt_cspr_3_1_1_cspr_706_4_1_plant_panel"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            OUT / f"{stem}.{extension}",
            facecolor="black", bbox_inches="tight", pad_inches=0,
        )
    plt.close(fig)


if __name__ == "__main__":
    make_panel()
