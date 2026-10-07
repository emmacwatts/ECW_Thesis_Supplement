"""Create the CSPR TRV2::GUS and TRV2::CORE plant-photo panel."""

from pathlib import Path
import os

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))

import matplotlib.pyplot as plt
from PIL import Image

IMAGES = [
    (HERE / "7046_gus.png", "TRV2::GUS"),
    (HERE / "704_6_trv2core.png", "TRV2::CORE"),
]
OUT = HERE / "202607_cspr_trv2_photo_panel_figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "savefig.dpi": 300,
})


def trim_to_visible(image):
    image = image.convert("RGBA")
    bounds = image.getchannel("A").getbbox()
    return image.crop(bounds) if bounds else image


def make_panel():
    fig, ax = plt.subplots(figsize=(9.6, 6.6), facecolor="black")
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_facecolor("black")
    ax.set_axis_off()

    image_axes = [
        ax.inset_axes([0.035, 0.10, 0.44, 0.73]),
        ax.inset_axes([0.555, 0.15, 0.38, 0.63]),
    ]
    for image_ax, (path, _) in zip(image_axes, IMAGES):
        image_ax.set_facecolor("black")
        image_ax.imshow(trim_to_visible(Image.open(path)))
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
        "fontstyle": "italic",
    }
    ax.text(0.5, 0.965, "cspr-706-4-1", **label_style)
    ax.text(0.255, 0.865, IMAGES[0][1], **label_style)
    ax.text(0.745, 0.865, IMAGES[1][1], **label_style)

    stem = "202607_cspr_trv2_gus_vs_trv2_core_photo_panel"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            OUT / f"{stem}.{extension}",
            facecolor="black", bbox_inches="tight", pad_inches=0,
        )
    plt.close(fig)


if __name__ == "__main__":
    make_panel()
