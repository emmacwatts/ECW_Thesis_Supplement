"""Combine the CSPR phenotype photographs and immunity model as panels A and B."""

from pathlib import Path
import os

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))

import matplotlib.pyplot as plt
from PIL import Image

PHOTO = (
    HERE / "202607_cspr_trv2_photo_panel_figures"
    / "202607_cspr_trv2_gus_vs_trv2_core_photo_panel.png"
)
MODEL = HERE / "202607_core_cspr_immunity_model.png"
OUT = HERE / "202607_cspr_photo_and_model_composite_figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "savefig.dpi": 300,
})


def make_figure():
    fig = plt.figure(figsize=(14.2, 6.8), facecolor="white")
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1, 1],
        left=0.045, right=0.985, bottom=0.075, top=0.90, wspace=0.055,
    )
    axes = [fig.add_subplot(gs[0, index]) for index in range(2)]

    for ax, path, background in zip(axes, [PHOTO, MODEL], ["black", "white"]):
        ax.set_facecolor(background)
        ax.imshow(Image.open(path).convert("RGB"), aspect="auto")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Equal axes from the gridspec give the two visual panels exactly equal height.
    label_y = 0.955
    for ax, letter in zip(axes, ["A", "B"]):
        box = ax.get_position()
        fig.text(
            box.x0 - 0.012, label_y, letter,
            ha="left", va="top", fontsize=21, fontweight="bold",
            fontfamily="Helvetica Neue", color="#4A4D4F",
        )

    stem = "202607_cspr_trv2_photo_and_immunity_model"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            OUT / f"{stem}.{extension}",
            bbox_inches="tight", facecolor="white", pad_inches=0.08,
        )
    plt.close(fig)


if __name__ == "__main__":
    make_figure()
