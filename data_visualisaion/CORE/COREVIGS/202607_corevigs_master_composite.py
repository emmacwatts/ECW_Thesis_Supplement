"""Master composite of the COREVIGS GUS and PVX-GFP experiments."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import os

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib-cache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import MultipleLocator, StrMethodFormatter


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gus = load_module("corevigs_gus", HERE / "202607_corevigs_gus_5dpi.py")
pvx = load_module(
    "corevigs_pvx", HERE / "PVX_GFP" / "202607_corevigs_pvx_gfp_batch2.py"
)

OUT = HERE / "202607_corevigs_master_composite_figures"
OUT.mkdir(parents=True, exist_ok=True)


def make_master():
    gus_data = gus.load_data()
    gus_summary, gus_test = gus.analyse(gus_data)
    pvx_data = pvx.load_data()
    pvx_summary, pvx_tests = pvx.analyse(pvx_data)

    fig = plt.figure(figsize=(13.6, 12.4), facecolor="white")
    gs = fig.add_gridspec(
        2, 2, width_ratios=[0.88, 1.12], height_ratios=[1, 1],
        left=0.085, right=0.985, bottom=0.075, top=0.945,
        wspace=0.12, hspace=0.25,
    )
    axes = [
        fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]),
    ]

    pvx.grouped_plot(axes[0], pvx_data, pvx_summary, pvx_tests)
    pvx.photo_panel(axes[1])
    gus.add_plot(axes[2], gus_data, gus_summary, gus_test)
    gus.add_photo_panel(axes[3])

    axes[2].yaxis.set_major_locator(MultipleLocator(1))
    axes[2].yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))

    # The grouped PVX-GFP bars in panel A retain their genotype colour key;
    # panel C is labelled directly on its x-axis.
    panel_legend = axes[0].get_legend()
    if panel_legend is not None:
        for legend_text in panel_legend.get_texts():
            legend_text.set_fontsize(15)

    # Increase every existing text element by two points for the master figure.
    for text in fig.findobj(match=Text):
        text.set_fontsize(text.get_fontsize() + 2)

    for ax, letter in zip(axes, "ABCD"):
        box = ax.get_position()
        fig.text(
            box.x0 - 0.012, box.y1 + 0.035, letter,
            ha="left", va="top", fontsize=22, fontweight="bold",
            fontfamily="Helvetica Neue", color=gus.TEXT_COLOR,
        )

    stem = "202607_corevigs_gus_pvx_gfp_master_composite"
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{extension}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    make_master()
