"""Assemble the ELISA dilution-response and western-verification panels."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont, ImageOps, JpegImagePlugin  # noqa: F401
import numpy as np


HERE = Path(__file__).resolve().parent
PLOT = (
    HERE / "202607_cova_elisa_protocol_verification_figures/"
    "202607_cova_elisa_protocol_verification_dilution_response.png"
)
WESTERN = HERE / "Western_result.png"
OUTPUT = HERE / "202607_cova_elisa_protocol_verification_composite_figures"
OUTPUT.mkdir(parents=True, exist_ok=True)

CANVAS = (3200, 3900)
WHITE = (255, 255, 255)
TEXT_COLOR = (74, 74, 74)  # #4A4A4A, matching the plot style used in A
FONT = ImageFont.truetype("/System/Library/Fonts/HelveticaNeue.ttc", 76, index=0)
WESTERN_FONT = ImageFont.truetype("/System/Library/Fonts/HelveticaNeue.ttc", 74, index=0)
WESTERN_ITALIC_FONT = ImageFont.truetype("/System/Library/Fonts/HelveticaNeue.ttc", 74, index=2)


def open_rgb(path):
    image = Image.open(path)
    if image.mode == "RGBA":
        base = Image.new("RGBA", image.size, WHITE + (255,))
        image = Image.alpha_composite(base, image)
    return image.convert("RGB")


def trim_white(image, threshold=250, padding=10):
    mask = Image.new("L", image.size)
    mask.putdata([255 if min(pixel) < threshold else 0 for pixel in image.getdata()])
    bbox = mask.getbbox()
    if bbox is None:
        return image
    left, top, right, bottom = bbox
    return image.crop((max(0, left - padding), max(0, top - padding),
                       min(image.width, right + padding), min(image.height, bottom + padding)))


def place(canvas, image, box):
    left, top, right, bottom = box
    fitted = ImageOps.contain(image, (right - left, bottom - top), Image.Resampling.LANCZOS)
    x = left + (right - left - fitted.width) // 2
    y = top + (bottom - top - fitted.height) // 2
    canvas.paste(fitted, (x, y))


def text_center(draw, xy, text, font, fill=TEXT_COLOR):
    bbox = draw.textbbox((0, 0), text, font=font)
    draw.text((xy[0] - (bbox[2] - bbox[0]) / 2,
               xy[1] - (bbox[3] - bbox[1]) / 2), text, font=font, fill=fill)


def rebuild_western(source):
    """Retain the data strips but redraw all labels at one consistent size."""
    panel = Image.new("RGB", (2920, 1050), WHITE)
    draw = ImageDraw.Draw(panel)
    # Exact interiors of the two black boxes in the supplied 1226 x 356 image.
    western_strip = source.crop((160, 162, 1194, 212))
    loading_strip = source.crop((160, 266, 1194, 316))
    western_strip = western_strip.resize((2500, 120), Image.Resampling.LANCZOS)
    loading_strip = loading_strip.resize((2500, 120), Image.Resampling.LANCZOS)
    panel.paste(western_strip, (340, 300))
    panel.paste(loading_strip, (340, 700))
    draw.rectangle((340, 300, 2840, 420), outline=(0, 0, 0), width=5)
    draw.rectangle((340, 700, 2840, 820), outline=(0, 0, 0), width=5)

    lane_centers = np.linspace(340 + 2500 / 24, 2840 - 2500 / 24, 12)
    group_centers = [lane_centers[index:index + 3].mean() for index in range(0, 12, 3)]
    for center, group in zip(group_centers, ["EV", "LC", "HC", "LCHC"]):
        text_center(draw, (center, 80), group, WESTERN_FONT)
    for center, dilution in zip(lane_centers, ["1:10", "1:100", "1:1000"] * 4):
        text_center(draw, (center, 220), dilution, WESTERN_FONT)
    text_center(draw, (170, 760), "55 kDa", WESTERN_FONT)
    draw.text((2840, 455), "α-SARS-CoV-2", font=WESTERN_ITALIC_FONT,
              fill=TEXT_COLOR, anchor="ra")
    draw.text((2840, 855), "Coomassie", font=WESTERN_ITALIC_FONT,
              fill=TEXT_COLOR, anchor="ra")
    return panel


def main():
    canvas = Image.new("RGB", CANVAS, WHITE)
    draw = ImageDraw.Draw(canvas)
    plot = trim_white(open_rgb(PLOT))
    western = rebuild_western(open_rgb(WESTERN))

    # A: ELISA dilution-response plot. B: western verification beneath it.
    place(canvas, plot, (150, 120, 3070, 2530))
    place(canvas, western, (150, 2750, 3070, 3820))
    draw.text((30, 25), "A", font=FONT, fill=(0, 0, 0))
    draw.text((30, 2600), "B", font=FONT, fill=(0, 0, 0))

    stem = OUTPUT / "202607_cova_elisa_protocol_verification_composite"
    canvas.save(stem.with_suffix(".png"), dpi=(300, 300), optimize=True)
    canvas.save(stem.with_suffix(".tiff"), dpi=(300, 300), compression="tiff_lzw")
    canvas.save(stem.with_suffix(".pdf"), resolution=300)
    print(f"Saved composite to {OUTPUT}")


if __name__ == "__main__":
    main()
