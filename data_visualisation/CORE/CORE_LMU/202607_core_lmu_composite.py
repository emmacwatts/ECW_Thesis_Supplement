"""Assemble the CORE leaf-position composite figure without resampling panel content."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont, ImageOps, JpegImagePlugin  # noqa: F401


HERE = Path(__file__).resolve().parent
ROS = HERE / "LeafAge_ROS/202607_core_leaf_age_ros_figures/202607_core_leaf_age_ros_all_elicitors_leaf_position.png"
SAMPLING = HERE / "Leaf sampling simplified v3.png"
LMU = HERE / "202607_core_lmu_figures/202607_core_lmu_leaf_position_grouped.png"
OUTPUT = HERE / "202607_core_lmu_composite_figures"
OUTPUT.mkdir(parents=True, exist_ok=True)

CANVAS_SIZE = (6000, 4900)
WHITE = (255, 255, 255)
FONT_PATH = "/System/Library/Fonts/HelveticaNeue.ttc"
LABEL_FONT = ImageFont.truetype(FONT_PATH, 108, index=1)


def open_rgb(path):
    image = Image.open(path)
    if image.mode == "RGBA":
        background = Image.new("RGBA", image.size, WHITE + (255,))
        image = Image.alpha_composite(background, image).convert("RGB")
    else:
        image = image.convert("RGB")
    return image


def trim_white(image, threshold=249, padding=12):
    """Remove only blank outer margins; preserve white areas inside each panel."""
    mask = Image.new("L", image.size)
    mask.putdata([
        255 if min(pixel) < threshold else 0
        for pixel in image.getdata()
    ])
    bbox = mask.getbbox()
    if bbox is None:
        return image
    left, top, right, bottom = bbox
    return image.crop((max(0, left - padding), max(0, top - padding),
                       min(image.width, right + padding), min(image.height, bottom + padding)))


def place_contained(canvas, image, box, align=(0.5, 0.5)):
    left, top, right, bottom = box
    target = (right - left, bottom - top)
    fitted = ImageOps.contain(image, target, Image.Resampling.LANCZOS)
    x = left + round((target[0] - fitted.width) * align[0])
    y = top + round((target[1] - fitted.height) * align[1])
    canvas.paste(fitted, (x, y))
    return x, y, fitted.width, fitted.height


def main():
    canvas = Image.new("RGB", CANVAS_SIZE, WHITE)
    draw = ImageDraw.Draw(canvas)

    ros = trim_white(open_rgb(ROS))
    sampling = trim_white(open_rgb(SAMPLING), threshold=252, padding=8)
    lmu = trim_white(open_rgb(LMU))

    # Upper row: A (sampling diagram) at left and B (ROS) at right.
    place_contained(canvas, sampling, (145, 150, 1775, 2250), align=(0.5, 0.5))
    place_contained(canvas, ros, (1910, 150, 5860, 2250), align=(1.0, 0.5))

    # Lower row: C spans the composite width while retaining its aspect ratio.
    place_contained(canvas, lmu, (145, 2520, 5860, 4740), align=(0.5, 0.5))

    # Labels sit in the outer white gutter rather than over panel content.
    draw.text((25, 38), "A", font=LABEL_FONT, fill=(0, 0, 0))
    draw.text((1790, 38), "B", font=LABEL_FONT, fill=(0, 0, 0))
    draw.text((25, 2385), "C", font=LABEL_FONT, fill=(0, 0, 0))

    stem = OUTPUT / "202607_core_lmu_leaf_age_composite"
    canvas.save(stem.with_suffix(".png"), dpi=(300, 300), optimize=True)
    canvas.save(stem.with_suffix(".tiff"), dpi=(300, 300), compression="tiff_lzw")
    canvas.save(stem.with_suffix(".pdf"), resolution=300)
    print(f"Saved composite to {OUTPUT}")


if __name__ == "__main__":
    main()
