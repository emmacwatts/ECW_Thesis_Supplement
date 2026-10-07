"""Build the relabelled PRp27 plant panel with a calibrated 10 cm scale bar."""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "Images" / "prp27PlantPhotos_relabelled.png"
OUTPUT = HERE / "Images" / "prp27PlantPhotos_relabelled_10cm_scale.png"

# The horizontal side-to-side width (not the top-to-tip extent) of PRp27#1's
# bottom forward-facing leaf spans approximately 170 px in the source image.
# The requested 10 cm reference bar is drawn to that same displayed width.
SCALE_BAR_PX = 170
CANVAS_EXTENSION_PX = 130
RIGHT_MARGIN_PX = 55
BAR_Y_FROM_BOTTOM_PX = 34
BAR_HEIGHT_PX = 7


def load_font(size):
    for path in (
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ):
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


source = Image.open(SOURCE).convert("RGB")
canvas = Image.new("RGB", (source.width, source.height + CANVAS_EXTENSION_PX), "black")
canvas.paste(source, (0, 0))

draw = ImageDraw.Draw(canvas)
x_right = canvas.width - RIGHT_MARGIN_PX
x_left = x_right - SCALE_BAR_PX
y = canvas.height - BAR_Y_FROM_BOTTOM_PX
draw.rectangle((x_left, y, x_right, y + BAR_HEIGHT_PX), fill="white")

label = "10 cm"
font = load_font(34)
label_box = draw.textbbox((0, 0), label, font=font)
label_width = label_box[2] - label_box[0]
draw.text(
    ((x_left + x_right - label_width) / 2, y - 46),
    label,
    fill="white",
    font=font,
)

canvas.save(OUTPUT)
print(f"Saved {OUTPUT} ({canvas.width} x {canvas.height} px; bar={SCALE_BAR_PX} px)")
