"""Stack the BION workflow and RT-PCR results as aligned panels A and B."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
SCHEMATIC = HERE / "bion treatment.png"
RESULTS = HERE / "image_npr1_labels_corrected.png"
OUTPUT = HERE / "bion_treatment_and_results.png"

LEFT_MARGIN = 180
RIGHT_MARGIN = 80
TOP_MARGIN = 105
LABEL_GAP = 26
PANEL_GAP = 145
B_LEFT_SHIFT = 130


def proportional_width(image, target_width):
    target_height = round(image.height * target_width / image.width)
    return image.resize((target_width, target_height), Image.Resampling.LANCZOS)


def main():
    schematic = Image.open(SCHEMATIC).convert("RGB")
    results = Image.open(RESULTS).convert("RGB")
    # Remove the one-pixel export artefacts along the right and bottom edges.
    results = results.crop((0, 0, results.width - 2, results.height - 2))

    panel_width = schematic.width
    results = proportional_width(results, panel_width)
    canvas_width = LEFT_MARGIN + panel_width + RIGHT_MARGIN
    a_top = TOP_MARGIN
    b_label_y = a_top + schematic.height + PANEL_GAP
    b_top = b_label_y + TOP_MARGIN
    canvas_height = b_top + results.height + 80

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    canvas.paste(schematic, (LEFT_MARGIN, a_top))
    canvas.paste(results, (LEFT_MARGIN - B_LEFT_SHIFT, b_top))

    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 76)
    label_x = 42
    draw.text((label_x, a_top - TOP_MARGIN + LABEL_GAP), "A", font=font, fill="#3E3E3E")
    draw.text((label_x, b_label_y + LABEL_GAP), "B", font=font, fill="#3E3E3E")

    canvas.save(OUTPUT, dpi=(300, 300))
    print(OUTPUT)


if __name__ == "__main__":
    main()
