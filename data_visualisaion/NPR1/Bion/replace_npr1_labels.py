"""Replace only the two NPR1 plant labels in the supplied Bion figure."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
source = HERE / "image.png"
output = HERE / "image_npr1_labels_corrected.png"

image = Image.open(source).convert("RGBA")
draw = ImageDraw.Draw(image)
font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Italic.ttf", 18)

# These rectangles contain only the two original plant-type labels.
draw.rectangle((274, 250, 378, 279), fill="white")
draw.rectangle((420, 250, 526, 279), fill="white")

for centre_x, label in [(326, "npr1-1"), (473, "npr1-2")]:
    box = draw.textbbox((0, 0), label, font=font)
    width = box[2] - box[0]
    draw.text((centre_x - width / 2, 254), label, font=font, fill="black")

image.save(output)
print(output)
