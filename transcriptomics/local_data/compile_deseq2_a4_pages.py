#!/usr/bin/env python3
"""Arrange DESeq2 target figures into alphabetical A4 portrait PNG pages."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from PIL import Image, ImageOps

DPI = 300
PAGE_WIDTH = round(210 / 25.4 * DPI)
PAGE_HEIGHT = round(297 / 25.4 * DPI)
MARGIN = 100
SLOT_GAP = 40
PLOTS_PER_PAGE = 3


def arguments() -> argparse.Namespace:
    transcriptomics_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir", type=Path,
        default=transcriptomics_dir / "output" / "deseq2",
    )
    parser.add_argument(
        "--gene-groups", type=Path,
        default=transcriptomics_dir / "gene_id_groups.example.csv",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=transcriptomics_dir / "local_data" / "deseq2_a4_pages",
    )
    return parser.parse_args()


def target_images(gene_groups: Path, source_dir: Path) -> list[tuple[str, Path]]:
    targets: dict[str, str] = {}
    with gene_groups.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"target", "image_file"}.issubset(reader.fieldnames):
            raise ValueError(f"{gene_groups} must contain target and image_file columns")
        for row in reader:
            target = (row.get("target") or "").strip()
            image_file = (row.get("image_file") or "").strip()
            if not target or not image_file:
                raise ValueError(f"{gene_groups} contains a blank target or image_file")
            if Path(image_file).name != image_file:
                raise ValueError(f"Expected a filename, got {image_file!r}")
            stem = Path(image_file).stem
            filename = f"{stem}.png" if stem.endswith("_DESeq2") else f"{stem}_DESeq2.png"
            if target in targets and targets[target] != filename:
                raise ValueError(f"Target {target!r} maps to more than one image")
            targets[target] = filename

    if not targets:
        raise ValueError(f"No targets found in {gene_groups}")
    ordered = sorted(targets.items(), key=lambda item: (item[0].casefold(), item[0]))
    result = [(target, source_dir / filename) for target, filename in ordered]
    missing = [str(path) for _, path in result if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing DESeq2 figure(s):\n" + "\n".join(missing))
    return result


def compile_pages(
    figures: list[tuple[str, Path]], output_dir: Path
) -> list[tuple[Path, list[str]]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    for old_page in output_dir.glob("page_*.png"):
        old_page.unlink()

    slot_width = PAGE_WIDTH - 2 * MARGIN
    slot_height = (PAGE_HEIGHT - 2 * MARGIN - (PLOTS_PER_PAGE - 1) * SLOT_GAP) // PLOTS_PER_PAGE
    pages = []
    for page_index, start in enumerate(range(0, len(figures), PLOTS_PER_PAGE), 1):
        page_items = figures[start:start + PLOTS_PER_PAGE]
        canvas = Image.new("RGB", (PAGE_WIDTH, PAGE_HEIGHT), "white")
        for slot_index, (_, image_path) in enumerate(page_items):
            with Image.open(image_path) as source:
                image = ImageOps.contain(
                    source.convert("RGBA"),
                    (slot_width, slot_height),
                    method=Image.Resampling.LANCZOS,
                )
            x = (PAGE_WIDTH - image.width) // 2
            y = MARGIN + slot_index * (slot_height + SLOT_GAP)
            y += (slot_height - image.height) // 2
            canvas.paste(image, (x, y), image)

        output_path = output_dir / f"page_{page_index:02d}.png"
        canvas.save(output_path, format="PNG", dpi=(DPI, DPI))
        pages.append((output_path, [target for target, _ in page_items]))
    return pages


def main() -> None:
    args = arguments()
    figures = target_images(args.gene_groups, args.source_dir)
    pages = compile_pages(figures, args.output_dir)
    print(f"Created {len(pages)} A4 portrait PNG pages in {args.output_dir}")
    for page_path, targets in pages:
        print(f"{page_path.name}: {', '.join(targets)}")


if __name__ == "__main__":
    main()