#!/usr/bin/env python3
"""Design RT-qPCR primers for NbLab360 N. benthamiana genes.

The program designs 70--200 bp products with Primer3 and prefers primer pairs in
which at least one primer spans a CDS exon boundary by five or more bases.  An
optional VIGS fragment can be supplied; primer pairs overlapping that fragment
are excluded.

Example
-------
python primerDesignUtils0201.py NbL08g16530 --output-dir rtqPCR_test
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Iterable


DEFAULT_CDS = Path("databases/NbLab360CDSStripped.fasta")
DEFAULT_GFF = Path("databases/NbLab360.v103.gff3")

PRIMER3_GLOBAL_ARGS = {
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 30,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 59.0,
    "PRIMER_MAX_TM": 61.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_POLY_X": 3,
    "PRIMER_MAX_HAIRPIN_TH": 47.0,
    "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 47.0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[70, 200]],
    "PRIMER_NUM_RETURN": 10_000,
}

PAIR_FIELDS = (
    "PRIMER_PAIR_PENALTY",
    "PRIMER_LEFT_PENALTY",
    "PRIMER_RIGHT_PENALTY",
    "PRIMER_LEFT_SEQUENCE",
    "PRIMER_RIGHT_SEQUENCE",
    "PRIMER_LEFT",
    "PRIMER_RIGHT",
    "PRIMER_LEFT_TM",
    "PRIMER_RIGHT_TM",
    "PRIMER_LEFT_GC_PERCENT",
    "PRIMER_RIGHT_GC_PERCENT",
    "PRIMER_LEFT_SELF_ANY_TH",
    "PRIMER_RIGHT_SELF_ANY_TH",
    "PRIMER_LEFT_SELF_END_TH",
    "PRIMER_RIGHT_SELF_END_TH",
    "PRIMER_LEFT_HAIRPIN_TH",
    "PRIMER_RIGHT_HAIRPIN_TH",
    "PRIMER_LEFT_END_STABILITY",
    "PRIMER_RIGHT_END_STABILITY",
    "PRIMER_PAIR_COMPL_ANY_TH",
    "PRIMER_PAIR_COMPL_END_TH",
    "PRIMER_PAIR_PRODUCT_SIZE",
)

GLOBAL_FIELDS = (
    "PRIMER_LEFT_EXPLAIN",
    "PRIMER_RIGHT_EXPLAIN",
    "PRIMER_PAIR_EXPLAIN",
    "PRIMER_LEFT_NUM_RETURNED",
    "PRIMER_RIGHT_NUM_RETURNED",
    "PRIMER_INTERNAL_NUM_RETURNED",
    "PRIMER_PAIR_NUM_RETURNED",
)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""

    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_cds_records(fasta_path: Path = DEFAULT_CDS) -> dict[str, str]:
    """Load transcript-oriented CDS records, keyed by ID without isoform suffix."""

    from Bio import SeqIO

    records: dict[str, str] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        gene_id = record.id.split(".")[0]
        # Preserve the first (normally .1) isoform, matching the original workflow.
        records.setdefault(gene_id, str(record.seq).upper())
    return records


def load_cds_exons(gff_path: Path = DEFAULT_GFF) -> dict[str, tuple[str, list[tuple[int, int]]]]:
    """Return ``gene -> (strand, transcript-relative inclusive exon intervals)``."""

    genomic: dict[str, list[tuple[int, int]]] = {}
    strands: dict[str, str] = {}
    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = dict(
                item.split("=", 1) for item in fields[8].split(";") if "=" in item
            )
            raw_id = attrs.get("Parent", attrs.get("ID", ""))
            gene_id = raw_id.split(".")[0]
            if not gene_id:
                continue
            start, stop = int(fields[3]), int(fields[4])
            genomic.setdefault(gene_id, []).append((start, stop))
            strands[gene_id] = fields[6]

    result: dict[str, tuple[str, list[tuple[int, int]]]] = {}
    for gene_id, intervals in genomic.items():
        strand = strands[gene_id]
        ordered = sorted(intervals, reverse=(strand == "-"))
        offset = 0
        relative: list[tuple[int, int]] = []
        for start, stop in ordered:
            length = stop - start + 1  # GFF coordinates are one-based and inclusive.
            relative.append((offset, offset + length - 1))
            offset += length
        result[gene_id] = strand, relative
    return result


def design_primers(template: str):
    """Run Primer3 and return global and per-pair pandas data frames."""

    import pandas as pd
    import primer3

    # design_primers is the current API; it produces the same result as the old
    # designPrimers compatibility alias used by the notebook.
    raw = primer3.bindings.design_primers(
        {"SEQUENCE_TEMPLATE": template}, PRIMER3_GLOBAL_ARGS
    )
    global_results = pd.DataFrame(
        [{field: raw.get(field) for field in GLOBAL_FIELDS}]
    )
    count = int(raw.get("PRIMER_PAIR_NUM_RETURNED", 0))
    def indexed_key(field: str, index: int) -> str:
        if field in {"PRIMER_LEFT", "PRIMER_RIGHT"}:
            return f"{field}_{index}"
        for prefix in ("PRIMER_PAIR_", "PRIMER_LEFT_", "PRIMER_RIGHT_"):
            if field.startswith(prefix):
                return field.replace(prefix, f"{prefix}{index}_", 1)
        raise ValueError(f"Unexpected Primer3 field: {field}")

    rows = [
        {field: raw.get(indexed_key(field, index)) for field in PAIR_FIELDS}
        for index in range(count)
    ]
    return global_results, pd.DataFrame(rows, columns=PAIR_FIELDS)


def _containing_exon(position: int, exons: list[tuple[int, int]]) -> int:
    for number, (start, stop) in enumerate(exons, start=1):
        if start <= position <= stop:
            return number
    raise ValueError(f"Template position {position} is outside the annotated CDS")


def annotate_pairs(pairs, exons: list[tuple[int, int]], vigs_region: tuple[int, int] | None):
    """Annotate and filter pairs using exon-boundary and VIGS criteria."""

    cases: list[str] = []
    overlaps: list[str] = []
    avoids_vigs: list[bool] = []

    for _, row in pairs.iterrows():
        pair_cases: list[str] = []
        pair_overlaps: list[str] = []
        pair_avoids_vigs = True
        left_start, left_length = row["PRIMER_LEFT"]
        right_3prime, right_length = row["PRIMER_RIGHT"]
        primer_intervals = (
            (left_start, left_start + left_length - 1),
            (right_3prime - right_length + 1, right_3prime),
        )
        for start, stop in primer_intervals:
            start_exon = _containing_exon(start, exons)
            stop_exon = _containing_exon(stop, exons)
            if start_exon != stop_exon:
                boundary = exons[start_exon - 1][1]
                minimum_overlap = min(boundary - start + 1, stop - boundary)
                pair_cases.append(
                    f"E-E Boundary at exon{start_exon}, exon{stop_exon}"
                )
                pair_overlaps.append(str(minimum_overlap))
            else:
                pair_cases.append(f"exon{start_exon}")
                pair_overlaps.append("no overlap")
            if vigs_region is not None:
                vigs_start, vigs_stop = vigs_region
                if max(start, vigs_start) <= min(stop, vigs_stop):
                    pair_avoids_vigs = False
        cases.append("; ".join(pair_cases))
        overlaps.append("; ".join(pair_overlaps))
        avoids_vigs.append(pair_avoids_vigs)

    annotated = pairs.copy()
    annotated["PrimerCases"] = cases
    annotated["overlapMin"] = overlaps
    if vigs_region is not None:
        annotated["avoidsVIGS"] = avoids_vigs

    acceptable = []
    for case, overlap, avoids in zip(cases, overlaps, avoids_vigs):
        span_ok = any(
            "E-E Boundary" in item and int(value) >= 5
            for item, value in zip(case.split("; "), overlap.split("; "))
            if value != "no overlap"
        )
        acceptable.append(span_ok and avoids)
    return annotated.loc[acceptable].reset_index(drop=True)


def exact_cds_match_counts(primers: Iterable[str], cds_records: dict[str, str]) -> list[int]:
    """Count exact occurrences of each primer-binding sequence in all CDS records."""

    requested = list(primers)
    wanted = set(requested)
    lengths = {len(primer) for primer in wanted}
    counts: Counter[str] = Counter()
    for sequence in cds_records.values():
        for length in lengths:
            for start in range(len(sequence) - length + 1):
                candidate = sequence[start : start + length]
                if candidate in wanted:
                    counts[candidate] += 1
    return [counts[primer] for primer in requested]


def write_gene_line(exons: list[tuple[int, int]], winner, path: Path, vigs_region=None) -> None:
    """Write an interactive HTML map of exons, selected primers, and VIGS region."""

    import plotly.graph_objects as go

    elements = {f"exon{i}": interval for i, interval in enumerate(exons, start=1)}
    left_start, left_length = winner["PRIMER_LEFT"]
    right_3prime, right_length = winner["PRIMER_RIGHT"]
    elements["leftP"] = (left_start, left_start + left_length - 1)
    elements["rightP"] = (right_3prime - right_length + 1, right_3prime)
    if vigs_region is not None:
        elements["vigsRegion"] = vigs_region

    figure = go.Figure()
    for index, (name, (start, stop)) in enumerate(elements.items()):
        figure.add_trace(
            go.Scatter(
                x=[start, stop], y=[0, 0], mode="lines", name=name,
                opacity=0.45, line={"width": 20 if index % 2 == 0 else 15},
            )
        )
    figure.update_xaxes(showgrid=False, title="CDS position (0-based)")
    figure.update_yaxes(showgrid=False, showticklabels=False, zeroline=True)
    figure.update_layout(height=220, plot_bgcolor="white")
    figure.write_html(path)


def design_gene(
    gene_id: str,
    output_dir: Path,
    cds_records: dict[str, str],
    annotations: dict[str, tuple[str, list[tuple[int, int]]]],
    vigs_sequence: str | None = None,
):
    """Design primers for one gene and write its candidate table and gene map."""

    import pandas as pd

    if gene_id not in cds_records:
        raise KeyError(f"{gene_id!r} was not found in the CDS FASTA")
    if gene_id not in annotations:
        raise KeyError(f"{gene_id!r} has no CDS annotation in the GFF")
    template = cds_records[gene_id]
    _, exons = annotations[gene_id]
    if sum(stop - start + 1 for start, stop in exons) != len(template):
        raise ValueError(f"CDS FASTA and GFF lengths disagree for {gene_id}")

    vigs_region = None
    if vigs_sequence:
        vigs_sequence = vigs_sequence.upper()
        start = template.find(vigs_sequence)
        if start < 0:
            raise ValueError(f"VIGS sequence was not found in the CDS for {gene_id}")
        vigs_region = (start, start + len(vigs_sequence) - 1)

    global_results, raw_pairs = design_primers(template)
    candidates = annotate_pairs(raw_pairs, exons, vigs_region)
    if candidates.empty:
        raise RuntimeError(f"No primer pair met the boundary/VIGS criteria for {gene_id}")

    winner = candidates.iloc[0].copy()
    binding_sequences = [
        winner["PRIMER_LEFT_SEQUENCE"],
        reverse_complement(winner["PRIMER_RIGHT_SEQUENCE"]),
    ]
    winner["CDSMatches"] = ", ".join(
        map(str, exact_cds_match_counts(binding_sequences, cds_records))
    )

    other_dir = output_dir / "otherPrimers"
    line_dir = output_dir / "geneLine"
    other_dir.mkdir(parents=True, exist_ok=True)
    line_dir.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(other_dir / f"{gene_id}.xlsx") as writer:
        global_results.to_excel(writer, sheet_name="GlobalResults", index=False)
        candidates.to_excel(writer, sheet_name="LocalResults", index=False)
    write_gene_line(exons, winner, line_dir / f"{gene_id}.html", vigs_region)
    return winner


def run(
    gene_ids: list[str],
    output_dir: Path,
    cds_path: Path = DEFAULT_CDS,
    gff_path: Path = DEFAULT_GFF,
    vigs: dict[str, str] | None = None,
):
    """Design primers for one or more genes and write the selected-pair workbook."""

    import pandas as pd

    cds_records = load_cds_records(cds_path)
    annotations = load_cds_exons(gff_path)
    vigs = vigs or {}
    winners = {}
    for gene_id in gene_ids:
        try:
            winner = design_gene(
                gene_id, output_dir, cds_records, annotations, vigs.get(gene_id)
            )
            winner["Status"] = "designed"
            winners[gene_id] = winner
        except (KeyError, ValueError, RuntimeError) as error:
            winners[gene_id] = {"Status": f"failed: {error}"}
            print(f"WARNING: {gene_id}: {error}")
    result = pd.DataFrame.from_dict(winners, orient="index")
    result.index.name = "gene"
    result.to_excel(output_dir / "topPrimerSetPerNbL.xlsx")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genes", nargs="+", help="NbLab360 gene ID(s), e.g. NbL08g16530")
    parser.add_argument("--output-dir", type=Path, default=Path("rtqPCR"))
    parser.add_argument("--cds", type=Path, default=DEFAULT_CDS)
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF)
    parser.add_argument(
        "--vigs",
        action="append",
        default=[],
        metavar="GENE=SEQUENCE",
        help="VIGS fragment to avoid; repeat for multiple genes",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    vigs = dict(item.split("=", 1) for item in args.vigs)
    result = run(args.genes, args.output_dir, args.cds, args.gff, vigs)
    print(result[["PRIMER_LEFT_SEQUENCE", "PRIMER_RIGHT_SEQUENCE", "CDSMatches"]])


if __name__ == "__main__":
    main()
