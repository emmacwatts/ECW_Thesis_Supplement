# RT-qPCR primer design

`primer_design.py` is a command-line implementation of the RT-qPCR primer-
design workflow used for the thesis. It designs primers against transcript-
oriented coding sequences (CDSs), retains primer pairs in which at least one
primer spans a CDS exon boundary, optionally rejects primers that overlap a
known VIGS fragment, and exports both the selected pair and the qualifying
alternatives.

The program is general for gene IDs represented in the supplied FASTA and GFF,
although the included example data and retained outputs use the LAB 3.60
*Nicotiana benthamiana* annotation.

## Files in this directory

```text
primer_design/
├── README.md
├── primer_design.py
├── data/
│   └── vigs_fragments.csv
└── output/
    ├── topPrimerSetPerNbL.xlsx
    ├── otherPrimers/
    └── geneLine/
```

`data/vigs_fragments.csv` records the physical VIGS constructs, their target
names, and intended NbLab360 target IDs. `output/` contains the retained thesis
run. The reference genome files are deliberately excluded because they are
large and are maintained by their original providers.

## Reference files

The thesis workflow used the chromosome-scale LAB 3.60 resource described by
[Ranawaka et al. (2023)](https://doi.org/10.1038/s41477-023-01489-8), specifically:

- `NbLab360CDSStripped.fasta`, containing transcript-oriented CDS sequences;
- `NbLab360.v103.gff3`, the matching version 1.0.3 annotation.

The NbLab360 assembly page provides the associated assembly and annotation
files: [Plant GARDEN: NbLab360](https://plantgarden.jp/en/list/t4100/genome/t4100.G004).
Download or derive the CDS FASTA and GFF before running the script, keep the
files in a local reference-data directory, and pass their paths with `--cds`
and `--gff`. Do not mix releases: the FASTA identifiers and concatenated CDS
lengths must agree with the GFF annotation.

If a provider supplies a CDS FASTA under a slightly different filename, that is
acceptable provided its record identifiers and sequences correspond to the
same v1.0.3 annotation. The script removes an isoform suffix such as `.1` when
matching IDs and preserves the first FASTA record encountered for each base
gene ID; the intended primary isoform should therefore occur first.

## Installation

Python 3.12 is recommended. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The Python dependencies are Biopython, pandas, openpyxl, Plotly, and
`primer3-py`. The virtual environment and local reference files should remain
untracked.

## Basic use

Run the program from the repository root with one or more base NbLab360 gene
IDs:

```bash
python primer_design/primer_design.py \
  NbL08g16530 NbL12g17000 \
  --cds /path/to/NbLab360CDSStripped.fasta \
  --gff /path/to/NbLab360.v103.gff3 \
  --output-dir /path/to/primer_results
```

The positional arguments are gene IDs. The optional arguments are:

| Argument | Meaning |
|---|---|
| `--cds PATH` | Transcript-oriented CDS FASTA. Default: `databases/NbLab360CDSStripped.fasta`. |
| `--gff PATH` | Matching GFF3 annotation. Default: `databases/NbLab360.v103.gff3`. |
| `--output-dir PATH` | Destination directory. Default: `rtqPCR`. |
| `--vigs GENE=SEQUENCE` | Exact VIGS sequence to exclude for one gene; repeat the option for multiple genes. |

Use `python primer_design/primer_design.py --help` for the current command-line
help.

## VIGS-aware design

Where a VIGS construct is an exact substring of the target CDS, provide it as
follows:

```bash
python primer_design/primer_design.py \
  NbL12g17000 \
  --cds /path/to/NbLab360CDSStripped.fasta \
  --gff /path/to/NbLab360.v103.gff3 \
  --vigs NbL12g17000=ACTG...GCTA \
  --output-dir /path/to/primer_results
```

Pairs are rejected if either primer overlaps the located VIGS interval. The
sequence must be an exact, forward-orientation substring of that CDS; otherwise
the gene is reported as failed. A physical VIGS construct may target multiple
homologues without being an exact substring of every CDS. Consequently,
`data/vigs_fragments.csv` must not be converted blindly to one `--vigs` value
per intended gene. Align each construct to each target, verify the relevant
interval and orientation, and document any mismatch-aware decision separately.

## Design and selection rules

Primer3 generates up to 10,000 candidate pairs using these principal settings:

| Property | Setting |
|---|---|
| Primer length | optimum 20 nt; range 18–30 nt |
| Melting temperature | optimum 60 °C; range 59–61 °C |
| GC content | 40–60% |
| Product length | 70–200 bp |
| Maximum homopolymer | 3 nt |
| Maximum hairpin melting temperature | 47 °C |

The candidates retain Primer3's penalty ordering. A pair qualifies only when:

1. at least one primer crosses a CDS exon boundary;
2. that primer contributes at least five nucleotides on each side of the
   boundary; and
3. neither primer overlaps the supplied VIGS interval, if one was supplied.

The first qualifying pair is selected. Exons are reconstructed in transcript
order from GFF `CDS` features, including reverse ordering for minus-strand
genes. All positions reported in the generated maps are zero-based, inclusive,
and measured from the 5′ end of the CDS rather than genomic coordinates.

## Outputs

The output directory contains:

- `topPrimerSetPerNbL.xlsx`: one selected pair or an explicit failure status
  for every requested gene;
- `otherPrimers/GENE.xlsx`: Primer3 diagnostic fields and every pair that met
  the boundary and VIGS criteria;
- `geneLine/GENE.html`: an interactive transcript-oriented map of CDS exons,
  the selected primers, and the VIGS interval when supplied.

The summary workbook reports the forward primer as
`PRIMER_LEFT_SEQUENCE` and the reverse primer in its ordered 5′→3′ form as
`PRIMER_RIGHT_SEQUENCE`. `CDSMatches` gives the number of exact occurrences of
the forward binding sequence and reverse-primer binding sequence across all CDS
records in the supplied FASTA.

## Validation and limitations

The generated pairs are candidates for experimental validation, not a guarantee
of assay performance or specificity. Before ordering primers:

- inspect the selected pair and the exon/primer map;
- confirm the intended transcript isoform and amplicon sequence;
- check genome-wide specificity with an appropriate alignment or primer-BLAST
  workflow, including plausible mismatches and paralogues;
- consider genomic DNA, splice variants, polymorphisms, and known VIGS-derived
  transcripts; and
- validate amplification efficiency, melt-curve specificity, and product size
  experimentally.

`CDSMatches` is only an exact-string screen against the supplied CDS FASTA. It
does not search genomic, intergenic, unannotated, or mismatch-containing sites,
and it is not a replacement for a full PCR-specificity analysis. Genes for
which no pair passes the hard filters remain in the summary with a failure
message; relaxing constraints for such outliers is a manual design decision and
is not performed automatically.

## Reference

Ranawaka, B., An, J., Lorenc, M. T., *et al.* (2023). A multi-omic
*Nicotiana benthamiana* resource for fundamental research and biotechnology.
*Nature Plants*, 9, 1558–1571.
[https://doi.org/10.1038/s41477-023-01489-8](https://doi.org/10.1038/s41477-023-01489-8)
