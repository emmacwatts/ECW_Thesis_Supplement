# ECW Thesis Supplement

Code and supporting data accompanying Emma C. Watts' thesis.

## RT-qPCR primer design

`primer_design/primer_design.py` designs RT-qPCR primers for Lab3.6
*Nicotiana benthamiana* gene models. It uses Primer3 to generate 70--200 bp
products and prioritises pairs in which at least one primer overlaps a CDS exon
boundary by at least five nucleotides. When a VIGS fragment can be located in a
target CDS, primer pairs overlapping that region can be excluded.

The large NbLab360 reference FASTA and GFF files are not included
in this repository. Users must obtain them from the LAB 3.6 genome annotation [Ranawaka et al., 2023](https://www.nature.com/articles/s41477-023-01489-8) and provide their paths on
the command line.

### Installation

Python 3.12 is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### Required reference files

The primer-design workflow was tested with:

- `NbLab360CDSStripped.fasta`: transcript-oriented NbLab360 CDS sequences.
- `NbLab360.v103.gff3`: the matching NbLab360 v1.0.3 genome annotation.

The FASTA and GFF must describe the same annotation release. The program checks
that the concatenated CDS lengths in the GFF agree with each selected FASTA
record. Reference files may be stored anywhere and should not be copied into
this repository.

### Example

```bash
python primer_design/primer_design.py \
  NbL08g16530 \
  --cds /path/to/NbLab360CDSStripped.fasta \
  --gff /path/to/NbLab360.v103.gff3 \
  --output-dir results
```

The output directory contains:

- `topPrimerSetPerNbL.xlsx`: selected pairs and status for every requested gene.
- `otherPrimers/`: qualifying Primer3 candidate tables.
- `geneLine/`: interactive transcript-oriented exon and primer maps.

The `output/` directory contains the corrected thesis run, restricted to genes
listed in `vigs_fragments.csv`.

### VIGS metadata

`primer_design/data/vigs_fragments.csv` records the VIGS fragments,
their biological target names, and all intended NbLab360 targets. Several
constructs target two homologues.

Before rerunning VIGS-aware primer exclusion, each construct should be aligned
to every intended target and the matched target interval recorded. The current
`--vigs GENE=SEQUENCE` option is appropriate only where the supplied sequence is
an exact substring of that gene's CDS.

## Reproducibility notes

- Coordinates in generated plots are zero-based, inclusive, and oriented from
  the 5' end of the CDS, including for minus-strand genes.
- `CDSMatches` counts exact occurrences of the two primer-binding sequences in
  the supplied CDS FASTA. Exact counting is a useful screen but is not a full
  PCR-specificity analysis allowing mismatches.
- Genes for which no pair met the exon-boundary and VIGS criteria retain an
  explicit failure status in the summary workbook for manual redesign if desired.
