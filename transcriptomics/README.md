# Transcriptomics workflow and plot compiler

This directory contains a reusable three-stage RNA-seq workflow:

1. `rnaseq_preprocess.sh` downloads SRA runs, converts and trims reads, runs
   FastQC, and quantifies a transcript/CDS reference with Kallisto.
2. `summarise_expression.R` imports Kallisto estimates, aggregates transcripts
   to genes, normalises counts with DESeq2, and writes mean ± SEM tables.
3. `compile_transcriptomics_plots.py` extracts any requested gene IDs from two
   similarly formatted summary tables and creates per-target plots, a combined
   PDF, and a long-form CSV of the plotted values. Optional DESeq2 summary
   inputs create a second, within-gene comparison figure set.

The scripts are reusable; `gene_id_groups.example.csv` is the thesis-specific
example used to select 73 NbLab360 IDs in 38 named target groups. Replace that
file with another mapping to plot a different gene set.

## Directory contents

```text
transcriptomics/
├── README.md
├── compile_transcriptomics_plots.py
├── gene_id_groups.example.csv
├── requirements.txt
├── rnaseq_preprocess.sh
└── summarise_expression.R
```

Large reference sequences, SRA/FASTQ files, Kallisto output, and full generated
Excel workbooks are intentionally excluded from Git. Keep them in a working
directory outside the repository.

For the thesis run, full workbooks are retained under the ignored
`transcriptomics/local_data/expression_tables/` directory. Compact selected-
gene plots and plotted-value CSVs are retained under the version-controlled
`transcriptomics/output/` directory.

## Installation

Python 3.12 is recommended for the plotting stage. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r transcriptomics/requirements.txt
```

The preprocessing stage also requires Bash, gzip, SRA Toolkit, Trimmomatic,
FastQC, and Kallisto. Versions used in the reconstructed analysis were SRA
Toolkit 3.0.10, Trimmomatic 0.39, FastQC 0.12.0, and Kallisto 0.50.1.

The summary stage requires R plus `tximport`, `DESeq2`, and `writexl`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("tximport", "DESeq2"))
install.packages("writexl")
```

## 1. Preprocess and quantify reads

Create a text file containing one SRA run accession per line. Blank lines and
lines beginning with `#` are ignored. Kallisto needs a transcript or CDS FASTA,
not an unprocessed genomic FASTA.

Paired-end example:

```bash
transcriptomics/rnaseq_preprocess.sh \
  --accessions metadata/accessions.txt \
  --layout paired \
  --reference /path/to/transcripts.fasta \
  --output /path/to/results/dataset \
  --threads 8
```

For single-end libraries, also supply the library fragment-length mean and SD:

```bash
transcriptomics/rnaseq_preprocess.sh \
  --accessions metadata/accessions.txt \
  --layout single \
  --reference /path/to/transcripts.fasta \
  --fragment-length 68 \
  --fragment-sd 7 \
  --output /path/to/results/dataset
```

Use `--index` instead of `--reference` to supply an existing Kallisto index.
Run the script with `--help` for trimming and bootstrap options. The workflow is
restartable: complete FASTQ, trimmed-read, index, and abundance outputs are
detected and skipped.

## 2. Build expression-summary tables

For a custom dataset, supply CSV or TSV metadata with exactly these required
columns:

| run | group | replicate |
|---|---|---|
| SRR000001 | control | 1 |
| SRR000002 | control | 2 |
| SRR000003 | treatment | 1 |

Each group must contain at least two biological replicates. Supply either a
transcript-to-gene table with `transcript_id,gene_id` columns or
`--targets-are-genes` when every Kallisto target is already a gene.

```bash
Rscript transcriptomics/summarise_expression.R \
  --metadata metadata/samples.csv \
  --kallisto-dir /path/to/results/dataset/kallisto \
  --tx2gene /path/to/tx2gene.csv \
  --output-dir /path/to/results/dataset/expression_tables \
  --prefix dataset
```

Embedded run metadata are also available through `--dataset hamel` and
`--dataset grosse-holz`. These reproduce the two studies used for the thesis.

The R script writes:

- `PREFIX_mean_TPM.xlsx`, for descriptive abundance plots;
- `PREFIX_mean_DESeq2_normalised_counts.xlsx`, for within-gene comparisons.

Both contain `gene_id` followed by adjacent `GROUP_mean` and `GROUP_error`
columns. Error is the standard error of the biological-replicate mean. The two
studies are normalised separately and these files are not formal differential-
expression results; no fold changes or adjusted P-values are calculated.

The retained thesis workbooks are named:

- `grosse_holz_mean_TPM.xlsx`;
- `hamel_mean_TPM.xlsx`;
- `grosse_holz_mean_DESeq2_normalised_counts.xlsx`;
- `hamel_mean_DESeq2_normalised_counts.xlsx`.

## 3. Compile selected-gene plots

The plot compiler accepts CSV, TSV, or Excel summaries. Each summary needs:

- an identifier column named `gene_id` or `target_id`;
- one or more paired `CONDITION_mean` and `CONDITION_error` columns.

The selection mapping must contain:

| target | gene_id | image_file |
|---|---|---|
| Example target | Gene001 | 01_example.png |
| Example target | Gene002 | 01_example.png |

`target` groups one or more gene IDs into a figure. `image_file` is optional;
if omitted, filenames are generated from target names. Thus the example NbLab
ID list is data, not program logic.

The main table must use time-course condition names ending in a numeric time,
such as `wt_mock_2_mean` and `wt_ag_2_mean`. Different series prefixes can be
selected with `--mock-series` and `--treatment-series`. The second (Hamel)
summary may contain any condition names.

Thesis example:

```bash
python transcriptomics/compile_transcriptomics_plots.py \
  --main-summary /path/to/tpmFull_withError.xlsx \
  --hamel-summary /path/to/tpmFull_withError_Hamel.xlsx \
  --gene-groups transcriptomics/gene_id_groups.example.csv \
  --output-dir /path/to/results/selected_gene_plots
```

Outputs are one PNG per target, `gene_expression_plots.pdf`, and
`gene_expression_values.csv`. When a target contains multiple gene IDs, each
gene is retained separately as grouped bars and as labelled time-course lines.
The long-form CSV likewise retains every gene separately.

### DESeq2-normalised comparison plots

DESeq2 enters between Kallisto quantification and plotting:

1. `tximport` imports Kallisto estimates and aggregates them to genes.
2. `DESeqDataSetFromTximport(..., design = ~ group)` constructs each study's
   DESeq2 model.
3. `DESeq()` estimates size factors and dispersions and fits that model.
4. `counts(..., normalized = TRUE)` supplies normalised replicate values.
5. `summarise_expression.R` writes group mean and SEM columns to
   `PREFIX_mean_DESeq2_normalised_counts.xlsx`.
6. The compiler uses those workbooks when both DESeq2 options below are given.

Inferential bootstrap replicates from Kallisto HDF5 files are not imported for
this summary step (`dropInfReps = TRUE`); DESeq2 uses estimated counts from the
`abundance.tsv` files. This avoids unnecessary bootstrap arrays and supports
retaining only TSV estimates for this stage.

```bash
python transcriptomics/compile_transcriptomics_plots.py \
  --main-summary /path/to/grosse_holz_mean_TPM.xlsx \
  --hamel-summary /path/to/hamel_mean_TPM.xlsx \
  --main-deseq-summary /path/to/grosse_holz_mean_DESeq2_normalised_counts.xlsx \
  --hamel-deseq-summary /path/to/hamel_mean_DESeq2_normalised_counts.xlsx \
  --gene-groups transcriptomics/gene_id_groups.example.csv \
  --output-dir /path/to/results/selected_gene_plots
```

The existing TPM plots are still produced. A `deseq2/` subdirectory additionally
contains one PNG per named target, `gene_expression_plots_DESeq2.pdf`, and
`gene_expression_values_DESeq2.csv`.

Each DESeq2 figure has one row per NbL ID and two study panels per row. Thus a
target represented by two NbL IDs has four panels: Hamel and Grosse-Holz for
the first ID, followed by Hamel and Grosse-Holz for the second. This preserves
within-gene comparisons and avoids averaging normalised counts across homologues.

These are descriptive plots of DESeq2-normalised counts. They do not replace
formal contrasts and do not report log2 fold changes, P-values, or adjusted
P-values.

## Included thesis outputs

`output/tpm/` contains the 38 corrected TPM figures, their combined PDF, and
the long-form plotted values. `output/deseq2/` contains 38 corresponding
DESeq2-normalised figures, their combined PDF, and their long-form values.

The regenerated Hamel TPM summary reproduces the historical workbook to
numerical precision. Regenerating the Grosse-Holz summary corrected historical
P19 2-day/5-day assignments to follow the explicit SRA-to-condition metadata
embedded in `summarise_expression.R`. Both plot collections therefore use the
regenerated, consistently assigned summaries.

## Dataset-specific notes

The embedded Hamel dataset contains 15 single-end libraries: five conditions
with three biological replicates each (SRA study SRP439157). Historical
reprocessing used fragment length 68 and SD 7; these were analysis assumptions,
not generally applicable defaults.

The embedded Grosse-Holz dataset contains 36 paired-end libraries: three
treatments, four time points, and three biological replicates (SRA study
SRP109347). `WT` in its deposited aliases denotes the agroinfiltrated condition.

The original analysis quantified against the LAB 3.60 *Nicotiana benthamiana*
CDS reference (`NbLab360.v103.gff3.CDS.fasta`). Use a matching annotation and
transcript-to-gene mapping, and record their source/version. Reference data are
not redistributed here.

## Interpretation and reproducibility

- Inspect FastQC reports before interpreting expression estimates.
- TPM is suitable for descriptive comparisons within a consistently processed
  dataset, not formal differential-expression testing.
- DESeq2-normalised counts are suitable for comparing conditions for the same
  gene; they are not TPM and should not be used for absolute comparisons among
  genes.
- Hamel and Grosse-Holz outputs are not cross-study batch corrected.
- Biological sample groups must come from verified metadata, never accession
  order alone.
