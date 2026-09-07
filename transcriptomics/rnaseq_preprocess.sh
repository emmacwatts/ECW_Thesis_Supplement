#!/usr/bin/env bash

# Download SRA runs, trim reads, run FastQC, and quantify against a
# transcript/CDS reference with Kallisto. Expression-table compilation and
# DESeq2 normalisation are performed separately by summarise_expression.R.

set -Eeuo pipefail

PROGRAM_NAME="$(basename "$0")"

usage() {
    cat <<EOF
Usage:
  ${PROGRAM_NAME} --accessions FILE --layout single|paired --output DIR \\
    (--reference TRANSCRIPTS.fa | --index TRANSCRIPTS.idx) [OPTIONS]

Required arguments:
  --accessions FILE       One SRA run accession (SRR/ERR/DRR) per line.
                         Blank lines and lines beginning with # are ignored.
  --layout TYPE           Library layout: single or paired.
  --output DIR            Output directory.
  --reference FILE        Transcript/CDS FASTA used to build a Kallisto index.
  --index FILE            Existing Kallisto index (alternative to --reference).

Single-end arguments:
  --fragment-length N     Mean library fragment length required by Kallisto.
  --fragment-sd N         Fragment-length standard deviation required by Kallisto.

Optional arguments:
  --threads N             Threads used by fasterq-dump, Trimmomatic, FastQC,
                         and Kallisto (default: 8).
  --bootstraps N          Kallisto bootstrap samples (default: 100).
  --leading N             Trimmomatic LEADING value (default: 3).
  --trailing N            Trimmomatic TRAILING value (default: 3).
  --window-size N         Trimmomatic sliding-window size (default: 4).
  --window-quality N      Trimmomatic sliding-window quality (default: 15).
  --min-length N          Trimmomatic minimum read length (default: 36).
  -h, --help              Show this help message.

Example:
  ${PROGRAM_NAME} \\
    --accessions metadata/hamel_sra.txt \\
    --layout single \\
    --reference references/NbLab360.v103.gff3.CDS.fasta \\
    --fragment-length 68 \\
    --fragment-sd 7 \\
    --output results/hamel \\
    --threads 8
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

compress_fastq() {
    local fastq="$1"
    if [[ -s "${fastq}.gz" ]]; then
        return
    fi
    [[ -s "$fastq" ]] || die "Expected FASTQ file was not produced: $fastq"
    gzip "$fastq"
}

ACCESSIONS_FILE=""
LAYOUT=""
OUTPUT_DIR=""
REFERENCE_FASTA=""
PROVIDED_INDEX=""
FRAGMENT_LENGTH=""
FRAGMENT_SD=""
THREADS=8
BOOTSTRAPS=100
LEADING=3
TRAILING=3
WINDOW_SIZE=4
WINDOW_QUALITY=15
MIN_LENGTH=36

while [[ $# -gt 0 ]]; do
    case "$1" in
        --accessions) ACCESSIONS_FILE="${2:-}"; shift 2 ;;
        --layout) LAYOUT="${2:-}"; shift 2 ;;
        --output) OUTPUT_DIR="${2:-}"; shift 2 ;;
        --reference) REFERENCE_FASTA="${2:-}"; shift 2 ;;
        --index) PROVIDED_INDEX="${2:-}"; shift 2 ;;
        --fragment-length) FRAGMENT_LENGTH="${2:-}"; shift 2 ;;
        --fragment-sd) FRAGMENT_SD="${2:-}"; shift 2 ;;
        --threads) THREADS="${2:-}"; shift 2 ;;
        --bootstraps) BOOTSTRAPS="${2:-}"; shift 2 ;;
        --leading) LEADING="${2:-}"; shift 2 ;;
        --trailing) TRAILING="${2:-}"; shift 2 ;;
        --window-size) WINDOW_SIZE="${2:-}"; shift 2 ;;
        --window-quality) WINDOW_QUALITY="${2:-}"; shift 2 ;;
        --min-length) MIN_LENGTH="${2:-}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1 (use --help for usage)" ;;
    esac
done

[[ -n "$ACCESSIONS_FILE" ]] || die "--accessions is required"
[[ -r "$ACCESSIONS_FILE" ]] || die "Cannot read accession file: $ACCESSIONS_FILE"
[[ "$LAYOUT" == "single" || "$LAYOUT" == "paired" ]] || \
    die "--layout must be 'single' or 'paired'"
[[ -n "$OUTPUT_DIR" ]] || die "--output is required"

if [[ -n "$REFERENCE_FASTA" && -n "$PROVIDED_INDEX" ]]; then
    die "Use either --reference or --index, not both"
elif [[ -n "$REFERENCE_FASTA" ]]; then
    [[ -r "$REFERENCE_FASTA" ]] || die "Cannot read reference FASTA: $REFERENCE_FASTA"
elif [[ -n "$PROVIDED_INDEX" ]]; then
    [[ -r "$PROVIDED_INDEX" ]] || die "Cannot read Kallisto index: $PROVIDED_INDEX"
else
    die "Either --reference or --index is required"
fi

if [[ "$LAYOUT" == "single" ]]; then
    [[ -n "$FRAGMENT_LENGTH" ]] || die "--fragment-length is required for single-end reads"
    [[ -n "$FRAGMENT_SD" ]] || die "--fragment-sd is required for single-end reads"
fi

for value in "$THREADS" "$BOOTSTRAPS" "$LEADING" "$TRAILING" \
             "$WINDOW_SIZE" "$WINDOW_QUALITY" "$MIN_LENGTH"; do
    [[ "$value" =~ ^[0-9]+$ ]] || die "Integer option received a non-integer value: $value"
done
[[ "$THREADS" -gt 0 ]] || die "--threads must be greater than zero"

if [[ "$LAYOUT" == "single" ]]; then
    [[ "$FRAGMENT_LENGTH" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "Invalid --fragment-length"
    [[ "$FRAGMENT_SD" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "Invalid --fragment-sd"
fi

require_command prefetch
require_command fasterq-dump
require_command trimmomatic
require_command fastqc
require_command kallisto
require_command gzip

mapfile -t RUNS < <(
    sed 's/\r$//' "$ACCESSIONS_FILE" |
        sed '/^[[:space:]]*#/d; /^[[:space:]]*$/d; s/^[[:space:]]*//; s/[[:space:]]*$//'
)
[[ ${#RUNS[@]} -gt 0 ]] || die "No accessions were found in $ACCESSIONS_FILE"

for run in "${RUNS[@]}"; do
    [[ "$run" =~ ^(SRR|ERR|DRR)[0-9]+$ ]] || die "Invalid SRA run accession: $run"
done

SRA_DIR="${OUTPUT_DIR}/sra"
RAW_DIR="${OUTPUT_DIR}/fastq/raw"
TRIM_DIR="${OUTPUT_DIR}/fastq/trimmed"
UNPAIRED_DIR="${OUTPUT_DIR}/fastq/unpaired"
FASTQC_RAW_DIR="${OUTPUT_DIR}/fastqc/raw"
FASTQC_TRIM_DIR="${OUTPUT_DIR}/fastqc/trimmed"
KALLISTO_DIR="${OUTPUT_DIR}/kallisto"
REFERENCE_DIR="${OUTPUT_DIR}/reference"
TEMP_DIR="${OUTPUT_DIR}/tmp"

mkdir -p "$SRA_DIR" "$RAW_DIR" "$TRIM_DIR" "$UNPAIRED_DIR" \
         "$FASTQC_RAW_DIR" "$FASTQC_TRIM_DIR" "$KALLISTO_DIR" \
         "$REFERENCE_DIR" "$TEMP_DIR"

if [[ -n "$PROVIDED_INDEX" ]]; then
    KALLISTO_INDEX="$PROVIDED_INDEX"
else
    KALLISTO_INDEX="${REFERENCE_DIR}/transcripts.idx"
    if [[ ! -s "$KALLISTO_INDEX" ]]; then
        printf '\nBuilding Kallisto index\n'
        kallisto index -i "$KALLISTO_INDEX" "$REFERENCE_FASTA"
    fi
fi

printf 'Processing %d SRA run(s) as %s-end libraries.\n' "${#RUNS[@]}" "$LAYOUT"

for run in "${RUNS[@]}"; do
    printf '\n[%s] Download and FASTQ conversion\n' "$run"

    if [[ "$LAYOUT" == "single" ]]; then
        expected_raw=("${RAW_DIR}/${run}.fastq.gz")
    else
        expected_raw=("${RAW_DIR}/${run}_1.fastq.gz" "${RAW_DIR}/${run}_2.fastq.gz")
    fi

    raw_complete=true
    for file in "${expected_raw[@]}"; do
        [[ -s "$file" ]] || raw_complete=false
    done

    if [[ "$raw_complete" == false ]]; then
        prefetch --max-size u --output-directory "$SRA_DIR" "$run"

        if [[ "$LAYOUT" == "single" ]]; then
            fasterq-dump \
                --threads "$THREADS" \
                --progress \
                --outdir "$RAW_DIR" \
                --temp "${TEMP_DIR}/${run}" \
                "${SRA_DIR}/${run}"
            compress_fastq "${RAW_DIR}/${run}.fastq"
        else
            fasterq-dump \
                --split-files \
                --threads "$THREADS" \
                --progress \
                --outdir "$RAW_DIR" \
                --temp "${TEMP_DIR}/${run}" \
                "${SRA_DIR}/${run}"
            compress_fastq "${RAW_DIR}/${run}_1.fastq"
            compress_fastq "${RAW_DIR}/${run}_2.fastq"
        fi
    else
        printf '[%s] Existing compressed FASTQ file(s) found; skipping download.\n' "$run"
    fi

    printf '[%s] Trimming\n' "$run"

    if [[ "$LAYOUT" == "single" ]]; then
        trimmed_read="${TRIM_DIR}/${run}.trimmed.fastq.gz"
        if [[ ! -s "$trimmed_read" ]]; then
            trimmomatic SE \
                -threads "$THREADS" \
                "${RAW_DIR}/${run}.fastq.gz" \
                "$trimmed_read" \
                "LEADING:${LEADING}" \
                "TRAILING:${TRAILING}" \
                "SLIDINGWINDOW:${WINDOW_SIZE}:${WINDOW_QUALITY}" \
                "MINLEN:${MIN_LENGTH}"
        else
            printf '[%s] Existing trimmed FASTQ found; skipping trimming.\n' "$run"
        fi
    else
        trimmed_forward="${TRIM_DIR}/${run}_1.paired.fastq.gz"
        trimmed_reverse="${TRIM_DIR}/${run}_2.paired.fastq.gz"
        unpaired_forward="${UNPAIRED_DIR}/${run}_1.unpaired.fastq.gz"
        unpaired_reverse="${UNPAIRED_DIR}/${run}_2.unpaired.fastq.gz"

        if [[ ! -s "$trimmed_forward" || ! -s "$trimmed_reverse" ]]; then
            trimmomatic PE \
                -threads "$THREADS" \
                "${RAW_DIR}/${run}_1.fastq.gz" \
                "${RAW_DIR}/${run}_2.fastq.gz" \
                "$trimmed_forward" \
                "$unpaired_forward" \
                "$trimmed_reverse" \
                "$unpaired_reverse" \
                "LEADING:${LEADING}" \
                "TRAILING:${TRAILING}" \
                "SLIDINGWINDOW:${WINDOW_SIZE}:${WINDOW_QUALITY}" \
                "MINLEN:${MIN_LENGTH}"
        else
            printf '[%s] Existing paired trimmed FASTQs found; skipping trimming.\n' "$run"
        fi
    fi

    printf '[%s] Kallisto quantification\n' "$run"
    sample_output="${KALLISTO_DIR}/${run}"

    if [[ -s "${sample_output}/abundance.tsv" ]]; then
        printf '[%s] Existing Kallisto output found; skipping quantification.\n' "$run"
    elif [[ "$LAYOUT" == "single" ]]; then
        kallisto quant \
            --threads "$THREADS" \
            --index "$KALLISTO_INDEX" \
            --output-dir "$sample_output" \
            --bootstrap-samples "$BOOTSTRAPS" \
            --single \
            --fragment-length "$FRAGMENT_LENGTH" \
            --sd "$FRAGMENT_SD" \
            "${TRIM_DIR}/${run}.trimmed.fastq.gz"
    else
        kallisto quant \
            --threads "$THREADS" \
            --index "$KALLISTO_INDEX" \
            --output-dir "$sample_output" \
            --bootstrap-samples "$BOOTSTRAPS" \
            "${TRIM_DIR}/${run}_1.paired.fastq.gz" \
            "${TRIM_DIR}/${run}_2.paired.fastq.gz"
    fi
done

printf '\nRunning FastQC on raw reads\n'
if [[ "$LAYOUT" == "single" ]]; then
    fastqc --threads "$THREADS" --outdir "$FASTQC_RAW_DIR" \
        "${RAW_DIR}"/*.fastq.gz
    printf '\nRunning FastQC on trimmed reads\n'
    fastqc --threads "$THREADS" --outdir "$FASTQC_TRIM_DIR" \
        "${TRIM_DIR}"/*.trimmed.fastq.gz
else
    fastqc --threads "$THREADS" --outdir "$FASTQC_RAW_DIR" \
        "${RAW_DIR}"/*_[12].fastq.gz
    printf '\nRunning FastQC on retained paired reads\n'
    fastqc --threads "$THREADS" --outdir "$FASTQC_TRIM_DIR" \
        "${TRIM_DIR}"/*.paired.fastq.gz
fi

printf '\nPipeline complete. Results: %s\n' "$OUTPUT_DIR"
printf 'Next step: run summarise_expression.R on %s/kallisto\n' "$OUTPUT_DIR"
