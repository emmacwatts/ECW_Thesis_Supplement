#!/usr/bin/env Rscript

# Import per-sample Kallisto estimates, aggregate transcripts to genes, and
# write two gene-by-group Excel tables. Each group is represented by adjacent
# GROUP_mean and GROUP_error columns, where error is the biological-replicate SEM.

usage <- function() {
    cat(
"Usage:
  Rscript summarise_expression.R \\
    --dataset hamel|grosse-holz \\
    --kallisto-dir DIR \\
    --output-dir DIR \\
    (--tx2gene FILE | --targets-are-genes) [--prefix NAME]

For another dataset:
  Rscript summarise_expression.R \\
    --metadata metadata.csv \\
    --kallisto-dir DIR \\
    --output-dir DIR \\
    (--tx2gene FILE | --targets-are-genes) [--prefix NAME]

Arguments:
  --dataset NAME           Use embedded metadata for hamel or grosse-holz.
  --metadata FILE          Custom CSV/TSV with run, group, and replicate columns.
  --kallisto-dir DIR       Directory containing RUN/abundance.tsv subdirectories.
  --output-dir DIR         Directory in which Excel workbooks are written.
  --tx2gene FILE           CSV/TSV with transcript_id and gene_id columns.
  --targets-are-genes      Treat every Kallisto target_id as a gene_id.
                           Use only when each reference record represents one gene.
  --prefix NAME            Output filename prefix (default: dataset name).
  --counts-from-abundance  tximport method: no, scaledTPM, lengthScaledTPM,
                           or dtuScaledTPM (default: lengthScaledTPM).
  -h, --help               Show this message.

Outputs:
  PREFIX_mean_TPM.xlsx
  PREFIX_mean_DESeq2_normalised_counts.xlsx

Each group has paired GROUP_mean and GROUP_error (SEM) columns.
", sep = "")
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)

parse_arguments <- function(args) {
    options <- list(
        dataset = NULL,
        metadata = NULL,
        kallisto_dir = NULL,
        output_dir = NULL,
        tx2gene = NULL,
        targets_are_genes = FALSE,
        prefix = NULL,
        counts_from_abundance = "lengthScaledTPM"
    )

    key_map <- c(
        "--dataset" = "dataset",
        "--metadata" = "metadata",
        "--kallisto-dir" = "kallisto_dir",
        "--output-dir" = "output_dir",
        "--tx2gene" = "tx2gene",
        "--prefix" = "prefix",
        "--counts-from-abundance" = "counts_from_abundance"
    )

    i <- 1L
    while (i <= length(args)) {
        argument <- args[[i]]
        if (argument %in% c("-h", "--help")) {
            usage()
            quit(save = "no", status = 0)
        } else if (argument == "--targets-are-genes") {
            options$targets_are_genes <- TRUE
            i <- i + 1L
        } else if (argument %in% names(key_map)) {
            if (i == length(args)) stopf("Missing value after %s", argument)
            options[[unname(key_map[[argument]])]] <- args[[i + 1L]]
            i <- i + 2L
        } else {
            stopf("Unknown argument: %s", argument)
        }
    }
    options
}

read_table_auto <- function(path) {
    if (!file.exists(path)) stopf("File not found: %s", path)
    first_line <- readLines(path, n = 1L, warn = FALSE)
    if (length(first_line) == 0L) stopf("File is empty: %s", path)
    separator <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
    read.table(
        path,
        header = TRUE,
        sep = separator,
        quote = "\"",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

hamel_metadata <- function() {
    condition <- rep(
        c("Non_infiltrated", "Mock", "EV", "P19", "P19_HA"),
        each = 3L
    )
    data.frame(
        run = c(
            "SRR24709124", "SRR24709123", "SRR24709122",
            "SRR24709121", "SRR24709120", "SRR24709119",
            "SRR24709118", "SRR24709117", "SRR24709116",
            "SRR24709115", "SRR24709114", "SRR24709113",
            "SRR24709112", "SRR24709111", "SRR24709110"
        ),
        condition = condition,
        group = condition,
        replicate = rep(1:3, times = 5L),
        stringsAsFactors = FALSE
    )
}

grosse_holz_metadata <- function() {
    groups <- list(
        Mock_2dpi = c("SRR5691069", "SRR5691068", "SRR5691071"),
        Mock_5dpi = c("SRR5691070", "SRR5691067", "SRR5691066"),
        Mock_7dpi = c("SRR5691044", "SRR5691045", "SRR5691042"),
        Mock_10dpi = c("SRR5691043", "SRR5691040", "SRR5691041"),
        WT_2dpi = c("SRR5691060", "SRR5691061", "SRR5691062"),
        WT_5dpi = c("SRR5691063", "SRR5691056", "SRR5691057"),
        WT_7dpi = c("SRR5691058", "SRR5691059", "SRR5691064"),
        WT_10dpi = c("SRR5691065", "SRR5691046", "SRR5691047"),
        P19_2dpi = c("SRR5691048", "SRR5691049", "SRR5691053"),
        P19_5dpi = c("SRR5691050", "SRR5691051", "SRR5691052"),
        P19_7dpi = c("SRR5691054", "SRR5691055", "SRR5691073"),
        P19_10dpi = c("SRR5691072", "SRR5691075", "SRR5691074")
    )

    do.call(rbind, lapply(names(groups), function(group_name) {
        components <- strsplit(group_name, "_", fixed = TRUE)[[1L]]
        data.frame(
            run = groups[[group_name]],
            condition = components[[1L]],
            time = sub("dpi$", "", components[[2L]]),
            group = group_name,
            replicate = c("A", "B", "C"),
            stringsAsFactors = FALSE
        )
    }))
}

summarise_by_group <- function(expression_matrix, groups, group_levels, output_names) {
    columns <- list()
    for (i in seq_along(group_levels)) {
        group_name <- group_levels[[i]]
        output_name <- output_names[[i]]
        values <- expression_matrix[, groups == group_name, drop = FALSE]
        columns[[paste0(output_name, "_mean")]] <- rowMeans(values)
        columns[[paste0(output_name, "_error")]] <- apply(values, 1L, stats::sd) / sqrt(ncol(values))
    }
    result <- do.call(cbind, columns)
    rownames(result) <- rownames(expression_matrix)
    result
}

to_output_table <- function(expression_matrix) {
    ordered_genes <- order(rownames(expression_matrix))
    data.frame(
        gene_id = rownames(expression_matrix)[ordered_genes],
        expression_matrix[ordered_genes, , drop = FALSE],
        check.names = FALSE,
        row.names = NULL
    )
}

options <- parse_arguments(commandArgs(trailingOnly = TRUE))

if (is.null(options$kallisto_dir)) stopf("--kallisto-dir is required")
if (is.null(options$output_dir)) stopf("--output-dir is required")
if (!dir.exists(options$kallisto_dir)) {
    stopf("Kallisto directory not found: %s", options$kallisto_dir)
}
if (!xor(is.null(options$dataset), is.null(options$metadata))) {
    stopf("Use exactly one of --dataset or --metadata")
}
if (!xor(!is.null(options$tx2gene), options$targets_are_genes)) {
    stopf("Use exactly one of --tx2gene or --targets-are-genes")
}

allowed_count_methods <- c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM")
if (!(options$counts_from_abundance %in% allowed_count_methods)) {
    stopf(
        "--counts-from-abundance must be one of: %s",
        paste(allowed_count_methods, collapse = ", ")
    )
}

required_packages <- c("tximport", "DESeq2", "writexl")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
    stopf("Missing R package(s): %s", paste(missing_packages, collapse = ", "))
}

if (!is.null(options$dataset)) {
    dataset_key <- tolower(gsub("_", "-", options$dataset, fixed = TRUE))
    if (dataset_key == "hamel") {
        metadata <- hamel_metadata()
        default_prefix <- "hamel"
        output_group_names <- c("Non-infiltrated", "Mock", "EV", "p19", "p19HA")
    } else if (dataset_key %in% c("grosse-holz", "grosseholz", "gh")) {
        metadata <- grosse_holz_metadata()
        default_prefix <- "grosse_holz"
        output_group_names <- c(
            "wt_mock_2", "wt_mock_5", "wt_mock_7", "wt_mock_10",
            "wt_ag_2", "wt_ag_5", "wt_ag_7", "wt_ag_10",
            "p19o_ag_2", "p19o_ag_5", "p19o_ag_7", "p19o_ag_10"
        )
    } else {
        stopf("Unknown embedded dataset: %s", options$dataset)
    }
} else {
    metadata <- read_table_auto(options$metadata)
    default_prefix <- tools::file_path_sans_ext(basename(options$metadata))
    output_group_names <- unique(metadata$group)
}

required_metadata_columns <- c("run", "group", "replicate")
missing_columns <- setdiff(required_metadata_columns, names(metadata))
if (length(missing_columns) > 0L) {
    stopf(
        "Metadata is missing required column(s): %s",
        paste(missing_columns, collapse = ", ")
    )
}
if (anyDuplicated(metadata$run)) stopf("Every run must occur once in metadata")
if (anyNA(metadata[, required_metadata_columns]) ||
    any(metadata$run == "") || any(metadata$group == "")) {
    stopf("run, group, and replicate metadata values must not be missing")
}

group_levels <- unique(metadata$group)
group_counts <- table(metadata$group)
if (any(group_counts < 2L)) {
    stopf(
        "Every group requires biological replication; insufficient group(s): %s",
        paste(names(group_counts)[group_counts < 2L], collapse = ", ")
    )
}

kallisto_files <- file.path(options$kallisto_dir, metadata$run, "abundance.tsv")
missing_kallisto <- kallisto_files[!file.exists(kallisto_files)]
if (length(missing_kallisto) > 0L) {
    stopf(
        "Missing Kallisto abundance file(s):\n%s",
        paste(missing_kallisto, collapse = "\n")
    )
}
names(kallisto_files) <- metadata$run

if (options$targets_are_genes) {
    first_abundance <- read.delim(
        kallisto_files[[1L]],
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (!("target_id" %in% names(first_abundance))) {
        stopf("Kallisto file does not contain a target_id column")
    }
    tx2gene <- unique(data.frame(
        transcript_id = first_abundance$target_id,
        gene_id = first_abundance$target_id,
        stringsAsFactors = FALSE
    ))
} else {
    tx2gene <- read_table_auto(options$tx2gene)
    required_mapping_columns <- c("transcript_id", "gene_id")
    missing_mapping <- setdiff(required_mapping_columns, names(tx2gene))
    if (length(missing_mapping) > 0L) {
        stopf(
            "Transcript-to-gene mapping is missing column(s): %s",
            paste(missing_mapping, collapse = ", ")
        )
    }
    tx2gene <- unique(tx2gene[, required_mapping_columns, drop = FALSE])
    if (anyDuplicated(tx2gene$transcript_id)) {
        stopf("Each transcript_id must map to exactly one gene_id")
    }
}

message("Importing ", length(kallisto_files), " Kallisto samples with tximport")
txi <- tximport::tximport(
    kallisto_files,
    type = "kallisto",
    tx2gene = tx2gene,
    countsFromAbundance = options$counts_from_abundance,
    dropInfReps = TRUE
)

metadata$group <- factor(metadata$group, levels = group_levels)
rownames(metadata) <- metadata$run

message("Fitting DESeq2 model: ~ group")
dds <- DESeq2::DESeqDataSetFromTximport(
    txi,
    colData = metadata,
    design = ~ group
)
dds <- DESeq2::DESeq(dds)

mean_tpm <- summarise_by_group(
    txi$abundance, metadata$group, group_levels, output_group_names
)
normalised_counts <- DESeq2::counts(dds, normalized = TRUE)
mean_normalised_counts <- summarise_by_group(
    normalised_counts, metadata$group, group_levels, output_group_names
)

dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- if (is.null(options$prefix)) default_prefix else options$prefix
tpm_path <- file.path(options$output_dir, paste0(prefix, "_mean_TPM.xlsx"))
deseq_path <- file.path(
    options$output_dir,
    paste0(prefix, "_mean_DESeq2_normalised_counts.xlsx")
)

writexl::write_xlsx(
    list(mean_TPM = to_output_table(mean_tpm)),
    path = tpm_path
)
writexl::write_xlsx(
    list(mean_DESeq2_normalised_counts = to_output_table(mean_normalised_counts)),
    path = deseq_path
)

message("Wrote descriptive TPM table: ", tpm_path)
message("Wrote within-gene comparison table: ", deseq_path)
message(
    "DESeq2-normalised counts are not TPM and should not be used as a formal ",
    "test of differences between genes."
)
