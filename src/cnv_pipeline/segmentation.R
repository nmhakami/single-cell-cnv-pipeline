#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(DNAcopy)
})

option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character",
        help = "Input GC-corrected CNV profile TSV."
    ),
    make_option(
        c("-o", "--output"),
        type = "character",
        help = "Output segmented CNV profile TSV."
    ),
    make_option(
        c("--column"),
        type = "character",
        default = "gc_corrected",
        help = "Column containing the signal to segment."
    ),
    make_option(
        c("--alpha"),
        type = "double",
        default = 0.01,
        help = "CBS significance threshold."
    )
)

option_parser <- OptionParser(
    option_list = option_list,
    description = "Segment a single-cell CNV profile using circular binary segmentation."
)

args <- parse_args(option_parser)

if (is.null(args$input) || is.null(args$output)) {
    print_help(option_parser)
    stop("Input and output files are required.")
}

profile <- read.delim(
    args$input,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

required_columns <- c(
    "chrom",
    "start",
    "end",
    args$column
)

missing_columns <- setdiff(required_columns, colnames(profile))

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Missing required columns:",
            paste(missing_columns, collapse = ", ")
        )
    )
}

profile <- profile[
    is.finite(profile[[args$column]]),
    ,
    drop = FALSE
]

if (nrow(profile) < 2) {
    stop("Not enough valid bins for segmentation.")
}

profile$chrom <- as.character(profile$chrom)

cna_data <- CNA(
    genomdat = profile[[args$column]],
    chrom = profile$chrom,
    maploc = profile$start,
    data.type = "logratio"
)

segmented <- segment(
    cna_data,
    alpha = args$alpha,
    min.width = 2,
    verbose = 0
)

segments <- segmented$output

segments <- segments[
    ,
    c(
        "chrom",
        "loc.start",
        "loc.end",
        "num.mark",
        "seg.mean"
    ),
    drop = FALSE
]

colnames(segments) <- c(
    "chrom",
    "start",
    "end",
    "num_bins",
    "segment_mean"
)

write.table(
    segments,
    file = args$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)