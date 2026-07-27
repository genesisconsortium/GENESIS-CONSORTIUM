#!/usr/bin/env Rscript

# Refit dreamlet and regenerate Pearson residuals from the saved Step 3
# final_processAssays objects, retaining only samples that pass final sex QC.
#
# Matching rule:
#   colnames(final_processAssays) <-> sex_qc_samples.csv$individualID
#
# Usage:
#   Rscript 03_rerun_pearson_residuals_sex_qc_pass.R \
#     <sex_qc_samples.csv> <processassays_dir> <output_dir> <file_region>
#
# Inputs:
#   sex_qc_samples.csv
#     Must contain individualID and final_sex_qc.
#
#   processassays_dir
#     Directory containing both class and subclass files ending in
#     _final_processAssays.RDS.
#
# Output:
#   output_dir/residuals/
#     Pearson residual matrices with primary genotype/output_id column names.
#
#   output_dir/qc/
#     Per-object sample matching reports and a combined run summary.

.libPaths(c("/sc/arion/projects/psychAD/aging/kiran/RLib_4_5", .libPaths()))

suppressPackageStartupMessages({
  library(dreamlet)
  library(variancePartition)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(R.utils)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript 03_rerun_pearson_residuals_sex_qc_pass.R ",
    "<sex_qc_samples.csv> <processassays_dir> <output_dir> <file_region>"
  )
}

SEX_QC_FILE <- args[1]
PROCESSASSAYS_DIR <- args[2]
OUTDIR <- args[3]
FILE_REGION <- args[4]

RESIDUAL_DIR <- file.path(OUTDIR, "residuals")
QC_DIR <- file.path(OUTDIR, "qc")

dir.create(RESIDUAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

form <- ~ sex + scale(age) +
  log(n_genes) + percent_mito +
  mito_genes + mito_ribo + ribo_genes

required_sex_qc_columns <- c(
  "file_region",
  "individualID",
  "final_sex_qc"
)

sex_qc <- read.csv(
  SEX_QC_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

missing_sex_qc_columns <- setdiff(required_sex_qc_columns, colnames(sex_qc))
if (length(missing_sex_qc_columns) > 0) {
  stop(
    "Missing required sex-QC columns: ",
    paste(missing_sex_qc_columns, collapse = ", ")
  )
}

sex_qc$individualID <- trimws(sex_qc$individualID)
sex_qc$final_sex_qc <- trimws(sex_qc$final_sex_qc)
sex_qc$file_region <- trimws(sex_qc$file_region)
sex_qc <- sex_qc[
  sex_qc$file_region == FILE_REGION,
  ,
  drop = FALSE
]

if (nrow(sex_qc) == 0) {
  stop("No sex-QC rows found for file_region: ", FILE_REGION)
}

sex_qc <- sex_qc[
  !is.na(sex_qc$individualID) & sex_qc$individualID != "",
  ,
  drop = FALSE
]

# An individual is eligible only when it has at least one QC record and every
# record for that individual is "pass". This prevents a duplicated individual
# with conflicting pass/fail records from being retained accidentally.
qc_status_by_individual <- aggregate(
  x = list(final_sex_qc = sex_qc$final_sex_qc),
  by = list(individualID = sex_qc$individualID),
  FUN = function(x) {
    x <- unique(x[!is.na(x) & x != ""])
    if (length(x) == 1 && identical(x, "pass")) {
      "pass"
    } else if (length(x) == 0) {
      "missing"
    } else {
      paste(sort(x), collapse = ";")
    }
  }
)

pass_individual_ids <- qc_status_by_individual$individualID[
  qc_status_by_individual$final_sex_qc == "pass"
]

cat("Sex-QC file:", SEX_QC_FILE, "\n")
cat("File region:", FILE_REGION, "\n")
cat("Individuals represented:", nrow(qc_status_by_individual), "\n")
cat("Individuals passing final sex QC:", length(pass_individual_ids), "\n")

processassays_files <- list.files(
  PROCESSASSAYS_DIR,
  pattern = "_(class|subclass)_final_processAssays\\.RDS$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(processassays_files) == 0) {
  stop(
    "No class/subclass *_final_processAssays.RDS files found in: ",
    PROCESSASSAYS_DIR
  )
}

get_level <- function(path) {
  filename <- basename(path)
  if (grepl("_subclass_final_processAssays\\.RDS$", filename)) {
    return("subclass")
  }
  if (grepl("_class_final_processAssays\\.RDS$", filename)) {
    return("class")
  }
  NA_character_
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

rename_cols_to_output_id <- function(mat, id_map) {
  mat <- as.matrix(mat)

  if (is.null(colnames(mat)) && !is.null(rownames(mat))) {
    row_match <- match(rownames(mat), id_map$processassays_colname)
    if (any(!is.na(row_match))) {
      mat <- t(mat)
    }
  }

  if (is.null(colnames(mat))) {
    stop("Pearson residual matrix has no column names.")
  }

  idx <- match(colnames(mat), id_map$processassays_colname)

  if (any(is.na(idx)) && !is.null(rownames(mat))) {
    row_match <- match(rownames(mat), id_map$processassays_colname)
    if (sum(!is.na(row_match)) > sum(!is.na(idx))) {
      mat <- t(mat)
      idx <- match(colnames(mat), id_map$processassays_colname)
    }
  }

  if (any(is.na(idx))) {
    stop(
      "Could not map residual columns to output_id: ",
      paste(colnames(mat)[is.na(idx)], collapse = ", ")
    )
  }

  output_ids <- as.character(id_map$output_id[idx])

  if (any(is.na(output_ids)) || any(output_ids == "")) {
    stop("Missing output_id values after residual-column matching.")
  }

  if (anyDuplicated(output_ids)) {
    duplicated_ids <- unique(output_ids[duplicated(output_ids)])
    stop(
      "Duplicate output_id values after renaming: ",
      paste(duplicated_ids, collapse = ", ")
    )
  }

  colnames(mat) <- output_ids
  mat
}

run_summaries <- list()

for (processassays_file in processassays_files) {
  level <- get_level(processassays_file)
  object_stem <- sub(
    "_final_processAssays\\.RDS$",
    "",
    basename(processassays_file)
  )

  cat("\n============================================================\n")
  cat("Processing:", basename(processassays_file), "\n")
  cat("Level:", level, "\n")

  res.proc <- readRDS(processassays_file)

  object_coldata <- as.data.frame(colData(res.proc))
  processassays_ids <- rownames(object_coldata)

  if (is.null(processassays_ids) || any(processassays_ids == "")) {
    stop(
      "processAssays colData has missing row names: ",
      processassays_file
    )
  }

  if (!"output_id" %in% colnames(object_coldata)) {
    stop(
      "processAssays colData does not contain output_id: ",
      processassays_file
    )
  }

  full_id_map <- data.frame(
    processassays_colname = processassays_ids,
    output_id = as.character(object_coldata$output_id),
    stringsAsFactors = FALSE
  )

  assay_names <- assayNames(res.proc)
  if (length(assay_names) == 0) {
    stop("No assays remain in: ", basename(processassays_file))
  }

  res.proc.filtered <- res.proc
  assay_audits <- list()

  # A dreamletProcessedData object is a named list of assay-specific ELists.
  # Each assay can retain a different set of individuals, so filter columns
  # separately within every list element.
  for (cell_type in assay_names) {
    assay_ids <- colnames(res.proc[[cell_type]])

    if (is.null(assay_ids)) {
      stop(
        "Processed assay has no column names: ",
        cell_type,
        " in ",
        basename(processassays_file)
      )
    }

    assay_qc_idx <- match(assay_ids, qc_status_by_individual$individualID)
    assay_matched <- !is.na(assay_qc_idx)

    if (!all(assay_matched)) {
      stop(
        "Not all samples in assay ",
        cell_type,
        " matched file_region ",
        FILE_REGION,
        " in sex_qc$individualID. Unmatched IDs: ",
        paste(assay_ids[!assay_matched], collapse = ", ")
      )
    }

    assay_keep <- assay_ids %in% pass_individual_ids

    assay_audit <- data.frame(
      cell_type = cell_type,
      processassays_colname = assay_ids,
      output_id = full_id_map$output_id[
        match(assay_ids, full_id_map$processassays_colname)
      ],
      individualID_in_sex_qc = assay_matched,
      final_sex_qc = qc_status_by_individual$final_sex_qc[
        assay_qc_idx
      ],
      keep_final_sex_qc_pass = assay_keep,
      stringsAsFactors = FALSE
    )
    assay_audits[[length(assay_audits) + 1L]] <- assay_audit

    cat(
      "Assay:",
      cell_type,
      "| before:",
      length(assay_ids),
      "| retained:",
      sum(assay_keep),
      "| removed:",
      sum(!assay_keep),
      "\n"
    )

    if (!all(assay_ids %in% full_id_map$processassays_colname)) {
      stop(
        "Some samples in assay ",
        cell_type,
        " are absent from rownames(colData(res.proc))."
      )
    }

    if (sum(assay_keep) == 0) {
      stop(
        "No samples pass final sex QC for assay ",
        cell_type
      )
    }

    res.proc.filtered[[cell_type]] <-
      res.proc.filtered[[cell_type]][, assay_keep]
  }

  audit_file <- file.path(
    QC_DIR,
    paste0(object_stem, "_final_sex_qc_sample_audit.csv")
  )
  write.csv(
    do.call(rbind, assay_audits),
    audit_file,
    row.names = FALSE,
    quote = FALSE
  )

  retained_ids <- unique(unlist(lapply(res.proc.filtered, colnames)))
  data_filtered <- object_coldata[
    rownames(object_coldata) %in% retained_ids,
    ,
    drop = FALSE
  ]

  cat("Running dreamlet on the filtered processAssays object...\n")
  fit <- dreamlet(
    res.proc.filtered,
    form,
    data = data_filtered
  )

  n_written <- 0L

  for (cell_type in names(fit)) {
    cat("Pearson residuals:", cell_type, "\n")

    pearson_residuals <- residuals(
      fit[[cell_type]],
      res.proc.filtered[[cell_type]],
      type = "pearson"
    )

    pearson_residuals <- rename_cols_to_output_id(
      pearson_residuals,
      full_id_map
    )

    residual_file <- file.path(
      RESIDUAL_DIR,
      paste0(
        object_stem,
        "_",
        safe_name(cell_type),
        "_residualsPearson.tsv"
      )
    )

    write.table(
      pearson_residuals,
      file = residual_file,
      quote = FALSE,
      sep = "\t",
      row.names = TRUE,
      col.names = NA,
      na = "NA"
    )

    R.utils::gzip(residual_file, overwrite = TRUE)
    n_written <- n_written + 1L

    rm(pearson_residuals)
    gc()
  }

  all_assay_ids <- unique(unlist(lapply(res.proc, colnames)))
  all_retained_ids <- unique(unlist(lapply(res.proc.filtered, colnames)))
  n_before <- length(all_assay_ids)
  n_matched <- sum(all_assay_ids %in% qc_status_by_individual$individualID)
  n_keep <- length(all_retained_ids)
  n_remove <- length(setdiff(all_assay_ids, all_retained_ids))

  run_summaries[[length(run_summaries) + 1L]] <- data.frame(
    processassays_file = basename(processassays_file),
    level = level,
    n_samples_before = n_before,
    n_samples_matched = n_matched,
    n_samples_retained = n_keep,
    n_samples_removed = n_remove,
    n_assays_in_object = length(assay_names),
    n_pearson_files_written = n_written,
    audit_file = audit_file,
    stringsAsFactors = FALSE
  )

  rm(
    res.proc,
    res.proc.filtered,
    object_coldata,
    full_id_map,
    assay_audits,
    data_filtered,
    fit
  )
  gc()
}

run_summary <- do.call(rbind, run_summaries)
summary_file <- file.path(QC_DIR, "final_sex_qc_pearson_residual_run_summary.csv")
write.csv(run_summary, summary_file, row.names = FALSE, quote = FALSE)

cat("\nCompleted class and subclass Pearson residual regeneration.\n")
cat("Residual directory:", RESIDUAL_DIR, "\n")
cat("QC summary:", summary_file, "\n")

