# Create pseudobulk expression and residual matrices from one harmonized h5ad.
#
# Usage:
#   Rscript 01_create_pseudobulk_residuals_from_h5ad.R <config_csv> <row_index>
#
# Required inputs:
#   1. config_csv
#      One row per h5ad dataset/region. The selected row gives the h5ad path,
#      donor/sample column, region column/default region, metadata cohort, and
#      output directory.
#
#   2. samples_single_cell_primary_GT.csv
#      Row-level single-cell metadata with curated WGS selection. This file is
#      produced outside this script and must include primary_genotype, which is
#      the chosen WGS ID per single-cell sample.
#
# Pipeline:
#   1. Read the selected config row and load the h5ad.
#   2. Read samples_single_cell_primary_GT.csv as character data so IDs such as
#      0054 are not converted to 54.
#   3. Join h5ad cell metadata to the primary-GT metadata. The script compares
#      possible metadata ID columns and uses the one with the best h5ad overlap.
#   4. Set output_id = primary_genotype.
#   5. Keep only cells with metadata, sex, age, and primary_genotype/WGS ID.
#   6. Aggregate counts to pseudobulk by donor and cell subclass/class.
#   7. Keep pseudobulk column names as the aggregateToPseudoBulk() sample IDs
#      while running dreamlet. dreamlet uses these IDs to match internal cell
#      counts, so changing them before processAssays() breaks bookkeeping.
#   8. Run processAssays(), detect outlier donors, remove outliers, rerun
#      processAssays()/dreamlet, and write residual matrices with final column
#      names converted to primary genotype IDs.
#
# Main outputs under output_dir:
#   pseudobulk/      Pseudobulk RDS files with original aggregateToPseudoBulk()
#                    sample IDs and primary genotype stored in colData$output_id.
#   processassays/   dreamlet processAssays RDS files.
#   residuals/       Standard and Pearson residual matrices, gzipped TSV.
#   outliers/        Outlier scores, summaries, and exclude ID lists.
#   sex_qc/          XIST/UTY sex QC plots.
#   id_maps/         Mapping from original h5ad donor/sample IDs to primary GT.
#   wgs_qc/          Primary genotype/WGS resolution reports.

.libPaths(c("/sc/arion/projects/psychAD/aging/kiran/RLib_4_4.1", .libPaths()))

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(HDF5Array)
  library(S4Vectors)
  library(dreamlet)
  library(crumblr)
  library(RNOmni)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(R.utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_create_pseudobulk_residuals_from_h5ad.R <config_csv> <row_index>")
}

config_csv <- args[1]
row_index <- as.integer(args[2])

METADATA_PATH <- "/sc/arion/projects/CommonMind/genesis/metadata/outputs/samples_single_cell_primary_GT.csv"

config_all <- read.csv(config_csv, stringsAsFactors = FALSE)
if (row_index < 1 || row_index > nrow(config_all)) {
  stop("row_index out of range: ", row_index, " for config with ", nrow(config_all), " rows")
}

config <- config_all[row_index, ]

PROJECT <- config$project
COHORT <- config$cohort
h5ad_file <- config$h5ad_path
col_donor <- config$donor_col
col_region <- config$region_col
default_region <- config$default_region
OUT <- config$output_dir

metadata_cohort <- if ("metadata_cohort" %in% colnames(config_all)) config$metadata_cohort else ""
metadata_cohort <- ifelse(is.na(metadata_cohort), "", metadata_cohort)

cluster_id_options <- c("class", "subclass")
fdr_cutoff <- 1e-5
nPC <- 2
min_donors_for_dreamlet <- 10

form <- ~ sex + scale(age) +
  log(n_genes) + percent_mito +
  mito_genes + mito_ribo + ribo_genes

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir_pb <- file.path(OUT, "pseudobulk"); dir.create(dir_pb, recursive = TRUE, showWarnings = FALSE)
dir_pa <- file.path(OUT, "processassays"); dir.create(dir_pa, recursive = TRUE, showWarnings = FALSE)
dir_resid <- file.path(OUT, "residuals"); dir.create(dir_resid, recursive = TRUE, showWarnings = FALSE)
dir_outl <- file.path(OUT, "outliers"); dir.create(dir_outl, recursive = TRUE, showWarnings = FALSE)
dir_qc <- file.path(OUT, "sex_qc"); dir.create(dir_qc, recursive = TRUE, showWarnings = FALSE)
dir_idmap <- file.path(OUT, "id_maps"); dir.create(dir_idmap, recursive = TRUE, showWarnings = FALSE)
dir_wgs <- file.path(OUT, "wgs_qc"); dir.create(dir_wgs, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else as.character(x[1])
}

make_prefix <- function(base_dir, project, cohort, donor_col, region, cluster_id, suffix = "") {
  file.path(
    base_dir,
    paste0(
      safe_name(project), "_",
      safe_name(cohort), "_",
      safe_name(donor_col), "_",
      safe_name(region), "_",
      safe_name(cluster_id),
      suffix
    )
  )
}

resolve_wgs_ids <- function(df) {
  wgs_cols <- grep("^WGS_id_", colnames(df), value = TRUE)
  if (length(wgs_cols) == 0) {
    return(data.frame(
      output_id = rep(NA_character_, nrow(df)),
      n_wgs_ids = rep(0L, nrow(df)),
      wgs_columns_present = rep("", nrow(df)),
      wgs_values_present = rep("", nrow(df)),
      stringsAsFactors = FALSE
    ))
  }

  output_id <- character(nrow(df))
  n_wgs_ids <- integer(nrow(df))
  wgs_columns_present <- character(nrow(df))
  wgs_values_present <- character(nrow(df))

  for (i in seq_len(nrow(df))) {
    vals <- as.character(df[i, wgs_cols, drop = TRUE])
    present <- !is.na(vals) & vals != ""
    n_wgs_ids[i] <- sum(present)
    wgs_columns_present[i] <- paste(wgs_cols[present], collapse = ";")
    wgs_values_present[i] <- paste(vals[present], collapse = ";")
    output_id[i] <- if (any(present)) vals[which(present)[1]] else NA_character_
  }

  data.frame(
    output_id = output_id,
    n_wgs_ids = n_wgs_ids,
    wgs_columns_present = wgs_columns_present,
    wgs_values_present = wgs_values_present,
    stringsAsFactors = FALSE
  )
}

summarize_wgs_sources <- function(df, output_id) {
  wgs_cols <- grep("^WGS_id_", colnames(df), value = TRUE)
  if (length(wgs_cols) == 0) {
    return(data.frame(wgs_column = character(), n_cells = integer(), n_donors = integer()))
  }

  out <- lapply(wgs_cols, function(col) {
    present <- !is.na(df[[col]]) & df[[col]] != ""
    data.frame(
      wgs_column = col,
      n_cells = sum(present),
      n_donors = length(unique(df[[col]][present])),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

rename_cols_to_output_id <- function(mat, id_map) {
  mat <- as.matrix(mat)

  if (is.null(colnames(mat)) && !is.null(rownames(mat))) {
    row_idx <- match(rownames(mat), id_map$pb_colname)
    if (any(!is.na(row_idx))) {
      mat <- t(mat)
    }
  }

  if (is.null(colnames(mat))) {
    stop("Residual matrix has no column names.")
  }

  idx <- match(colnames(mat), id_map$pb_colname)
  if (any(is.na(idx)) && !is.null(rownames(mat))) {
    row_idx <- match(rownames(mat), id_map$pb_colname)
    if (sum(!is.na(row_idx)) > sum(!is.na(idx))) {
      mat <- t(mat)
      idx <- match(colnames(mat), id_map$pb_colname)
    }
  }

  if (any(is.na(idx))) {
    stop("Missing output_id for columns: ", paste(colnames(mat)[is.na(idx)], collapse = ", "))
  }

  new_names <- id_map$output_id[idx]
  if (any(is.na(new_names)) || any(new_names == "")) {
    stop("Some output IDs are missing.")
  }
  if (anyDuplicated(new_names)) {
    stop("Duplicate output IDs after renaming.")
  }

  colnames(mat) <- new_names
  mat
}

drop_aggr_mean_vars_from_coldata <- function(pb) {
  # These variables are cell-type-specific after aggregateToPseudoBulk() and
  # live in metadata(pb)$aggr_means. Leaving donor-level copies in colData(pb)
  # can confuse dreamlet formula resolution.
  vars <- c("n_genes", "percent_mito", "mito_genes", "mito_ribo", "ribo_genes")
  keep <- setdiff(colnames(colData(pb)), vars)
  colData(pb) <- colData(pb)[, keep, drop = FALSE]
  pb
}

plot_xist_uty <- function(pb, project, cohort, region, cluster_id, out_path) {
  sex_expr <- lapply(assayNames(pb), function(cell_type) {
    mat <- assay(pb, cell_type)

    xist <- if ("XIST" %in% rownames(mat)) as.numeric(mat["XIST", ]) else rep(NA_real_, ncol(mat))
    uty <- if ("UTY" %in% rownames(mat)) as.numeric(mat["UTY", ]) else rep(NA_real_, ncol(mat))
    sex <- if ("sex" %in% colnames(colData(pb))) as.character(colData(pb)$sex) else rep(NA_character_, ncol(pb))

    data.frame(
      donor_id = colnames(mat),
      cell_type = cell_type,
      sex = sex,
      XIST = xist,
      UTY = uty,
      stringsAsFactors = FALSE
    )
  })

  sex_expr <- do.call(rbind, sex_expr)

  p <- ggplot(sex_expr, aes(x = XIST, y = UTY, color = sex)) +
    geom_point(alpha = 0.75, size = 1.8) +
    facet_wrap(~ cell_type, scales = "free") +
    theme_classic() +
    labs(
      title = paste(project, cohort, region, cluster_id, "XIST vs UTY"),
      x = "XIST pseudobulk counts",
      y = "UTY pseudobulk counts",
      color = "Metadata sex"
    )

  ggsave(out_path, p, width = 12, height = 8, bg = "white")
}

print_processassays_details <- function(res.proc, label) {
  cat("\nprocessAssays details:", label, "\n")
  det <- tryCatch(
    details(res.proc),
    error = function(e) {
      cat("details(res.proc) failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(det)) {
    print(det)
  }
  cat("\n")
}

filter_assays_by_n_retain <- function(res.proc, min_donors, report_path, label) {
  det <- tryCatch(
    as.data.frame(details(res.proc)),
    error = function(e) {
      cat("Could not inspect processAssays details for donor-count filter:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(det) || !"assay" %in% colnames(det) || !"n_retain" %in% colnames(det)) {
    cat("Skipping donor-count filter because details(res.proc) did not return assay/n_retain.\n")
    return(res.proc)
  }

  det$keep_for_dreamlet <- det$n_retain >= min_donors
  det$reason <- ifelse(
    det$keep_for_dreamlet,
    "retained",
    paste0("n_retain < ", min_donors)
  )

  write.csv(det, report_path, row.names = FALSE)

  skipped <- det$assay[!det$keep_for_dreamlet]
  if (length(skipped) > 0) {
    cat("\nSkipping assays before dreamlet:", label, "\n")
    cat("Minimum retained donors:", min_donors, "\n")
    print(det[!det$keep_for_dreamlet, c("assay", "n_retain", "reason"), drop = FALSE])
  }

  keep_assays <- det$assay[det$keep_for_dreamlet]
  if (length(keep_assays) == 0) {
    stop("No assays remain after donor-count filter for ", label)
  }

  keep_idx <- match(keep_assays, assayNames(res.proc))
  res.proc[keep_idx]
}

cat("\n============================\n")
cat("Row:", row_index, "\n")
cat("Project:", PROJECT, "\n")
cat("Cohort:", COHORT, "\n")
cat("H5AD:", h5ad_file, "\n")
cat("Donor column:", col_donor, "\n")
cat("Region column:", col_region, "\n")
cat("Default region:", default_region, "\n")
cat("Metadata cohort:", metadata_cohort, "\n")
cat("Output:", OUT, "\n")
cat("============================\n\n")

sample_meta <- read.csv(
  METADATA_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)
metadata_sample_col <- if ("sample_id" %in% colnames(sample_meta)) {
  "sample_id"
} else if ("genesis_sample" %in% colnames(sample_meta)) {
  "genesis_sample"
} else {
  stop("Metadata must contain either sample_id or genesis_sample: ", METADATA_PATH)
}

if (!"primary_genotype" %in% colnames(sample_meta)) {
  stop("Metadata must contain primary_genotype: ", METADATA_PATH)
}

if (metadata_cohort != "" && "cohort" %in% colnames(sample_meta)) {
  sample_meta <- sample_meta[sample_meta$cohort == metadata_cohort, , drop = FALSE]
}
sample_meta$metadata_donor_id <- if ("donor_id" %in% colnames(sample_meta)) {
  sample_meta$donor_id
} else {
  NA_character_
}
sample_meta <- sample_meta[!duplicated(sample_meta[[metadata_sample_col]]), , drop = FALSE]
cat("Metadata path:", METADATA_PATH, "\n")
cat("Metadata sample ID column:", metadata_sample_col, "\n")
cat("Metadata rows after cohort filter:", nrow(sample_meta), "\n")
cat("Metadata rows with primary_genotype:", sum(!is.na(sample_meta$primary_genotype) & sample_meta$primary_genotype != ""), "\n")

sce <- readH5AD(h5ad_file, use_hdf5 = TRUE, verbose = TRUE)
cat("Loaded:", ncol(sce), "cells,", nrow(sce), "genes\n")

if (!col_donor %in% colnames(colData(sce))) {
  stop("Missing donor column in colData: ", col_donor)
}

if (is.na(col_region) || col_region == "" || !col_region %in% colnames(colData(sce))) {
  cat("Region column missing. Setting brain_region =", default_region, "\n")
  colData(sce)$brain_region <- default_region
  col_region <- "brain_region"
} else {
  region_values <- as.character(colData(sce)[[col_region]])
  region_values[is.na(region_values) | region_values == ""] <- default_region
  colData(sce)[[col_region]] <- region_values
}

if (!"counts" %in% assayNames(sce)) {
  cat("No counts assay found. Loading raw/X as counts...\n")
  raw_counts <- HDF5Array::H5SparseMatrix(h5ad_file, "raw/X")
  dimnames(raw_counts) <- dimnames(sce)
  assay(sce, "counts", withDimnames = FALSE) <- raw_counts
}

cd <- as.data.frame(colData(sce))
cd$cell_barcode <- rownames(cd)

metadata_join_candidates <- unique(c(metadata_sample_col, "donor_id"))
metadata_join_candidates <- metadata_join_candidates[metadata_join_candidates %in% colnames(sample_meta)]
h5ad_ids <- unique(as.character(cd[[col_donor]]))
join_overlap <- vapply(metadata_join_candidates, function(join_col) {
  sum(h5ad_ids %in% unique(as.character(sample_meta[[join_col]])))
}, integer(1))
metadata_join_col <- metadata_join_candidates[which.max(join_overlap)]

cat("Metadata join candidates:\n")
print(data.frame(join_col = metadata_join_candidates, n_h5ad_ids_matched = join_overlap))
cat("Using metadata join column:", metadata_join_col, "\n")

sample_meta_join <- sample_meta[!duplicated(sample_meta[[metadata_join_col]]), , drop = FALSE]

cd2 <- merge(
  cd,
  sample_meta_join,
  by.x = col_donor,
  by.y = metadata_join_col,
  all.x = TRUE,
  suffixes = c("", "_meta")
)

cd2 <- cd2[match(rownames(colData(sce)), cd2$cell_barcode), , drop = FALSE]
rownames(cd2) <- cd2$cell_barcode

# Output IDs are curated outside this script in samples_single_cell_primary_GT.csv.
# primary_genotype is the selected WGS ID with the highest available WGS coverage.
cd2$output_id <- as.character(cd2$primary_genotype)
cd2$output_id[is.na(cd2$output_id) | cd2$output_id == ""] <- NA_character_

wgs_cols <- grep("^WGS_id_", colnames(cd2), value = TRUE)
if (length(wgs_cols) > 0) {
  wgs_resolved <- resolve_wgs_ids(cd2)
  cd2$n_wgs_ids <- wgs_resolved$n_wgs_ids
  cd2$wgs_columns_present <- wgs_resolved$wgs_columns_present
  cd2$wgs_values_present <- wgs_resolved$wgs_values_present
} else {
  cd2$n_wgs_ids <- 0L
  cd2$wgs_columns_present <- ""
  cd2$wgs_values_present <- ""
}

metadata_missing <- is.na(cd2$metadata_donor_id) | cd2$metadata_donor_id == ""
cat("Cells missing metadata:", sum(metadata_missing), "\n")
cat("Donors in H5AD:", length(unique(cd2[[col_donor]])), "\n")
cat("Donors missing metadata:", length(unique(cd2[[col_donor]][metadata_missing])), "\n")
cat("Donors missing primary_genotype/output_id:", length(unique(cd2[[col_donor]][is.na(cd2$output_id) | cd2$output_id == ""])), "\n")

wgs_summary <- summarize_wgs_sources(cd2, cd2$output_id)
if (nrow(wgs_summary) > 0) {
  cat("WGS source columns after metadata merge:\n")
  print(wgs_summary)
}

report_cols <- c(
  col_donor,
  "metadata_donor_id",
  "donor_id",
  "sex",
  "age",
  "output_id",
  "primary_genotype",
  "primary_genotype_coverage",
  "primary_genotype_status",
  "WGS_coverage",
  "n_wgs_ids",
  "wgs_columns_present",
  "wgs_values_present"
)
report_cols <- report_cols[report_cols %in% colnames(cd2)]
wgs_report <- unique(cd2[, report_cols, drop = FALSE])
colnames(wgs_report)[1] <- "pb_colname"
wgs_report <- wgs_report[order(wgs_report$n_wgs_ids, wgs_report$pb_colname), , drop = FALSE]

wgs_report_path <- file.path(
  dir_wgs,
  paste0(safe_name(PROJECT), "_", safe_name(COHORT), "_", safe_name(default_region), "_wgs_id_resolution_report.tsv")
)
write.table(wgs_report, file = wgs_report_path, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Primary genotype report:", wgs_report_path, "\n")
cat("Donors with 0 WGS IDs:", sum(wgs_report$n_wgs_ids == 0), "\n")
cat("Donors with >1 WGS IDs:", sum(wgs_report$n_wgs_ids > 1), "\n")
if ("primary_genotype_status" %in% colnames(wgs_report)) {
  cat("Primary genotype status:\n")
  print(table(wgs_report$primary_genotype_status, useNA = "ifany"))
}

keep <- !is.na(cd2$sex) &
  !is.na(cd2$age) &
  !is.na(cd2$output_id) &
  cd2$output_id != ""

cat("Cells retained after metadata/WGS filter:", sum(keep), "of", length(keep), "\n")
sce <- sce[, keep]
colData(sce) <- S4Vectors::DataFrame(cd2[keep, , drop = FALSE])
colData(sce)$sex <- factor(as.character(colData(sce)$sex))
colData(sce)$age <- as.numeric(colData(sce)$age)

region_options <- sort(unique(as.character(colData(sce)[[col_region]])))
cat("Regions:", paste(region_options, collapse = ", "), "\n")

# =============================================================================
# STEP 1: PSEUDOBULK + PROCESSASSAYS + SEX QC PLOT
# =============================================================================
cat("\n========== STEP 1: Pseudobulk + processAssays + sex QC ==========\n")

for (region in region_options) {
  for (cluster_id in cluster_id_options) {
    cat("\n---", region, cluster_id, "---\n")

    if (!cluster_id %in% colnames(colData(sce))) {
      warning("Missing cluster column: ", cluster_id, ". Skipping.")
      next
    }

    sce_region <- sce[, as.character(colData(sce)[[col_region]]) == region]
    cat("Cells:", ncol(sce_region), "\n")
    cat("Donors:", length(unique(colData(sce_region)[[col_donor]])), "\n")

    pb <- aggregateToPseudoBulk(
      sce_region,
      assay = "counts",
      cluster_id = cluster_id,
      sample_id = col_donor
    )

    meta_region <- as.data.frame(colData(sce_region))
    donor_meta <- meta_region |>
      group_by(.data[[col_donor]]) |>
      summarize(
        sex = first_non_na(sex),
        age = median(as.numeric(age), na.rm = TRUE),
        output_id = first_non_na(output_id),
        .groups = "drop"
      )

    donor_meta <- as.data.frame(donor_meta)
    rownames(donor_meta) <- donor_meta[[col_donor]]
    donor_meta <- donor_meta[colnames(pb), , drop = FALSE]

    colData(pb)$sex <- factor(donor_meta$sex)
    colData(pb)$age <- donor_meta$age
    colData(pb)$output_id <- donor_meta$output_id
    original_pb_colnames <- colnames(pb)
    primary_genotype_colnames <- as.character(colData(pb)$output_id)

    if (any(is.na(primary_genotype_colnames)) || any(primary_genotype_colnames == "")) {
      stop("Some primary genotype IDs are missing in ", PROJECT, " ", COHORT, " ", region)
    }
    if (anyDuplicated(primary_genotype_colnames)) {
      stop("Duplicate primary genotype IDs in ", PROJECT, " ", COHORT, " ", region)
    }

    prefix_pb <- make_prefix(dir_pb, PROJECT, COHORT, col_donor, region, cluster_id, "_PB")
    prefix_pa <- make_prefix(dir_pa, PROJECT, COHORT, col_donor, region, cluster_id, "_processAssays")
    prefix_qc <- make_prefix(dir_qc, PROJECT, COHORT, col_donor, region, cluster_id, "_XIST_vs_UTY_by_sex.pdf")
    prefix_id <- make_prefix(dir_idmap, PROJECT, COHORT, col_donor, region, cluster_id, "_id_map")

    id_map <- data.frame(
      pb_colname = colnames(pb),
      donor_id = original_pb_colnames,
      output_id = as.character(colData(pb)$output_id),
      project = PROJECT,
      cohort = COHORT,
      donor_col = col_donor,
      brain_region = region,
      cluster_id = cluster_id,
      stringsAsFactors = FALSE
    )

    saveRDS(pb, paste0(prefix_pb, ".RDS"))
    write.csv(id_map, paste0(prefix_id, ".csv"), row.names = FALSE, quote = FALSE)

    plot_xist_uty(pb, PROJECT, COHORT, region, cluster_id, prefix_qc)

    pb <- drop_aggr_mean_vars_from_coldata(pb)
    res.proc <- processAssays(pb, form)
    print_processassays_details(
      res.proc,
      paste(PROJECT, COHORT, region, cluster_id, "STEP1", sep = " | ")
    )
    saveRDS(res.proc, paste0(prefix_pa, ".RDS"))

    cat("Saved:", paste0(prefix_pb, ".RDS"), "\n")
    cat("Saved:", paste0(prefix_pa, ".RDS"), "\n")
    cat("Saved:", paste0(prefix_id, ".csv"), "\n")
    cat("Saved:", prefix_qc, "\n")

    rm(sce_region, pb, res.proc, id_map, donor_meta, meta_region)
    gc()
  }
}

# =============================================================================
# STEP 2: OUTLIER DETECTION
# =============================================================================
cat("\n========== STEP 2: Outlier Detection ==========\n")

for (region in region_options) {
  cluster_id <- "class"

  cat("\n---", region, cluster_id, "---\n")

  prefix_pb <- make_prefix(dir_pb, PROJECT, COHORT, col_donor, region, cluster_id, "_PB")
  prefix_outl <- make_prefix(dir_outl, PROJECT, COHORT, col_donor, region, cluster_id)

  pb <- readRDS(paste0(prefix_pb, ".RDS"))
  pb <- drop_aggr_mean_vars_from_coldata(pb)

  cat("Fitting dreamlet for outlier residuals...\n")
  res.proc.qc <- processAssays(pb, form)
  print_processassays_details(
    res.proc.qc,
    paste(PROJECT, COHORT, region, cluster_id, "STEP2_OUTLIER_QC", sep = " | ")
  )
  fit.qc <- dreamlet(res.proc.qc, form)
  residList <- residuals(fit.qc)
  rm(fit.qc, res.proc.qc)
  gc()

  residList_rn <- lapply(residList, function(Y) t(apply(Y, 1, RankNorm)))
  df_scores <- as.data.frame(outlierByAssay(residList_rn, nPC = nPC))

  df_outliers <- dplyr::as_tibble(df_scores) |>
    dplyr::group_by(ID) |>
    dplyr::summarize(
      n = dplyr::n(),
      x = sum(chisq, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pValue = stats::pchisq(x, n * nPC, lower.tail = FALSE),
      FDR = stats::p.adjust(pValue, method = "fdr")
    ) |>
    dplyr::arrange(pValue)

  excludeIDs <- df_outliers |>
    dplyr::filter(FDR < fdr_cutoff) |>
    dplyr::pull(ID)

  df_scores$excludeIDs <- df_scores$ID %in% excludeIDs

  cat("Outliers (FDR <", fdr_cutoff, "):", length(excludeIDs), "\n")
  if (length(excludeIDs) > 0) print(excludeIDs)

  p <- ggplot(df_scores, aes(PC1, PC2, colour = excludeIDs)) +
    geom_point(size = 1.5, alpha = 0.7) +
    ggrepel::geom_text_repel(
      data = subset(df_scores, excludeIDs),
      aes(label = ID),
      size = 2.5,
      max.overlaps = 20,
      segment.colour = "grey50"
    ) +
    scale_colour_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
      labels = c("FALSE" = "Retained", "TRUE" = "Outlier"),
      name = NULL
    ) +
    facet_wrap(~ assay, nrow = 2, scales = "free") +
    labs(
      title = paste0("Outlier PCA | ", PROJECT, " | ", COHORT, " | ", region, " | ", cluster_id),
      caption = paste0(
        "n outliers FDR < ", fdr_cutoff, " = ", length(excludeIDs),
        " | n donors = ", length(unique(df_scores$ID))
      )
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.caption = element_text(size = 7, colour = "grey40"),
      legend.position = "top"
    )

  ggsave(paste0(prefix_outl, "_outlier_PCA.pdf"), plot = p, width = 14, height = 10, bg = "white")
  write.csv(df_scores, paste0(prefix_outl, "_scores.csv"), row.names = FALSE)
  write.csv(df_outliers, paste0(prefix_outl, "_outlier_summary.csv"), row.names = FALSE)
  writeLines(as.character(excludeIDs), paste0(prefix_outl, "_excludeIDs.txt"))

  rm(pb, residList, residList_rn, df_scores, df_outliers, p)
  gc()
}

# =============================================================================
# STEP 3: FINAL PROCESSASSAYS + RESIDUALS AFTER OUTLIER REMOVAL
# =============================================================================
cat("\n========== STEP 3: Final processAssays + residuals ==========\n")

for (region in region_options) {
  for (cluster_id in cluster_id_options) {
    cat("\n---", region, cluster_id, "---\n")

    prefix_pb <- make_prefix(dir_pb, PROJECT, COHORT, col_donor, region, cluster_id, "_PB")
    prefix_outl_class <- make_prefix(dir_outl, PROJECT, COHORT, col_donor, region, "class")
    prefix_pa_final <- make_prefix(dir_pa, PROJECT, COHORT, col_donor, region, cluster_id, "_final_processAssays")
    prefix_id <- make_prefix(dir_idmap, PROJECT, COHORT, col_donor, region, cluster_id, "_id_map")
    prefix_filter <- make_prefix(dir_pa, PROJECT, COHORT, col_donor, region, cluster_id, "_dreamlet_assay_filter")

    pb <- readRDS(paste0(prefix_pb, ".RDS"))
    pb <- drop_aggr_mean_vars_from_coldata(pb)
    id_map <- read.csv(paste0(prefix_id, ".csv"), stringsAsFactors = FALSE)

    exclude_path <- paste0(prefix_outl_class, "_excludeIDs.txt")
    excludeIDs <- if (file.exists(exclude_path)) {
      scan(exclude_path, what = character(), quiet = TRUE)
    } else {
      character(0)
    }
    excludeIDs <- excludeIDs[nchar(trimws(excludeIDs)) > 0]

    cat("Class-level outliers applied to", cluster_id, ":", length(excludeIDs), "\n")
    if (length(excludeIDs) > 0) {
      pb <- pb[, !colnames(pb) %in% excludeIDs]
    }

    id_map_clean <- id_map[id_map$pb_colname %in% colnames(pb), , drop = FALSE]
    id_map_clean <- id_map_clean[match(colnames(pb), id_map_clean$pb_colname), , drop = FALSE]
    stopifnot(nrow(id_map_clean) == ncol(pb))

    res.proc <- processAssays(pb, form)
    print_processassays_details(
      res.proc,
      paste(PROJECT, COHORT, region, cluster_id, "STEP3_FINAL", sep = " | ")
    )
    res.proc <- filter_assays_by_n_retain(
      res.proc,
      min_donors_for_dreamlet,
      paste0(prefix_filter, ".csv"),
      paste(PROJECT, COHORT, region, cluster_id, sep = " | ")
    )
    saveRDS(res.proc, paste0(prefix_pa_final, ".RDS"))

    fit <- dreamlet(res.proc, form)

    for (CT in assayNames(res.proc)) {
      cat("Cell type:", CT, "\n")

      prefix_ct <- make_prefix(dir_resid, PROJECT, COHORT, col_donor, region, cluster_id, paste0("_", CT))

      data_std <- residuals(fit[[CT]])
      data_std <- rename_cols_to_output_id(data_std, id_map_clean)
      std_file <- paste0(prefix_ct, "_residuals.tsv")
      write.table(format(data_std, digits = 5), file = std_file, quote = FALSE, sep = "\t")
      gzip(std_file, overwrite = TRUE)

      data_prs <- residuals(fit[[CT]], res.proc[[CT]], type = "pearson")
      data_prs <- rename_cols_to_output_id(data_prs, id_map_clean)
      prs_file <- paste0(prefix_ct, "_residualsPearson.tsv")
      write.table(format(data_prs, digits = 5), file = prs_file, quote = FALSE, sep = "\t")
      gzip(prs_file, overwrite = TRUE)

      gene_mean <- rowMeans(data_prs, na.rm = TRUE)
      gene_sd <- apply(data_prs, 1, sd, na.rm = TRUE)
      plot_df <- data.frame(gene = rownames(data_prs), mean = gene_mean, sd = gene_sd)
      sd_thresh <- quantile(gene_sd, 0.99, na.rm = TRUE)
      plot_df$label <- ifelse(plot_df$sd >= sd_thresh, plot_df$gene, NA)
      plot_df$is_high <- !is.na(plot_df$label)

      p_sdmean <- ggplot(plot_df, aes(x = mean, y = sqrt(sd))) +
        geom_point(aes(colour = is_high), size = 0.8, alpha = 0.6) +
        geom_smooth(method = "loess", se = TRUE, colour = "firebrick", linewidth = 0.8, fill = "firebrick", alpha = 0.15) +
        geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 20, segment.colour = "grey50", na.rm = TRUE) +
        scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkorange"), name = NULL) +
        labs(
          title = paste0("sqrt(SD) vs Mean -- Pearson Residuals\n", PROJECT, " | ", COHORT, " | ", region, " | ", CT),
          x = "Mean across donors (WGS/output_id)",
          y = "sqrt(SD) across donors (WGS/output_id)"
        ) +
        theme_classic()

      ggsave(paste0(prefix_ct, "_PearsonResid_SDvsMean.png"), plot = p_sdmean, width = 7, height = 6, dpi = 180, bg = "white")

      rm(data_std, data_prs, gene_mean, gene_sd, plot_df, p_sdmean)
      gc()
    }

    rm(pb, id_map, id_map_clean, res.proc, fit)
    gc()
  }
}

cat("\nPipeline complete:", PROJECT, COHORT, "\n")
