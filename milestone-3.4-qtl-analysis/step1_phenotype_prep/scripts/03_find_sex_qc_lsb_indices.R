#!/usr/bin/env Rscript

# Find pseudobulk_inputs.csv row indices for datasets containing samples that
# fail final sex QC. The reported config row numbers are the LSB_JOBINDEX values
# used by 03_submit_sex_qc_pass_residual_jobs.sh.
#
# Matching strategy:
#   sex_qc$file_region identifies the source dataset (for example GEN_A4_MEC).
#   Its GEN project component is matched to config$project, and any region
#   suffix is matched to config$default_region. project_for_plot is deliberately
#   not used because it can refer to a different plotting dataset.
#
# Usage:
#   Rscript 03_find_sex_qc_lsb_indices.R \
#     <sex_qc_samples.csv> <pseudobulk_inputs.csv>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript 03_find_sex_qc_lsb_indices.R ",
    "<sex_qc_samples.csv> <pseudobulk_inputs.csv>"
  )
}

sex_qc_file <- args[1]
config_file <- args[2]

sex_qc <- read.csv(
  sex_qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

config <- read.csv(
  config_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

required_sex_qc <- c("file_region", "individualID", "final_sex_qc")
required_config <- c("project", "default_region", "output_dir")

missing_sex_qc <- setdiff(required_sex_qc, colnames(sex_qc))
missing_config <- setdiff(required_config, colnames(config))

if (length(missing_sex_qc) > 0) {
  stop(
    "Sex-QC file is missing columns: ",
    paste(missing_sex_qc, collapse = ", ")
  )
}

if (length(missing_config) > 0) {
  stop(
    "Config file is missing columns: ",
    paste(missing_config, collapse = ", ")
  )
}

outlier_statuses <- c(
  "snRNAseq_QC_outlier",
  "RNA_sex_QC_outlier",
  "WGS_sex_QC_outlier",
  "RNA_and_WGS_sex_QC_outlier"
)

outliers <- sex_qc[
  sex_qc$final_sex_qc %in% outlier_statuses &
    !is.na(sex_qc$file_region) &
    trimws(sex_qc$file_region) != "",
  ,
  drop = FALSE
]

if (nrow(outliers) == 0) {
  stop("No final sex-QC outliers were found.")
}

outliers$file_region <- trimws(outliers$file_region)
config$project <- trimws(config$project)
config$default_region <- trimws(config$default_region)
config$lsb_jobindex <- seq_len(nrow(config))

match_file_region <- function(file_region, config) {
  # Prefer the longest project name that prefixes file_region.
  project_candidates <- unique(
    config$project[
      file_region == config$project |
        startsWith(file_region, paste0(config$project, "_"))
    ]
  )

  if (length(project_candidates) == 0) {
    stop("No config project matches file_region: ", file_region)
  }

  project_name <- project_candidates[
    which.max(nchar(project_candidates))
  ]

  candidate_rows <- which(config$project == project_name)
  region_suffix <- sub(
    paste0("^", project_name, "_?"),
    "",
    file_region
  )

  if (region_suffix != "") {
    region_rows <- candidate_rows[
      config$default_region[candidate_rows] == region_suffix
    ]

    if (length(region_rows) != 1) {
      stop(
        "Expected one config row for ",
        file_region,
        " using project=",
        project_name,
        " and default_region=",
        region_suffix,
        "; found ",
        length(region_rows)
      )
    }

    return(region_rows)
  }

  if (length(candidate_rows) != 1) {
    stop(
      "Ambiguous config match for file_region ",
      file_region,
      ": project ",
      project_name,
      " has ",
      length(candidate_rows),
      " config rows and file_region has no region suffix."
    )
  }

  candidate_rows
}

affected_file_regions <- sort(unique(outliers$file_region))
matched_rows <- vapply(
  affected_file_regions,
  match_file_region,
  integer(1),
  config = config
)

outlier_counts <- aggregate(
  x = list(n_outlier_individuals = outliers$individualID),
  by = list(file_region = outliers$file_region),
  FUN = function(x) length(unique(x))
)

result <- data.frame(
  file_region = affected_file_regions,
  lsb_jobindex = matched_rows,
  project = config$project[matched_rows],
  cohort = if ("cohort" %in% colnames(config)) {
    config$cohort[matched_rows]
  } else {
    NA_character_
  },
  default_region = config$default_region[matched_rows],
  output_dir = config$output_dir[matched_rows],
  stringsAsFactors = FALSE
)

result <- merge(
  result,
  outlier_counts,
  by = "file_region",
  all.x = TRUE,
  sort = FALSE
)
result <- result[order(result$lsb_jobindex), , drop = FALSE]

indices <- sort(unique(result$lsb_jobindex))
array_expression <- paste(indices, collapse = ",")

cat("\nDatasets containing final sex-QC outliers:\n\n")
print(result, row.names = FALSE)

cat("\nLSB indices:\n")
cat(array_expression, "\n")

cat("\nBSUB array directive:\n")
cat('#BSUB -J "SEXQC[', array_expression, ']"\n', sep = "")

