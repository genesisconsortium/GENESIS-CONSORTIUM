(scanpy310) [girdhk01@li04e04 scripts]$ cat create_samples_single_cell_primary_GT.R
#!/usr/bin/env Rscript

# Create samples_single_cell_primary_GT.csv.
#
# For each row in samples_single_cell.csv:
#   1. Read WGS IDs from all WGS_id_* columns, in column order.
#   2. Look up wgs_MEAN_COVERAGE from Picard technical covariates.
#   3. Write WGS_coverage in the same order as the WGS IDs appear.
#   4. Write primary_genotype as the WGS ID with the highest coverage.
#   5. If a WGS ID has no Picard row, its coverage is NA.

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- "/sc/arion/projects/CommonMind/genesis/metadata/outputs/samples_single_cell.csv"
genotype_root <- "/sc/arion/projects/CommonMind/genesis/Genotype_files"
output_file <- "/sc/arion/projects/CommonMind/genesis/metadata/outputs/samples_single_cell_primary_GT.csv"

if (length(args) >= 1) {
  metadata_file <- args[1]
}
if (length(args) >= 2) {
  genotype_root <- args[2]
}
if (length(args) >= 3) {
  output_file <- args[3]
}

message("Metadata file: ", metadata_file)
message("Genotype root: ", genotype_root)
message("Output file: ", output_file)

samples <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

technical_files <- data.frame(
  project_id = c(
    "CMC_NYGC",
    "CMC_Novogene",
    "AMP_PD",
    "GEN_A2_ROSMAP_DeJager",
    "ROSMAP_RUSH_RADC_Diversity",
    "ROSMAP_RUSH_Diversity_360Samples",
    "NDA_other",
    "BD2",
    "NDA_MSSM",
    "SEAAD"
  ),
  wgs_column = c(
    "WGS_id_CMC",
    "WGS_id_CMC",
    "WGS_id_AMPPD",
    "WGS_id_ROSMAP",
    "WGS_id_ROSMAP_DIVERSITY",
    "WGS_id_ROSMAP_DIVERSITY",
    "WGS_id_NDA",
    "WGS_id_NDA",
    "WGS_id_NDA",
    "WGS_id_SEAAD"
  ),
  relative_path = c(
    "CommonMind/CommonMind_TRs/CMC_NYGC/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "CommonMind/CommonMind_TRs/CMC_Novogene/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "AMP_PD/AMP_PD_TRs/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "GEN_A2_ROSMAP_DeJager/GEN_A2_TRs/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "ROSMAP_RUSH_RADC_Diversity/ROSMAP_RUSH_RADC_Diversity_TRs/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "ROSMAP_RUSH_RADC_Diversity/ROSMAP_RUSH_Diversity_360Samples/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "NDA_other/NDA_other_TRs/Processed_TRs/covariates/technical/NDA_other_picard_technical_covariates.tsv",
    "BD2/BD2_TRs/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "NDA_MSSM/NDA_MSSM_TRs/Processed_TRs/covariates/technical/picard_technical_covariates.tsv",
    "SEA-AD/SEAAD_TRs/Processed_TRs/covariates/technical/picard_technical_covariates.tsv"
  ),
  stringsAsFactors = FALSE
)

all_coverage <- data.frame()

get_wgs_family <- function(wgs_column_name) {
  if (grepl("^WGS_id_NDA", wgs_column_name)) {
    return("NDA")
  }
  if (grepl("^WGS_id_CMC", wgs_column_name)) {
    return("CMC")
  }
  if (grepl("^WGS_id_AMPPD", wgs_column_name)) {
    return("AMPPD")
  }
  if (grepl("^WGS_id_ROSMAP_DIVERSITY", wgs_column_name)) {
    return("ROSMAP_DIVERSITY")
  }
  if (grepl("^WGS_id_ROSMAP", wgs_column_name)) {
    return("ROSMAP")
  }
  if (grepl("^WGS_id_SEAAD", wgs_column_name)) {
    return("SEAAD")
  }
  return(wgs_column_name)
}

for (i in seq_len(nrow(technical_files))) {
  this_file <- file.path(genotype_root, technical_files$relative_path[i])
  this_project <- technical_files$project_id[i]
  this_wgs_column <- technical_files$wgs_column[i]
  this_wgs_family <- get_wgs_family(this_wgs_column)

  if (!file.exists(this_file)) {
    message("Missing Picard covariate file: ", this_file)
    next
  }

  message("Reading Picard covariates: ", this_file)
  covars <- read.delim(
    this_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )

  if (!"sample_id" %in% colnames(covars)) {
    message("  Skipping because sample_id column is missing.")
    next
  }

  if (!"wgs_MEAN_COVERAGE" %in% colnames(covars)) {
    message("  Skipping because wgs_MEAN_COVERAGE column is missing.")
    next
  }

  covars_small <- data.frame(
    WGS_id = as.character(covars$sample_id),
    WGS_coverage_value = as.numeric(covars$wgs_MEAN_COVERAGE),
    WGS_project = this_project,
    WGS_metadata_column = this_wgs_column,
    WGS_source_family = this_wgs_family,
    stringsAsFactors = FALSE
  )

  all_coverage <- rbind(all_coverage, covars_small)
}

if (nrow(all_coverage) == 0) {
  stop("No Picard technical covariates were loaded.")
}

# If the same WGS ID appears more than once for the same metadata WGS column,
# keep the highest coverage row.
all_coverage <- all_coverage[order(
  all_coverage$WGS_metadata_column,
  all_coverage$WGS_source_family,
  all_coverage$WGS_id,
  -all_coverage$WGS_coverage_value
), ]
all_coverage <- all_coverage[!duplicated(
  paste(all_coverage$WGS_source_family, all_coverage$WGS_id, sep = "||")
), ]

wgs_columns <- grep("^WGS_id_", colnames(samples), value = TRUE)

samples$WGS_coverage <- NA_character_
samples$primary_genotype <- NA_character_
samples$primary_genotype_coverage <- NA_real_
samples$n_WGS_ids <- NA_integer_
samples$n_WGS_ids_with_coverage <- NA_integer_
samples$primary_genotype_status <- NA_character_

for (row_index in seq_len(nrow(samples))) {
  row_wgs_ids <- as.character(samples[row_index, wgs_columns])
  row_wgs_ids[row_wgs_ids == ""] <- NA_character_
  row_wgs_ids[row_wgs_ids == "NA"] <- NA_character_

  row_coverages <- rep(NA_real_, length(row_wgs_ids))
  present_wgs_index <- which(!is.na(row_wgs_ids))
  samples$n_WGS_ids[row_index] <- length(present_wgs_index)

  for (wgs_index in seq_along(wgs_columns)) {
    this_wgs_column <- wgs_columns[wgs_index]
    this_wgs_family <- get_wgs_family(this_wgs_column)
    this_wgs_id <- row_wgs_ids[wgs_index]

    if (is.na(this_wgs_id)) {
      next
    }

    match_index <- which(
      all_coverage$WGS_source_family == this_wgs_family &
        all_coverage$WGS_id == this_wgs_id
    )

    if (length(match_index) == 0) {
      next
    }

    row_coverages[wgs_index] <- all_coverage$WGS_coverage_value[match_index[1]]
  }

  if (length(present_wgs_index) == 0) {
    samples$WGS_coverage[row_index] <- NA_character_
    samples$n_WGS_ids_with_coverage[row_index] <- 0L
    samples$primary_genotype_status[row_index] <- "no_wgs_id_in_metadata"
    next
  }

  present_coverages <- row_coverages[present_wgs_index]
  samples$n_WGS_ids_with_coverage[row_index] <- sum(!is.na(present_coverages))
  coverage_strings <- ifelse(is.na(present_coverages), "NA", as.character(present_coverages))
  samples$WGS_coverage[row_index] <- paste(coverage_strings, collapse = ",")

  if (all(is.na(present_coverages))) {
    samples$primary_genotype_status[row_index] <- "wgs_id_without_picard_coverage"
    next
  }

  best_present_index <- which.max(present_coverages)
  best_index <- present_wgs_index[best_present_index]
  samples$primary_genotype[row_index] <- row_wgs_ids[best_index]
  samples$primary_genotype_coverage[row_index] <- row_coverages[best_index]
  samples$primary_genotype_status[row_index] <- "selected"
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(samples, output_file, row.names = FALSE, quote = TRUE, na = "NA")

summary_file <- sub("\\.csv$", "_summary.csv", output_file)
missing_summary_file <- sub("\\.csv$", "_missing_by_wgs_column.csv", output_file)

status_summary <- as.data.frame(
  table(samples$primary_genotype_status, useNA = "ifany"),
  stringsAsFactors = FALSE
)
colnames(status_summary) <- c("primary_genotype_status", "n_samples")
write.csv(status_summary, summary_file, row.names = FALSE, quote = TRUE, na = "NA")

missing_rows <- samples[is.na(samples$primary_genotype), , drop = FALSE]
missing_by_column <- data.frame()

for (this_wgs_column in wgs_columns) {
  has_id <- !is.na(missing_rows[[this_wgs_column]]) &
    missing_rows[[this_wgs_column]] != "" &
    missing_rows[[this_wgs_column]] != "NA"

  one_row <- data.frame(
    WGS_metadata_column = this_wgs_column,
    missing_primary_rows_with_this_WGS_id = sum(has_id),
    unique_WGS_ids = length(unique(missing_rows[[this_wgs_column]][has_id])),
    stringsAsFactors = FALSE
  )

  missing_by_column <- rbind(missing_by_column, one_row)
}

write.csv(missing_by_column, missing_summary_file, row.names = FALSE, quote = TRUE, na = "NA")

message("Rows written: ", nrow(samples))
message("Rows with primary genotype: ", sum(!is.na(samples$primary_genotype)))
message("Rows without primary genotype: ", sum(is.na(samples$primary_genotype)))
message("Primary genotype status:")
print(table(samples$primary_genotype_status, useNA = "ifany"))
message("Wrote: ", output_file)
message("Wrote: ", summary_file)
message("Wrote: ", missing_summary_file)
