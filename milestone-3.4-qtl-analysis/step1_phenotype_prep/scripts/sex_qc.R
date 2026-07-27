.libPaths(c("/sc/arion/projects/psychAD/aging/kiran/RLib_4_5", .libPaths()))

library(data.table)
library(dreamlet)
library(SummarizedExperiment)
library(openxlsx)

single_cell_file <- "/sc/arion/projects/CommonMind/genesis/metadata/outputs/samples_single_cell.csv"
base_dir <- "/sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/outputs"
depth_file <- "/sc/arion/projects/CommonMind/genesis/Genotype_files/sex_chr_depth_qc/WGS_mosdepth_chrX_chrY_autosome_ratios.tsv"

out_dir <- "/sc/arion/projects/CommonMind/genesis/pseudobulk_dreamlet_rc2/sex_marker_qc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_xlsx <- file.path(out_dir, "sex_QC_failed_samples_for_review.xlsx")
out_tsv <- file.path(
  out_dir,
  "samples_single_cell_with_RNA_median_and_WGS_sex_QC.tsv"
)
out_csv <- file.path(out_dir, "sex_qc_samples.csv")
out_outliers_csv <- file.path(out_dir, "sex_qc_outliers.csv")

clean_gt <- function(x) {
  x <- as.character(x)
  x <- sub("_vcpa1\\.1$", "", x)
  x[x == "" | x == "NA"] <- NA_character_
  x
}

clean_sex <- function(x) {
  x <- trimws(tolower(as.character(x)))
  fifelse(
    x == "male",
    "Male",
    fifelse(x == "female", "Female", NA_character_)
  )
}

clean_genesis_sample <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | x == "NA"] <- NA_character_

  # Harmonize CMC IDs between samples_single_cell.csv and RDS colData:
  # CMC_MSSM_190_PFC_A5 -> CMC-MSSM-190_PFC_A5
  x <- sub("^CMC_MSSM_", "CMC-MSSM-", x)
  x
}

first_nonmissing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[1]
}

sc <- fread(single_cell_file)
sc[, primary_clean := clean_gt(primary_gt)]
sc[, sex_clean := clean_sex(sex)]
sc[, has_primary_gt := !is.na(primary_clean) & primary_clean != ""]

if (!"sample_id" %in% names(sc) && "individualID" %in% names(sc)) {
  sc[, sample_id := individualID]
}

if (!"snRNAseq_QC_outlier" %in% names(sc)) {
  sc[, snRNAseq_QC_outlier := FALSE]
}

sc[, rna_join_id := clean_genesis_sample(genesis_sample)]

# Build a project- and region-specific lookup for RDS files whose expression
# columns contain individual IDs rather than full genesis_sample IDs.
sc[, individual_join_id := fifelse(
  !is.na(individualID) &
    trimws(as.character(individualID)) != "",
  trimws(as.character(individualID)),
  trimws(as.character(sample_id))
)]
sc[, brain_region_join := trimws(as.character(brain_region))]
sc[, genesis_id_join := trimws(as.character(genesis_id))]

sc_rna_lookup <- unique(
  sc[
    !is.na(individual_join_id) &
      individual_join_id != "" &
      !is.na(brain_region_join) &
      brain_region_join != "" &
      !is.na(genesis_id_join) &
      genesis_id_join != "" &
      !is.na(rna_join_id),
    .(
      individual_join_id,
      brain_region_join,
      genesis_id_join,
      rna_join_id_from_sc = rna_join_id
    )
  ]
)

ambiguous_sc_lookup <- sc_rna_lookup[
  ,
  .(
    n_genesis_samples = uniqueN(rna_join_id_from_sc),
    genesis_samples = paste(
      sort(unique(rna_join_id_from_sc)),
      collapse = ";"
    )
  ),
  by = .(
    individual_join_id,
    brain_region_join,
    genesis_id_join
  )
][n_genesis_samples > 1]

if (nrow(ambiguous_sc_lookup) > 0) {
  print(ambiguous_sc_lookup)
  stop("Ambiguous RNA-to-single-cell sample mappings detected.")
}

# Explicitly retain S17783 / NDAR_INVZR997PX1 as not being an
# snRNA-seq QC outlier.
sc[sample_id == "S17783", snRNAseq_QC_outlier := FALSE]

rds_files <- list.files(
  base_dir,
  pattern = "_class_final_processAssays\\.RDS$",
  recursive = TRUE,
  full.names = TRUE
)
rds_files <- rds_files[grepl("/processassays/", rds_files)]

extract_all_classes_xist_uty <- function(rds_file) {
  message("Reading: ", rds_file)

  obj <- readRDS(rds_file)
  project_folder <- basename(dirname(dirname(rds_file)))
  rds_name <- basename(rds_file)

  cd <- as.data.table(
    colData(obj),
    keep.rownames = "pseudobulk_sample_id"
  )

  cd_small <- cd[
    ,
    .(
      pseudobulk_sample_id,
      primary_genotype,
      coldata_sex = sex,
      coldata_cohort = cohort,
      coldata_subcohort = subcohort,
      coldata_ancestry = ancestry,
      genesis_sample,
      genesis_id,
      brain_region
    )
  ]

  class_dt <- lapply(names(obj), function(cell_class) {
    elist <- obj[[cell_class]]
    if (is.null(elist$E)) return(NULL)

    expr <- elist$E
    if (!all(c("XIST", "UTY") %in% rownames(expr))) return(NULL)

    dt <- data.table(
      pseudobulk_sample_id = colnames(expr),
      cell_class = cell_class,
      XIST = as.numeric(expr["XIST", ]),
      UTY = as.numeric(expr["UTY", ]),
      project_folder = project_folder,
      rds_file = rds_name
    )

    merge(
      dt,
      cd_small,
      by = "pseudobulk_sample_id",
      all.x = TRUE
    )
  })

  rbindlist(class_dt, fill = TRUE)
}

xist_uty_all_classes <- rbindlist(
  lapply(rds_files, extract_all_classes_xist_uty),
  fill = TRUE
)

saveRDS(
  xist_uty_all_classes,
  file.path(out_dir, "all_cell_classes_XIST_UTY_with_metadata.rds")
)

# First try the normalized genesis_sample supplied by the RDS colData.
xist_uty_all_classes[
  ,
  rna_join_id_exact := clean_genesis_sample(genesis_sample)
]

# Also construct the project- and region-specific individual lookup key. This
# handles RDS expression columns such as H19.33.004, which need to map to a
# full sample ID such as H19.33.004_MEC_A4.
xist_uty_all_classes[
  ,
  individual_join_id := trimws(as.character(pseudobulk_sample_id))
]
xist_uty_all_classes[
  ,
  brain_region_join := trimws(as.character(brain_region))
]
xist_uty_all_classes[
  ,
  genesis_id_join := trimws(as.character(genesis_id))
]

xist_uty_all_classes <- merge(
  xist_uty_all_classes,
  sc_rna_lookup,
  by = c(
    "individual_join_id",
    "brain_region_join",
    "genesis_id_join"
  ),
  all.x = TRUE
)

valid_sc_rna_ids <- unique(sc[!is.na(rna_join_id), rna_join_id])

# Use an exact genesis_sample match when it exists in the single-cell
# metadata. Otherwise use individualID + brain_region + genesis_id.
xist_uty_all_classes[, rna_join_id := NA_character_]
xist_uty_all_classes[
  !is.na(rna_join_id_exact) &
    rna_join_id_exact %chin% valid_sc_rna_ids,
  rna_join_id := rna_join_id_exact
]
xist_uty_all_classes[
  is.na(rna_join_id) & !is.na(rna_join_id_from_sc),
  rna_join_id := rna_join_id_from_sc
]

# Report expression records that still cannot be assigned to a sample.
unmatched_expression <- unique(
  xist_uty_all_classes[
    is.na(rna_join_id),
    .(
      pseudobulk_sample_id,
      genesis_sample,
      genesis_id,
      brain_region,
      primary_genotype,
      project_folder
    )
  ]
)

if (nrow(unmatched_expression) > 0) {
  message(
    "Unmatched expression records after RNA sample mapping: ",
    nrow(unmatched_expression)
  )
  print(unmatched_expression)
  fwrite(
    unmatched_expression,
    file.path(out_dir, "unmatched_RNA_expression_records.tsv"),
    sep = "\t"
  )
}

# Stop if normalization maps multiple source samples or genotypes to the same
# RNA join ID. This prevents accidental expression aggregation across samples.
rna_join_check <- xist_uty_all_classes[
  !is.na(rna_join_id),
  .(
    n_genesis_samples = uniqueN(
      as.character(genesis_sample),
      na.rm = TRUE
    ),
    n_primary_genotypes = uniqueN(
      as.character(primary_genotype),
      na.rm = TRUE
    )
  ),
  by = rna_join_id
]

problematic_rna_ids <- rna_join_check[
  n_genesis_samples > 1 | n_primary_genotypes > 1
]

if (nrow(problematic_rna_ids) > 0) {
  print(problematic_rna_ids)
  stop("RNA join-ID normalization produced ambiguous sample mappings.")
}

rna_median <- xist_uty_all_classes[
  !is.na(rna_join_id),
  .(
    median_XIST_all_classes = median(XIST, na.rm = TRUE),
    median_UTY_all_classes = median(UTY, na.rm = TRUE),
    n_cell_classes_with_XIST_UTY = uniqueN(cell_class),
    cell_classes_with_XIST_UTY = paste(
      sort(unique(cell_class)),
      collapse = ";"
    ),
    pseudobulk_sample_id_from_expr = first_nonmissing(
      pseudobulk_sample_id
    ),
    genesis_sample_from_expr = first_nonmissing(genesis_sample),
    primary_genotype_from_pseudobulk = first_nonmissing(
      primary_genotype
    ),
    project_folder_from_expr = first_nonmissing(project_folder)
  ),
  by = rna_join_id
]

sc_qc <- merge(
  sc,
  rna_median,
  by = "rna_join_id",
  all.x = TRUE
)

sc_qc[, RNA_expr_sex := fifelse(
  median_UTY_all_classes > 6 & median_XIST_all_classes < 7,
  "Male",
  fifelse(
    median_XIST_all_classes > 7 & median_UTY_all_classes < 6,
    "Female",
    "Review"
  )
)]

sc_qc[
  is.na(median_XIST_all_classes) | is.na(median_UTY_all_classes),
  RNA_expr_sex := NA_character_
]

sc_qc[, RNA_qc := "pass"]
sc_qc[is.na(RNA_expr_sex), RNA_qc := "not_tested"]
sc_qc[
  RNA_expr_sex == "Review",
  RNA_qc := "possible_sex_chromosome_abnormality"
]

sc_qc[
  RNA_expr_sex %in% c("Male", "Female") &
    !is.na(sex_clean) &
    RNA_expr_sex != sex_clean,
  RNA_qc := "potential_swap"
]

depth <- fread(depth_file)
depth[, primary_clean := clean_gt(sample_id)]
depth[, chrX_autosome_ratio := as.numeric(chrX_autosome_ratio)]
depth[, chrY_autosome_ratio := as.numeric(chrY_autosome_ratio)]

depth_small <- depth[
  ,
  .(
    primary_clean,
    mean_autosome_depth,
    mean_chrX_depth,
    mean_chrY_depth,
    chrX_autosome_ratio,
    chrY_autosome_ratio
  )
]

sc_qc <- merge(sc_qc, depth_small, by = "primary_clean", all.x = TRUE)

sc_qc[, WGS_depth_sex := NA_character_]

sc_qc[
  chrX_autosome_ratio >= 0.75 &
    chrX_autosome_ratio <= 1.10 &
    chrY_autosome_ratio >= 0.00 &
    chrY_autosome_ratio <= 0.25,
  WGS_depth_sex := "Female"
]

sc_qc[
  chrX_autosome_ratio >= 0.40 &
    chrX_autosome_ratio <= 0.65 &
    chrY_autosome_ratio >= 0.50 &
    chrY_autosome_ratio <= 1.20,
  WGS_depth_sex := "Male"
]

sc_qc[, WGS_qc := "pass"]
sc_qc[has_primary_gt == FALSE, WGS_qc := "no_primary_genotype"]

sc_qc[
  has_primary_gt == TRUE & is.na(chrX_autosome_ratio),
  WGS_qc := "not_tested"
]

sc_qc[
  has_primary_gt == TRUE &
    !is.na(chrX_autosome_ratio) &
    !is.na(chrY_autosome_ratio) &
    is.na(WGS_depth_sex),
  WGS_qc := "possible_sex_chromosome_abnormality"
]

sc_qc[
  !is.na(WGS_depth_sex) &
    !is.na(sex_clean) &
    WGS_depth_sex != sex_clean,
  WGS_qc := "potential_swap"
]

sc_qc[
  ,
  RNA_failed := RNA_qc %in%
    c("possible_sex_chromosome_abnormality", "potential_swap")
]

sc_qc[
  ,
  WGS_failed := WGS_qc %in%
    c("possible_sex_chromosome_abnormality", "potential_swap")
]

sc_qc[, final_sex_qc := fifelse(
  has_primary_gt == FALSE,
  "no_primary_genotype",
  fifelse(
    snRNAseq_QC_outlier == TRUE,
    "snRNAseq_QC_outlier",
    fifelse(
      RNA_failed == TRUE & WGS_failed == TRUE,
      "RNA_and_WGS_sex_QC_outlier",
      fifelse(
        RNA_failed == TRUE,
        "RNA_sex_QC_outlier",
        fifelse(
          WGS_failed == TRUE,
          "WGS_sex_QC_outlier",
          "pass"
        )
      )
    )
  )
)]

summary_final <- sc_qc[, .N, by = final_sex_qc][order(-N)]

summary_by_project <- sc_qc[
  final_sex_qc != "pass",
  .N,
  by = .(
    GENESIS_study,
    genesis_id,
    cohort,
    subcohort,
    final_sex_qc
  )
][order(GENESIS_study, genesis_id, final_sex_qc)]

if (!"sample_id" %in% names(sc_qc) &&
    "individualID" %in% names(sc_qc)) {
  sc_qc[, sample_id := individualID]
}

review_cols <- c(
  "genesis_sample",
  "genesis_sample_from_expr",
  "rna_join_id",
  "sample_id",
  "primary_gt",
  "primary_clean",
  "sex",
  "sex_clean",
  "cohort",
  "subcohort",
  "ancestry",
  "GENESIS_study",
  "genesis_id",
  "brain_region",
  "snRNAseq_QC_outlier",
  "median_XIST_all_classes",
  "median_UTY_all_classes",
  "RNA_expr_sex",
  "RNA_qc",
  "WGS_depth_sex",
  "WGS_qc",
  "final_sex_qc",
  "mean_autosome_depth",
  "mean_chrX_depth",
  "mean_chrY_depth",
  "chrX_autosome_ratio",
  "chrY_autosome_ratio",
  "n_cell_classes_with_XIST_UTY",
  "cell_classes_with_XIST_UTY"
)

review_cols <- intersect(review_cols, names(sc_qc))

failed_for_review <- sc_qc[
  final_sex_qc != "pass",
  ..review_cols
][order(final_sex_qc, GENESIS_study, genesis_sample)]

if ("sex" %in% names(failed_for_review)) {
  setnames(failed_for_review, "sex", "sex_metadata")
}

readme <- data.table(
  field = c(
    "RNA_qc",
    "WGS_qc",
    "final_sex_qc",
    "RNA join ID",
    "RNA male rule",
    "RNA female rule",
    "WGS female rule",
    "WGS male rule",
    "S17783"
  ),
  description = c(
    paste(
      "RNA sex QC from median XIST/UTY across all available",
      "pseudobulk cell classes."
    ),
    paste(
      "WGS sex QC from chrX/autosome and chrY/autosome",
      "mosdepth depth ratios."
    ),
    "Final combined sex-QC decision used for exclusion/review.",
    paste(
      "RNA expression is joined using normalized genesis_sample;",
      "CMC_MSSM_ is harmonized to CMC-MSSM-."
    ),
    paste(
      "Male if median_UTY_all_classes > 6 and",
      "median_XIST_all_classes < 7."
    ),
    paste(
      "Female if median_XIST_all_classes > 7 and",
      "median_UTY_all_classes < 6."
    ),
    paste(
      "Female WGS depth if chrX/autosome is 0.75-1.10 and",
      "chrY/autosome is 0.00-0.25."
    ),
    paste(
      "Male WGS depth if chrX/autosome is 0.40-0.65 and",
      "chrY/autosome is 0.50-1.20."
    ),
    paste(
      "S17783 / NDAR_INVZR997PX1 is explicitly retained with",
      "snRNAseq_QC_outlier = FALSE."
    )
  )
)

# Create a CSV with names compatible with the earlier sex_qc_samples.csv
# output while retaining all columns from the revised sample-level analysis.
sex_qc_samples <- copy(sc_qc)

sex_qc_samples[, primary_gt_clean := primary_clean]
sex_qc_samples[
  ,
  has_XIST_UTY_any_class :=
    !is.na(median_XIST_all_classes) &
    !is.na(median_UTY_all_classes)
]
sex_qc_samples[
  ,
  RNA_expr_sex_median := RNA_expr_sex
]
sex_qc_samples[
  ,
  RNA_median_sex_qc := RNA_qc
]
sex_qc_samples[
  ,
  metadata_sex_clean := sex_clean
]
sex_qc_samples[
  ,
  WGS_sex_qc := WGS_qc
]
sex_qc_samples[
  ,
  project_for_plot := project_folder_from_expr
]

# Match the file_region convention used by the earlier sex_qc_samples.csv:
# A3 and A4 are region-specific; other GENESIS projects use only genesis_id.
sex_qc_samples[
  ,
  file_region := fifelse(
    is.na(genesis_id) | trimws(as.character(genesis_id)) == "",
    NA_character_,
    fifelse(
      genesis_id %in% c("A3", "A4"),
      paste0(
        "GEN_",
        as.character(genesis_id),
        "_",
        as.character(brain_region)
      ),
      paste0("GEN_", as.character(genesis_id))
    )
  )
]

# Explicit RNA sex-QC grouping for plots and review.
sex_qc_samples[
  ,
  RNA_sex_QC_group := fifelse(
    RNA_qc == "potential_swap",
    "Potential swap",
    fifelse(
      RNA_qc == "possible_sex_chromosome_abnormality",
      "Possible abnormality",
      fifelse(
        RNA_qc == "not_tested",
        "Not tested",
        "Pass"
      )
    )
  )
]

# Broader mutually exclusive group based on the final combined QC call.
sex_qc_samples[
  ,
  final_sex_QC_group := fifelse(
    final_sex_qc == "RNA_sex_QC_outlier",
    "RNA outlier",
    fifelse(
      final_sex_qc == "WGS_sex_QC_outlier",
      "WGS outlier",
      fifelse(
        final_sex_qc == "RNA_and_WGS_sex_QC_outlier",
        "RNA + WGS outlier",
        fifelse(
          final_sex_qc == "snRNAseq_QC_outlier",
          "snRNA-seq outlier",
          fifelse(
            final_sex_qc == "no_primary_genotype",
            "No primary genotype",
            "Pass"
          )
        )
      )
    )
  )
]

outlier_categories <- c(
  "RNA_sex_QC_outlier",
  "WGS_sex_QC_outlier",
  "RNA_and_WGS_sex_QC_outlier",
  "snRNAseq_QC_outlier"
)

sex_qc_outliers <- sex_qc_samples[
  final_sex_qc %in% outlier_categories
]

fwrite(sc_qc, out_tsv, sep = "\t")
fwrite(sex_qc_samples, out_csv)
fwrite(sex_qc_outliers, out_outliers_csv)

wb <- createWorkbook()

addWorksheet(wb, "failed_samples_for_review")
writeData(wb, "failed_samples_for_review", failed_for_review)
freezePane(wb, "failed_samples_for_review", firstRow = TRUE)
addFilter(
  wb,
  "failed_samples_for_review",
  rows = 1,
  cols = seq_len(ncol(failed_for_review))
)

addWorksheet(wb, "summary_final")
writeData(wb, "summary_final", summary_final)

addWorksheet(wb, "summary_by_project")
writeData(wb, "summary_by_project", summary_by_project)

addWorksheet(wb, "README")
writeData(wb, "README", readme)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

print(summary_final)
message("Saved full TSV: ", out_tsv)
message("Saved full sex-QC CSV: ", out_csv)
message("Saved outlier-only CSV: ", out_outliers_csv)
message("Saved review Excel: ", out_xlsx)

