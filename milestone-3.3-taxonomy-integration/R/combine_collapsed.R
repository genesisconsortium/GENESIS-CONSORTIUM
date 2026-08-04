library(ComplexHeatmap)
library(grid)
library(readxl)
library(stringr)
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(tidyr)

# Folders containing the collapsed CSVs - need to run five_PFC_python_metaneighbor_plots.R on this first
base_folder_A1 <- here("/sbgenomics/project-files/python_Metaneighbor_runs/GEN_A1_pass2.five_PFC_himem.level3.1768273838/")
base_folder_A2 <- here("/sbgenomics/project-files/python_Metaneighbor_runs/GEN_A2_five_PFC_himem.pass3.level3.1768319663/")
base_folder_A3 <- here("/sbgenomics/project-files/python_Metaneighbor_runs/GEN_A3_five_PFC_himem.pass2.level3.1768321612/")

# Paths to the collapsed files
#seems I need to rerun five_PFC_python... and have it delete the old collapsed files before writing new ones! Due to the _1_ issues
file_A1 <- file.path(base_folder_A1, "_1_top_hits.0.8.collapsed.csv")
file_A2 <- file.path(base_folder_A2, "_1_top_hits.0.8.collapsed.csv")
file_A3 <- file.path(base_folder_A3, "_1_top_hits.0.8.collapsed.csv")

# Read tables 
a1 <- read_csv(file_A1, show_col_types = FALSE)
a2 <- read_csv(file_A2, show_col_types = FALSE)
a3 <- read_csv(file_A3, show_col_types = FALSE)

prep_collapsed <- function(df, label) {
  id_cols <- names(df)[grepl("^subclass_hit_", names(df))]
  
  auroc_cols <- names(df)[grepl("^mean_AUROC_", names(df)) & !grepl("across_5", names(df))]
  cells_col  <- "total_cells"
  
  # count col is whatever is left after removing ids + total_cells + AUROCs (+ across_5 if present)
  drop_cols <- c(id_cols, cells_col, auroc_cols, "mean_AUROC_across_5")
  count_candidates <- setdiff(names(df), drop_cols)
  if (length(count_candidates) != 1) {
    stop(
      "Couldn't uniquely identify the per-dataset count column in ", label,
      ". Candidates: ", paste(count_candidates, collapse = ", ")
    )
  }
  count_col <- count_candidates[[1]]
  
  df %>%
    # drop the across-5 column if present (per your request)
    select(-any_of("mean_AUROC_across_5")) %>%
    # rename key numeric columns
    rename(
      !!label := all_of(count_col),
      !!paste0("n_cells_", label) := all_of(cells_col)
    ) %>%
    # suffix AUROC columns so they don't collide in joins
    rename_with(~ paste0(., "_", label), all_of(auroc_cols))
}

a1p <- prep_collapsed(a1, "A1")
a2p <- prep_collapsed(a2, "A2")
a3p <- prep_collapsed(a3, "A3")

id_cols <- names(a1p)[grepl("^subclass_hit_", names(a1p))]

joined <- full_join(a1p, a2p, by = id_cols) %>%
  full_join(a3p, by = id_cols) %>%
  mutate(
    A1 = coalesce(A1, 0),
    A2 = coalesce(A2, 0),
    A3 = coalesce(A3, 0),
    
    n_cells_A1 = coalesce(n_cells_A1, 0),
    n_cells_A2 = coalesce(n_cells_A2, 0),
    n_cells_A3 = coalesce(n_cells_A3, 0),
    
    total   = A1 + A2 + A3,
    n_cells = n_cells_A1 + n_cells_A2 + n_cells_A3
  )

# Average the mean_AUROC_* columns across A1/A2/A3 (excluding any "across_5")
base_auroc <- names(a1)[grepl("^mean_AUROC_", names(a1)) & !grepl("across_5", names(a1))]
for (col in base_auroc) {
  joined[[paste0(col, "_avg")]] <- rowMeans(
    cbind(
      joined[[paste0(col, "_A1")]],
      joined[[paste0(col, "_A2")]],
      joined[[paste0(col, "_A3")]]
    ),
    na.rm = TRUE
  )
}


joined <- joined %>%
  arrange(desc(total)) %>%
  rename_with(~ str_remove_all(., "^subclass_hit_"))


GEN_A1_long_names <- read_csv("/sbgenomics/project-files/GEN_A1/PsychAD_SupplementaryTable2.csv",
                              show_col_types = FALSE) %>%
  select(cxg_GEN_A1 = subclass, long_name) %>%
  distinct() %>%
  group_by(cxg_GEN_A1) %>%
  summarize(PsychAD_long_name = paste(long_name, collapse = ", "), .groups = "drop")

joined %<>% left_join(GEN_A1_long_names)

# Add averaged AUROC (across A1/A2/A3) into the subclass-name columns
joined <- joined %>%
  {
    auroc_avg_cols <- grep("^mean_AUROC_.*_avg$", names(.), value = TRUE)
    dataset_names <- sub("^mean_AUROC_", "", sub("_avg$", "", auroc_avg_cols))
    
    for (ds in dataset_names) {
      subclass_col <- ds
      auroc_col <- paste0("mean_AUROC_", ds, "_avg")
      
      if (subclass_col %in% names(.) && auroc_col %in% names(.)) {
        . <- mutate(
          .,
          !!subclass_col := dplyr::case_when(
            is.na(.data[[subclass_col]]) ~ NA_character_,
            is.na(.data[[auroc_col]]) ~ as.character(.data[[subclass_col]]),
            TRUE ~ paste0(
              .data[[subclass_col]],
              " (",
              format(round(.data[[auroc_col]], 2), nsmall = 2),
              ")"
            )
          )
        )
      }
    }
    .
  }

joined <- joined %>%
  select(-matches("^(AUROC_|mean_AUROC_)"))

joined %<>% select(Ma_2022, SEA_AD, Siletti, freeze2_GEN_A3_ctx, cxg_GEN_A1, PsychAD_long_name, everything())

print(joined, n = 30)
print(joined %>% filter(total == 1), n = 30)

# --- write outputs to a timestamped folder ---
parent_out <- here("/sbgenomics/output-files/")

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(parent_out, paste0("combined_five_PFC.", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(out_dir, "top_hits.0.8.collapsed.combined.csv")
write_csv(joined, out_csv)

out_txt <- file.path(out_dir, "inputs.txt")
writeLines(
  c(
    "Input files:",
    paste0("A1: ", normalizePath(file_A1, winslash = "/", mustWork = FALSE)),
    paste0("A2: ", normalizePath(file_A2, winslash = "/", mustWork = FALSE)),
    paste0("A3: ", normalizePath(file_A3, winslash = "/", mustWork = FALSE)),
    "",
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  con = out_txt
)

message("Wrote combined CSV: ", out_csv)
message("Wrote inputs TXT:   ", out_txt)

#######################################
#######################################
#######################################
#code to make 1v1 heatmaps for rows we want to look at in detail
#######################################
#######################################
#######################################
#get the rows with issues
annotations <- read_xlsx("/sbgenomics/project-files/annotations/GENA1-3_PFC/combined_five_PFC.20260113_210034/GEN-A1-3.top_hits.0.8.collapsed.combined.xlsx")
annotations %<>% filter(Note == "second best?") %>% select(-Note)

annotations %<>%
  # 1) drop numeric columns
  select(where(~ !is.numeric(.))) %>%
  # 2) remove trailing " (0.99)" style AUROC text from all remaining character columns
  mutate(across(
    where(is.character),
    ~ str_trim(str_remove(.x, "\\s*\\([^)]*\\)\\s*$"))
  ))
annotations %<>% rename_with(~ paste0("subclass_hit_", .x))
annotations %<>% select(-subclass_hit_PsychAD_long_name)
annotations


#######################################
#######################################
#set target here
#######################################
#######################################

target_row <- 8
#find the subtypes with those mappings
singe_annotation <- annotations[target_row,]


# If your object is named `singe_annotation`, keep this; otherwise set it to `single_annotation`
# single_annotation <- singe_annotation

load_one_v_one_long <- function(base_folder, dataset_prefix, single_annotation) {
  subtype_mapping <- read_csv(
    file.path(base_folder, "top_hits.0.8.combined.csv"),
    show_col_types = FALSE
  )
  
  target_subtypes <- subtype_mapping %>%
    inner_join(single_annotation) %>%   # uses common columns by default
    pull(subtype) %>%
    unique()
  
  # If none, return an empty tibble with the expected columns
  if (length(target_subtypes) == 0) {
    return(tibble::tibble())
  }
  
  target_cols <- paste0(dataset_prefix, "|", target_subtypes)
  
  one_v_one <- read_csv(
    file.path(base_folder, "aurocs_1v1.csv.gz"),
    show_col_types = FALSE
  ) %>%
    rename(target = "...1")
  
  keep_cols <- intersect(c("target", target_cols), names(one_v_one))
  
  one_v_one %>%
    select(all_of(keep_cols)) %>%
    # drop rows where ALL subtype columns are NA (all columns except target)
    filter(!if_all(-target, is.na)) %>%
    # drop rows where target starts with this dataset's prefix
    filter(!startsWith(target, dataset_prefix)) %>%
    # long format for easy combining
    pivot_longer(
      cols = -target,
      names_to = "query",
      values_to = "AUROC",
      values_drop_na = TRUE
    ) %>%
    mutate(source_dataset = dataset_prefix) %>%
    separate(target, into = c("target_dataset", "target_subtype"), sep = "\\|", remove = FALSE, fill = "right") %>%
    separate(query,  into = c("query_dataset",  "query_subtype"),  sep = "\\|", remove = FALSE, fill = "right")
}

datasets <- list(
  GEN_A1 = list(base_folder = base_folder_A1, prefix = "GEN_A1"),
  GEN_A2 = list(base_folder = base_folder_A2, prefix = "GEN_A2"),
  GEN_A3 = list(base_folder = base_folder_A3, prefix = "GEN_A3")
)

one_v_one_AUROCs_long <- bind_rows(lapply(datasets, \(d)
                                          load_one_v_one_long(d$base_folder, d$prefix, singe_annotation)
))

#load_one_v_one_long(base_folder_A1, "GEN_A1", singe_annotation)
# If you want a combined WIDE matrix back (optional):
one_v_one_AUROCs <- one_v_one_AUROCs_long %>%
  select(target, query, AUROC) %>%
  pivot_wider(names_from = query, values_from = AUROC)

#plot it 
# Build matrix from wide tibble (target column + numeric AUROCs)
one_v_one_AUROCs_matrix <- one_v_one_AUROCs %>%
  rename(target = target) %>%
  {
    mat <- as.matrix(select(., -target))
    rownames(mat) <- pull(., target)
    mat
  }




#########
row_means <- rowMeans(one_v_one_AUROCs_matrix, na.rm = TRUE)
col_means <- colMeans(one_v_one_AUROCs_matrix, na.rm = TRUE)

one_v_one_AUROCs_matrix <- one_v_one_AUROCs_matrix[
  order(row_means, decreasing = TRUE, na.last = TRUE),
  order(col_means, decreasing = TRUE, na.last = TRUE),
  drop = FALSE
]

# Font sizes that still show ALL labels
n_rows <- nrow(one_v_one_AUROCs_matrix)
n_cols <- ncol(one_v_one_AUROCs_matrix)
row_fs <- max(12, min(10, round(140 / max(1, n_rows))))
col_fs <- max(12, min(10, round(140 / max(1, n_cols))))

mat_plot <- one_v_one_AUROCs_matrix
mat_plot[is.na(mat_plot)] <- 0.5  # neutral AUROC baseline

ht <- Heatmap(
  mat_plot,
  name = "AUROC",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = row_fs),
  column_names_gp = gpar(fontsize = col_fs),
  row_names_max_width = unit(20, "cm"),
  column_names_max_height = unit(14, "cm"),
  column_names_rot = 0,
  column_title = singe_annotation[1,] %>% paste0(collapse = ", ")
)
draw(ht, heatmap_legend_side = "right")



#######################################
#######################################
#######################################
#write out mappings
#######################################
#######################################
#######################################

#get the rows with issues
#annotation_filename <- "/sbgenomics/project-files/annotations/GENA1-3_PFC/combined_five_PFC.20260113_210034/GEN-A1-3.top_hits.0.8.collapsed.combined.xlsx"
annotation_filename <- "/sbgenomics/project-files/annotations/GENA1-3_PFC/combined_five_PFC.20260113_210034/GEN-A1-3.top_hits.0.8.collapsed.combined.micro_PVM_merge.xlsx"
annotations <- read_xlsx(annotation_filename)
annotations %<>% filter(!grepl("^x", Note))
annotations %>% pull(Note) %>% unique()

annotations %<>%
  # 1) drop numeric columns
  select(where(~ !is.numeric(.))) %>%
  # 2) remove trailing " (0.99)" style AUROC text from all remaining character columns
  mutate(across(
    where(is.character),
    ~ str_trim(str_remove(.x, "\\s*\\([^)]*\\)\\s*$"))
  ))
annotations %<>% rename_with(~ paste0("subclass_hit_", .x))
annotations %<>% select(-subclass_hit_PsychAD_long_name)
annotations %<>% rename(subclass_annotation_v1 = subclass_hit_Note)
annotations %<>%
  mutate(across(where(is.character), ~ str_replace_all(.x, fixed("NA"), "NA_STRING")))

map_one_folder <- function(base_folder_target) {
  subtype_mapping <- read_csv(
    file.path(base_folder_target, "top_hits.0.8.combined.csv"),
    show_col_types = FALSE
  ) %>%
    select(where(is.character)) %>%
    mutate(across(everything(), ~ replace_na(.x, "NA_STRING")))
  
  joined <- inner_join(subtype_mapping, annotations) %>%
    select(subtype, subclass_annotation_v1) %>%
    mutate(
      base_folder = base_folder_target,
      annotation_file = annotation_filename
    )
  
  # robust parsing for both:
  # - GEN_A1_pass2.five_...
  # - GEN_A2_five_....pass3.level3....
  folder_name <- basename(sub("/$", "", base_folder_target))
  
  joined <- joined %>%
    mutate(
      dataset = str_extract(folder_name, "GEN_[A-Z0-9]+"),
      pass_num = str_extract(folder_name, "(?<=_pass)\\d+|(?<=\\.pass)\\d+"),
      dataset_pass = if_else(
        !is.na(dataset) & !is.na(pass_num),
        paste0(dataset, "_pass", pass_num),
        NA_character_
      )
    ) %>%
    select(subtype, subclass_annotation_v1, dataset, dataset_pass, everything()) %>%
    arrange(subclass_annotation_v1)
  
  joined
}
all_mapped <- bind_rows(map_one_folder(base_folder_A1),
          map_one_folder(base_folder_A2),
          map_one_folder(base_folder_A3))
all_mapped %>% group_by(subclass_annotation_v1) %>% count() %>% as.data.frame()

parent_out <- here("/sbgenomics/output-files/")
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
map_out_dir <- file.path(parent_out, paste0("mappings_for_combined_five_PFC.", ts))
dir.create(map_out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(map_out_dir, "subtype_mappings.csv")
write_csv(all_mapped, out_csv)
out_csv
