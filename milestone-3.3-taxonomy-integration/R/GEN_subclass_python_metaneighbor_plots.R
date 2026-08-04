library(circlize)
library(ComplexHeatmap)
library(grid)
library(igraph)
library(tidygraph)
library(ggraph)   
library(ggplot2)  
library(scales) 
library(igraph)
library(tidyr)
library(reshape2)
library(readr)
library(ComplexHeatmap)
library(dplyr)
library(MetaNeighbor)
library(magrittr)
library(here)

# This R script processes pyMN/MetaNeighbor outputs to evaluate the
# replicability of cell-type mappings across studies. It identifies the
# primary target dataset, applies concordance filters based on all-vs-all
# and best-vs-next AUROC results, summarizes passing and failing subtypes,
# and exports mapping tables and summary statistics.
#
# It also contains several heatmap, distribution, and network
# visualizations. Many versions of these visualizations are retained
# because their design and ordering methods have changed over time.
#
# The script is currently set up for interactive analysis and inspection,
# but could be refactored into modular functions or a more automated
# workflow.

####
base_folder <- NULL
#this is set to the current pyMN output script location
base_folder <- here("/sbgenomics/project-files/python_Metaneighbor_runs/GEN_A6_vrs_GEN_A1-2.cortex.subclass.1782242432/")

output_folder <- gsub("project-files", "output-files", base_folder)
output_folder <- file.path(output_folder, format(Sys.time(), "R_results_%Y-%m-%d_%H-%M/"))
dir.create(paste0(output_folder), recursive = TRUE)

obs <- read_csv(paste0(base_folder, "merged.obs.csv.gz"))



#select the one with the most cell types (subtypes)
genesis_study_name <- obs %>% group_by(study_id) %>% select(cell_type)%>% distinct() %>% count()%>% arrange(-n) %>% pull(study_id) %>% first()
genesis_study_name

#old appraoch
#(genesis_study_name <- obs %>% filter(grepl("^GEN_",study_id)) %>% pull(study_id) %>% unique() )
obs %>% select(cell_type, study_id) %>% distinct() %>% group_by(study_id) %>% count()
obs %>% group_by(study_id) %>% count()
#obs %>% group_by(cell_type) %>% count()

studies_keep <- obs %>% pull(study_id) %>% unique()
other_studies <- obs %>% filter(study_id != genesis_study_name) %>% pull(study_id) %>% unique()
## ---- study filter ----
#studies_keep <- c(genesis_study_name, "GEN_A4_PFC", "GEN_A5", "GEN_A13", "GEN_A16") 
other_studies <- setdiff(other_studies, genesis_study_name)
obs %<>% filter(study_id %in% studies_keep)
target_types <- obs %>% filter(study_id == genesis_study_name) %>% pull(cell_type) %>% unique()
other_subclasses <- obs %>% filter(study_id != genesis_study_name) %>% mutate(combined = paste0(study_id, "|", cell_type)) %>% pull(combined) %>% unique()
obs %>% filter(study_id == genesis_study_name) %>% nrow()

occurances_in_refs <- obs %>% filter(study_id != genesis_study_name) %>% select(cell_type,study_id) %>% distinct() %>% group_by(cell_type) %>% count()
occurances_in_refs %<>% arrange(n)

nrow(obs)
var <- read_csv(paste0(base_folder, "merged.var.csv.gz"))
nrow(var)

fmt_n <- function(x) { format(x, big.mark = ",", scientific = FALSE, trim = TRUE) } #commas to numbers for readability

top_hits <- read_csv(paste0(base_folder, "top_hits.0.8.csv"))
top_hits %<>% select(-`...1`)
top_hits %<>% mutate(study_1 = getStudyId(`Study_ID|Celltype_1`))
top_hits %<>% mutate(study_2 = getStudyId(`Study_ID|Celltype_2`))

#study filter
top_hits %<>% filter(study_1 %in% studies_keep, study_2 %in% studies_keep)

top_hits %>% filter( Match_type == "Reciprocal_top_hit")
top_hits %>% filter( Match_type != "Reciprocal_top_hit")
top_hits %>% filter(study_2 == genesis_study_name | study_1 == genesis_study_name) %>% filter( Match_type == "Reciprocal_top_hit")

top_hits %<>% 
  # mark rows that need to be swapped
  mutate(swap = study_2 == genesis_study_name & study_1 != genesis_study_name) %>% 
  mutate(
    cell1_tmp  = `Study_ID|Celltype_1`,
    study1_tmp = study_1,
    `Study_ID|Celltype_1` = if_else(swap, `Study_ID|Celltype_2`, `Study_ID|Celltype_1`),
    study_1               = if_else(swap, study_2,               study_1),
    `Study_ID|Celltype_2` = if_else(swap, cell1_tmp,             `Study_ID|Celltype_2`),
    study_2               = if_else(swap, study1_tmp,            study_2)
  ) %>% 
  select(-swap, -cell1_tmp, -study1_tmp)  
top_hits %<>% rename(Celltype_1 = `Study_ID|Celltype_1`, Celltype_2 = `Study_ID|Celltype_2`)



top_hits_GEN <- top_hits %>% filter(study_1 == genesis_study_name) %>% select(subtype = Celltype_1, subclass_hit = Celltype_2, study_2, Mean_AUROC)
top_hits_GEN %>% group_by(subtype, study_2) %>% count() %>% filter(n!=1)
#due to reference directions
top_hits_GEN %<>% 
  # keep only the row with the highest AUROC for each subtype–study_2 pair
  group_by(subtype, study_2) %>% 
  slice_max(Mean_AUROC, n = 1, with_ties = FALSE) %>% 
  ungroup() 

top_hits_GEN %<>% mutate(subclass_hit = getCellType(subclass_hit))

for_filter <- top_hits_GEN
for_filter %<>% inner_join(occurances_in_refs %>% rename(subclass_hit = cell_type, subclass_occurances = n))
for_filter %<>% group_by(subtype) %>%
  arrange(desc(Mean_AUROC), .by_group = TRUE) %>%
  group_modify(~{
    dat <- .x
    
    # best row for this subtype
    best_row <- dat %>% slice(1)
    best_hit <- best_row$subclass_hit
    n_hits <- best_row$subclass_occurances
    
    # protect against asking for more rows than exist
    top_rows <- dat %>% slice_head(n = min(n_hits, nrow(dat)))
    passes <- all(top_rows$subclass_hit == best_hit)
    tibble(
      best_subclass_hit = best_hit,
      best_Mean_AUROC = best_row$Mean_AUROC,
      subclass_occurances = n_hits,
      passes_all_vs_all_filter = passes
    )
  }) %>%
  ungroup()
for_filter %>% arrange(subclass_occurances)
for_filter %>% filter(!passes_all_vs_all_filter) %>% count()
for_filter %>% group_by(passes_all_vs_all_filter) %>% count()
top_hits_GEN %>% filter(subtype == "GEN_A13|2_5_1" )
top_hits_GEN %>% filter(subtype == "GEN_A13|10_3_4" )
top_hits_GEN %>% filter(subtype == "GEN_A13|2_3_1" )

#run same filter on best_vs_next, then combine and filter downstream? 
one_v_one_AUROCs <- read_csv(paste0(base_folder, "aurocs_1v1.csv.gz"))
one_v_one_AUROCs %<>% rename(target = "...1")
one_v_one_AUROCs_long <- one_v_one_AUROCs %>%
  pivot_longer(
    cols = -target,
    names_to = "reference",
    values_to = "AUROC"
  ) %>%
  filter(!is.na(AUROC))

one_v_one_AUROCs_long %<>% filter(getStudyId(reference) == genesis_study_name) %>% 
  rename(subclass_hit = target, subtype = reference) %>% 
  mutate(study_2 = getStudyId(subclass_hit), subclass_hit = getCellType(subclass_hit))
#one_v_one_AUROCs_long %>% filter(subtype == "GEN_A13|4_7_1") 

one_v_one_AUROCs_long %<>% inner_join(occurances_in_refs %>% rename(subclass_hit = cell_type, subclass_occurances = n))
one_v_one_AUROCs_long %<>% filter(AUROC > 0.6)
#code dupe from above
one_v_one_AUROCs_long_summary <- one_v_one_AUROCs_long %>% group_by(subtype) %>%
  arrange(desc(AUROC), .by_group = TRUE) %>%
  group_modify(~{
    dat <- .x
    
    # best row for this subtype
    best_row <- dat %>% slice(1)
    best_hit <- best_row$subclass_hit
    n_hits <- best_row$subclass_occurances
    
    # protect against asking for more rows than exist
    top_rows <- dat %>% slice_head(n = min(n_hits, nrow(dat)))
    passes <- all(top_rows$subclass_hit == best_hit)
    tibble(
      best_subclass_hit = best_hit,
      best_AUROC = best_row$AUROC,
      subclass_occurances = n_hits,
      passes_best_vs_next_filter = passes
    )
  }) %>%
  ungroup()

one_v_one_AUROCs_long_summary %>% filter(subtype == "GEN_A13|4_7_1") 
top_hits_GEN %>% filter(subtype == "GEN_A13|6_6_5") %>% arrange(-Mean_AUROC)
one_v_one_AUROCs_long %>% filter(subtype == "GEN_A13|6_6_8") %>% arrange(-AUROC)
top_hits_GEN %>% filter(subtype == "GEN_A13|6_6_8") %>% arrange(-Mean_AUROC)
#require global hit is same as local
one_v_one_AUROCs_long_summary %<>% inner_join(for_filter %>% select(subtype, best_subclass_hit) %>% unique())
for_filter %<>% left_join(one_v_one_AUROCs_long_summary %>% select(subtype, passes_best_vs_next_filter))
for_filter %<>% mutate(passes_both_filters = passes_all_vs_all_filter & passes_best_vs_next_filter)

for_filter %>% filter(!passes_both_filters) %>% count() #this is NOT the full count of rejects
for_filter %>% filter(passes_both_filters) %>% count() 

subtypes_passing_filter_with_studyID <- for_filter %>% filter(passes_both_filters) %>% pull(subtype)
subtypes_passing_filter <- subtypes_passing_filter_with_studyID %>% getCellType()
subtypes_failing_filter <- setdiff(target_types, subtypes_passing_filter)

write_subtype_file <- function(x, file) {
  lines <- c("subtype", as.character(x))
  writeLines(lines, con = file)
}
write_subtype_file(subtypes_passing_filter, paste0(output_folder, "passing_subtypes.txt"))
write_subtype_file(subtypes_failing_filter, paste0(output_folder, "failing_subtypes.txt"))

for_filter %>% filter(passes_both_filters) %>% select(subtype, best_subclass_hit) %>% mutate(subtype = getCellType(subtype)) %>% 
  write_csv(paste0(output_folder, "assigned_subclasses_for_passing_subtypes.csv"))

missing_subclasses <- setdiff( other_subclasses %>% getCellType() %>% unique(), for_filter %>% filter(passes_both_filters) %>% select(subtype, best_subclass_hit) %>% mutate(subtype = getCellType(subtype)) %>% pull(best_subclass_hit) %>% unique())

subtypes_passing_filter

#get proportion passing filter
length(subtypes_passing_filter)
proportion_subtypes_passing_filter <- length(subtypes_passing_filter)/length(target_types)
proportion_cells_passing_filter <- (obs %>% filter(cell_type %in% subtypes_passing_filter) %>% nrow()) / (obs %>% filter(cell_type %in% target_types) %>% nrow())
top_hits_wide <- top_hits_GEN %>% 
  pivot_wider(
    id_cols    = subtype,
    names_from  = study_2, 
    values_from = c(Mean_AUROC, subclass_hit),
    names_glue  = "{.value}_{study_2}",   # e.g. Mean_AUROC_human, subclass_hit_human
    values_fill = list(Mean_AUROC = NA_real_, subclass_hit = NA_character_)
  )

top_hits_wide %<>% 
  rowwise() %>%                                # work row-by-row
  mutate(
    Avg_AUROC   = mean(c_across(starts_with("Mean_AUROC_")), na.rm = TRUE),
  ) %>% 
  ungroup()

top_hits_wide %>% arrange(Avg_AUROC) %>% head(20)
#calc average AUROC
top_hits_wide %>% pull(Avg_AUROC) %>% median()

top_hits_wide <- top_hits_wide %>% 
  ## 1. Drop the "Mean_" prefix from AUROC columns
  rename_with(~ sub("^Mean_", "", .x), starts_with("Mean_"))

top_hits_wide %>% arrange(Avg_AUROC)

#write out with parantheses
top_hits_wide
#missing hits
top_hits_wide %<>% mutate(subtype = getCellType(subtype))
subtypes_without_tophits <- setdiff(target_types, top_hits_wide %>% pull(subtype))
top_hits_wide %<>% arrange(Avg_AUROC)
#colnames(top_hits_wide) <- colnames(top_hits_wide) %>% gsub("subclass_hit_", "", .) #probably breaks the combined below
top_hits_wide %>% mutate(Manual_agreement = "") %>% write_csv(paste0(output_folder, "top_hits.0.8.combined.csv"))
print(top_hits_wide,n=30)



top_hits_combined <- top_hits_wide %>%
  # Get dataset names by extracting from column names
  {
    # Extract unique dataset names from AUROC columns
    auroc_cols <- grep("^AUROC_", names(.), value = TRUE)
    dataset_names <- sub("^AUROC_", "", auroc_cols)
    
    # Create combined columns
    for (dataset in dataset_names) {
      auroc_col <- paste0("AUROC_", dataset)
      subclass_col <- paste0("subclass_hit_", dataset)
      
      # Create new combined column
      . <- mutate(., 
                  !!dataset := paste0(.[[subclass_col]], 
                                      " (", 
                                      format(round(.[[auroc_col]], 2), nsmall = 2), 
                                      ")"))
    }
    .
  } %>%
  # Select only the subtype, new combined columns, and Avg_AUROC. 
  select(subtype, 
         all_of(other_studies), 
         Avg_AUROC)

# View the result
top_hits_combined %<>% arrange(Avg_AUROC)
print(top_hits_combined, n =30) 
top_hits_combined %>% mutate(Manual_agreement = "") %>% write_csv(paste0(output_folder, "top_hits.0.8.combined.slim.csv"))

###########
#show filtered out:
print(top_hits_combined %>% filter(!(subtype %in% subtypes_passing_filter)), n =30) 
####################

cell_counts <- obs  %>% filter(study_id == !!genesis_study_name) %>%
  dplyr::count(cell_type, name = "n_cells") %>% rename(subtype = cell_type)
top_hits_wide %<>% inner_join(cell_counts)
# 
# with_agreement <- top_hits_wide %>% mutate(
#   all_subclass_same = {
#     sc <- pick(contains("subclass"))
#     if (ncol(sc) <= 1) {
#       TRUE
#     } else {
#       out <- Reduce(`&`, lapply(sc[-1], function(x) x == sc[[1]]))
#       tidyr::replace_na(out, FALSE)
#     }
#   }) %>% filter(!all_subclass_same)
# #proportion of cells
# with_agreement %>% pull(n_cells) %>% sum() / nrow(obs %>% filter(study_id== genesis_study_name))

# top_hits_wide %>%
#   select(contains("subclass"), n_cells) %>%
#   group_by(across(contains("subclass"))) %>%
#   summarise(
#     n = dplyr::n(),                            # number of rows in this grouping (same as prior count())
#     total_cells = sum(n_cells, na.rm = TRUE),  # total cells across those rows
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n)) %>%
#   rename(!!genesis_study_name := n) %>%
#   write_csv(paste0(output_folder, "top_hits.0.8.collapsed.csv"))


# top_hits_wide %>%
#   select(contains("subclass"), starts_with("AUROC_"), n_cells) %>%
#   group_by(across(contains("subclass"))) %>%
#   summarise(
#     n = dplyr::n(),                          # number of rows in this grouping
#     total_cells = sum(n_cells, na.rm = TRUE),
#     
#     # mean AUROC per reference within this collapsed grouping
#     mean_AUROC_Ma_2022          = mean(AUROC_Ma_2022, na.rm = TRUE),
#     mean_AUROC_SEA_AD           = mean(AUROC_SEA_AD, na.rm = TRUE),
#     mean_AUROC_Siletti          = mean(AUROC_Siletti, na.rm = TRUE),
#     mean_AUROC_cxg_GEN_A1       = mean(AUROC_cxg_GEN_A1, na.rm = TRUE),
#     mean_AUROC_freeze2_GEN_A3_ctx = mean(AUROC_freeze2_GEN_A3_ctx, na.rm = TRUE),
#     
# 
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n)) %>%
#   rename(!!genesis_study_name := n) %>%
#   write_csv(paste0(output_folder, "top_hits.0.8.collapsed.csv"))

auroc_cols <- paste0("AUROC_", other_studies)

top_hits_collapsed <- top_hits_wide %>%
  select(subtype, contains("subclass"), any_of(auroc_cols), n_cells) %>%
  group_by(across(contains("subclass"))) %>%
  summarise(
    n = dplyr::n(),
    total_cells = sum(n_cells, na.rm = TRUE),
    
    # subtypes collapsed into this grouped row
    subtypes_collapsed = paste(unique(stats::na.omit(subtype)), collapse = ";"),
    
    across(
      any_of(auroc_cols),
      ~ mean(.x, na.rm = TRUE),
      .names = "mean_{.col}"
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    all_subclass_same = {
      sc <- pick(contains("subclass"))
      if (ncol(sc) <= 1) {
        TRUE
      } else {
        out <- Reduce(`&`, lapply(sc[-1], function(x) x == sc[[1]]))
        tidyr::replace_na(out, FALSE)
      }
    }
  ) %>%
  arrange(desc(n)) %>%
  rename(!!genesis_study_name := n) %>%
  select(-subtypes_collapsed, everything(), subtypes_collapsed)

top_hits_collapsed %>% filter(!all_subclass_same)

readr::write_csv(
  top_hits_collapsed,
  paste0(output_folder, "top_hits.0.8.collapsed.csv")
)


#######################
### AUROCs
#######################

all_AUROCs <- read_csv(paste0(base_folder, "aurocs_full.csv.gz"))
all_AUROCs %<>% rename(target = "...1")
dim(all_AUROCs)
#FILTERING
keep_ids <- union(
  as.character(subtypes_passing_filter_with_studyID),
  as.character(other_subclasses)
)

all_AUROCs %<>% filter(target %in% keep_ids) %>% select(target, any_of(keep_ids))
dim(all_AUROCs)
print("Filtered out subtypes")
all_AUROCs_matrix <- as.matrix(all_AUROCs %>% select(-target))
dim(all_AUROCs_matrix)
rownames(all_AUROCs_matrix) <- all_AUROCs %>% pull(target)
all_AUROCs_matrix <- all_AUROCs_matrix[sort(rownames(all_AUROCs_matrix)), sort(rownames(all_AUROCs_matrix))]
dim(all_AUROCs_matrix)

###########################
###########################
###########################
###########################
###########################



all_studies <- getStudyId(colnames(all_AUROCs_matrix))

mat0 <- all_AUROCs_matrix


#color
gamma <- 2.5  # >1 increases contrast near 0 and 1, while keeping a broader middle

warp_mid_01 <- function(x, gamma = 2) {
  x <- pmin(pmax(x, 0), 1)            # clamp to [0, 1]
  y <- 2 * x - 1                      # map to [-1, 1]
  y <- sign(y) * (abs(y) ^ gamma)     # warp around the midpoint
  (y + 1) / 2                         # map back to [0, 1]
}

col_base <- colorRamp2(
  c(0, 0.5, 1),
  c("blue", "white", "red")
)

col_fun <- function(x) col_base(warp_mid_01(x, gamma = gamma))


## ------------------------------------------------------------
## 1. Clustering on the filtered matrix
## ------------------------------------------------------------
hc   <- orderCellTypes(mat0)
dend <- as.dendrogram(hc)

## ------------------------------------------------------------
## 2. Define study IDs and cell-type labels
##    Use all studies EXCEPT GEN_A13 for block detection
## ------------------------------------------------------------
labs <- colnames(mat0)
study_id  <- getStudyId(labs)
cell_type <- sub("^[^|]+\\|", "", labs)   # everything after "study|"

is_non_a13 <- study_id != genesis_study_name

# These are the labels used for defining blocks
ref_cell_type <- ifelse(is_non_a13, cell_type, NA_character_)
names(ref_cell_type) <- labs

## ------------------------------------------------------------
## 3. Annotate dendrogram with the set of unique non-GEN_A13 cell types
##    contained in each subtree
## ------------------------------------------------------------
annotate_ref_types <- function(d, ref_labels) {
  if (is.leaf(d)) {
    lab <- attr(d, "label")
    types <- ref_labels[lab]
    if (is.na(types) || length(types) == 0L) {
      types <- character(0)
    } else {
      types <- unique(types)
    }
  } else {
    d[[1]] <- annotate_ref_types(d[[1]], ref_labels)
    d[[2]] <- annotate_ref_types(d[[2]], ref_labels)
    types <- unique(c(attr(d[[1]], "ref_types"), attr(d[[2]], "ref_types")))
  }
  
  attr(d, "ref_types") <- types
  attr(d, "n_ref_types") <- length(types)
  d
}

dend <- annotate_ref_types(dend, ref_cell_type)

## ------------------------------------------------------------
## 4. Collect maximal clades with exactly one unique non-GEN_A13 cell type
##    (parent must not also have exactly one)
## ------------------------------------------------------------
get_single_reftype_clades <- function(d, parent_n = NULL) {
  n <- attr(d, "n_ref_types")
  res <- list()
  
  if (n == 0L) {
    return(res)
  }
  
  if (n == 1L && (is.null(parent_n) || parent_n != 1L)) {
    return(list(labels(d)))
  }
  
  if (!is.leaf(d)) {
    res <- c(
      get_single_reftype_clades(d[[1]], n),
      get_single_reftype_clades(d[[2]], n)
    )
  }
  
  res
}

clade_leaves <- get_single_reftype_clades(dend)

## ------------------------------------------------------------
## 5. Turn clades into diagonal block ranges + block labels
## ------------------------------------------------------------
leaf_order <- labels(dend)
N          <- length(leaf_order)
leaf_index <- structure(seq_len(N), names = leaf_order)

block_ranges <- vector("list", length(clade_leaves))
block_labels <- character(length(clade_leaves))

for (i in seq_along(clade_leaves)) {
  seg <- clade_leaves[[i]]
  idx <- sort(leaf_index[seg])
  
  block_ranges[[i]] <- c(start = min(idx), end = max(idx))
  
  # unique non-GEN_A13 cell type inside this clade
  ct <- unique(na.omit(ref_cell_type[seg]))
  block_labels[i] <- if (length(ct) == 1L) ct else ""
}

## ------------------------------------------------------------
## 6. Reorder matrix to dendrogram order
## ------------------------------------------------------------
ord_names <- leaf_order
mat_ord   <- mat0[ord_names, ord_names, drop = FALSE]

## ------------------------------------------------------------
## 7. Study colours
## ------------------------------------------------------------

#study_colors <- c("blue", "red", "purple", "orange", "cyan")
study_colors <- c("#F27EA5", "#81CEF3", "#FBB454", "#C1E5AD", "cyan")

# Create named vector for study colors
study_color_map <- c(
  setNames("#C4E6B1", genesis_study_name),
  setNames(study_colors[1:length(other_studies)], other_studies)
)


studies <- names(study_color_map)

row_anno <- rowAnnotation(
  study = getStudyId(rownames(mat_ord)),
  col = list(study = study_color_map),
  show_annotation_name = FALSE
)

col_anno <- HeatmapAnnotation(
  study = getStudyId(colnames(mat_ord)),
  col = list(study = study_color_map), 
  show_legend = FALSE
)

## ------------------------------------------------------------
## 8. Heatmap
##    Turn row names off since block labels are drawn outside boxes
## ------------------------------------------------------------
ht <- Heatmap(
  mat_ord,
  name = "Cluster replicability\n(AUROC)",
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  
  show_row_names    = FALSE,
  show_column_names = FALSE,
  col = col_fun,
  top_annotation  = col_anno,
  left_annotation = row_anno
)

ht <- draw(ht)

## ------------------------------------------------------------
## 9. Draw diagonal boxes and place cell-type label to the right
## ------------------------------------------------------------
decorate_heatmap_body("Cluster replicability\n(AUROC)", {
  nr <- nrow(mat_ord)
  nc <- ncol(mat_ord)
  
  offset <- 0.01
  
  for (i in seq_along(block_ranges)) {
    r <- block_ranges[[i]]
    lab_text <- block_labels[i]
    
    if (is.na(lab_text) || lab_text == "") next
    
    r1 <- r["start"]
    r2 <- r["end"]
    
    c1 <- r1
    c2 <- r2
    
    x_left   <- (c1 - 1) / nc
    x_right  <- c2 / nc
    y_top    <- 1 - (r1 - 1) / nr
    y_bottom <- 1 - r2 / nr
    
    ## block border
    grid.rect(
      x = unit((x_left + x_right) / 2, "npc"),
      y = unit((y_bottom + y_top) / 2, "npc"),
      width  = unit(x_right - x_left, "npc"),
      height = unit(y_top - y_bottom, "npc"),
      gp = gpar(fill = NA, col = "black", lwd = 1)
    )
    
    ## label position
    x_text <- unit(min(x_right + offset, 0.98), "npc")
    y_text <- unit((y_bottom + y_top) / 2, "npc")
    
    ## exact text measurement
    label_gp <- gpar(cex = .9)
    label_grob <- textGrob(
      lab_text,
      x = x_text,
      y = y_text,
      just = "left",
      gp = label_gp
    )
    
    pad_x <- unit(1.1, "mm")
    pad_y <- unit(0.8, "mm")
    
    lab_w <- grobWidth(label_grob)
    lab_h <- grobHeight(label_grob)
    
    ## rounded white label box
    grid.roundrect(
      x = x_text + lab_w / 2,
      y = y_text,
      width  = lab_w + 2 * pad_x,
      height = lab_h + 2 * pad_y,
      r = unit(1.2, "mm"),
      just = "center",
      gp = gpar(fill = "white", col = "black", lwd = 0.8)
    )
    
    ## label text
    grid.draw(label_grob)
  }
})


###########################
###########################
###########################
###########################
###########################

rownames(all_AUROCs_matrix) <- all_AUROCs %>% pull(target)
all_AUROCs_matrix <- all_AUROCs_matrix[sort(rownames(all_AUROCs_matrix)), sort(rownames(all_AUROCs_matrix))]
dim(all_AUROCs_matrix)


row_anno <- rowAnnotation(
  study = getStudyId(rownames(all_AUROCs_matrix)),
  col = list(study = study_color_map)
)

# Create column annotation
col_anno <- HeatmapAnnotation(
  study = getStudyId(colnames(all_AUROCs_matrix)),
  col = list(study = study_color_map),
  show_legend = FALSE
)

# Create the heatmap with annotations
plot(1)
plotHeatmap(all_AUROCs_matrix)

orderCellTypes(all_AUROCs_matrix)
Heatmap(all_AUROCs_matrix, 
        row_order = rownames(all_AUROCs_matrix)[orderCellTypes(all_AUROCs_matrix)$order],
        column_order = colnames(all_AUROCs_matrix)[orderCellTypes(all_AUROCs_matrix)$order],
        name = "AUROC", 
        top_annotation = col_anno, 
        show_row_names = F,
        show_column_names = F,
        left_annotation = row_anno,
        column_names_gp = gpar(fontsize = 10),  # Adjust column text size
        row_names_gp = gpar(fontsize = 8)      # Adjust row text size
)

#show tree
plot(1)
hc <- orderCellTypes(all_AUROCs_matrix)

# Generate heatmap with dendrograms
(ht <- Heatmap(all_AUROCs_matrix, 
               cluster_rows = as.hclust(hc),  # Convert to hclust object
               cluster_columns = as.hclust(hc),  # Convert to hclust object
               name = "AUROC", 
               show_row_names = FALSE,
               show_column_names = FALSE,
               top_annotation = col_anno, 
               left_annotation = row_anno,
               column_names_gp = gpar(fontsize = 10),  # Adjust column text size
               row_names_gp = gpar(fontsize = 8)      # Adjust row text size
))
clustered_column_names <- colnames(all_AUROCs_matrix)[column_order(ht)]
#to see the ordering
clustered_column_names[!grepl(genesis_study_name, clustered_column_names, fixed = TRUE)]


########################
#1 vs 1
########################
########################
one_v_one_AUROCs <- read_csv(paste0(base_folder, "aurocs_1v1.csv.gz"))
one_v_one_AUROCs %<>% rename(target = "...1")
dim(one_v_one_AUROCs)
one_v_one_AUROCs %<>% filter(target %in% keep_ids) %>% select(target, any_of(keep_ids))
dim(one_v_one_AUROCs)
print("Filtered out subtypes")

one_v_one_AUROCs_matrix <- as.matrix(one_v_one_AUROCs %>% select(-target))

rownames(one_v_one_AUROCs_matrix) <- one_v_one_AUROCs %>% pull(target)
one_v_one_AUROCs_matrix <- one_v_one_AUROCs_matrix[sort(rownames(one_v_one_AUROCs_matrix)), sort(rownames(one_v_one_AUROCs_matrix))]
dim(one_v_one_AUROCs_matrix)

one_v_one_AUROCs

row_anno <- rowAnnotation(
  study = getStudyId(rownames(one_v_one_AUROCs_matrix)),
  col = list(study = study_color_map)
)

# Create column annotation
col_anno <- HeatmapAnnotation(
  study = getStudyId(colnames(one_v_one_AUROCs_matrix)),
  col = list(study = study_color_map),
  show_legend = FALSE
)

orderCellTypes(one_v_one_AUROCs_matrix)
plotHeatmap(one_v_one_AUROCs_matrix)
plot(1)
Heatmap(one_v_one_AUROCs_matrix, 
        row_order = rownames(one_v_one_AUROCs_matrix)[orderCellTypes(one_v_one_AUROCs_matrix)$order],
        column_order = colnames(one_v_one_AUROCs_matrix)[orderCellTypes(one_v_one_AUROCs_matrix)$order],
        name = "AUROC", 
        show_row_names = F,
        show_column_names = F,
        top_annotation = col_anno, 
        left_annotation = row_anno,
        column_names_gp = gpar(fontsize = 10),  # Adjust column text size
        row_names_gp = gpar(fontsize = 8)      # Adjust row text size
)

one_v_one_AUROCs_matrix_flat <- one_v_one_AUROCs_matrix[(rownames(one_v_one_AUROCs_matrix) %>% getStudyId()) != genesis_study_name, ]

row_anno <- rowAnnotation(
  study = getStudyId(rownames(one_v_one_AUROCs_matrix_flat)),
  col = list(study = study_color_map)
)

filtered_names <- rownames(one_v_one_AUROCs_matrix)[
  orderCellTypes(one_v_one_AUROCs_matrix)$order
] 
filtered_names <- filtered_names[filtered_names %in% rownames(one_v_one_AUROCs_matrix_flat)] 

Heatmap(one_v_one_AUROCs_matrix_flat, 
        row_order = filtered_names,
        column_order = colnames(one_v_one_AUROCs_matrix)[orderCellTypes(one_v_one_AUROCs_matrix)$order],
        name = "AUROC", 
        show_row_names = TRUE,
        show_column_names = FALSE,
        top_annotation = col_anno, 
        left_annotation = row_anno,
        column_names_gp = gpar(fontsize = 10),  # Adjust column text size
        row_names_gp = gpar(fontsize = 8)      # Adjust row text size
)


one_v_one_AUROCs_matrix_flat <- one_v_one_AUROCs_matrix[(rownames(one_v_one_AUROCs_matrix) %>% getStudyId()) != genesis_study_name, (rownames(one_v_one_AUROCs_matrix) %>% getStudyId()) == genesis_study_name]
dim(one_v_one_AUROCs_matrix_flat)
row_anno <- rowAnnotation(
  study = getStudyId(rownames(one_v_one_AUROCs_matrix_flat)),
  col = list(study = study_color_map)
)

filtered_names <- rownames(one_v_one_AUROCs_matrix)[
  orderCellTypes(one_v_one_AUROCs_matrix)$order
] 
filtered_names <- filtered_names[filtered_names %in% rownames(one_v_one_AUROCs_matrix_flat)] 

col_anno <- HeatmapAnnotation(
  study = getStudyId(colnames(one_v_one_AUROCs_matrix_flat)),
  col = list(study = study_color_map),
  show_legend = FALSE
)

ht <- Heatmap(one_v_one_AUROCs_matrix_flat, 
        row_order = filtered_names,
        #column_order = colnames(one_v_one_AUROCs_matrix_flat)[orderCellTypes(one_v_one_AUROCs_matrix)$order],
        cluster_columns = F,
        name = "AUROC", 
        show_row_names = TRUE,
        show_column_names = FALSE,
        top_annotation = col_anno, 
        left_annotation = row_anno,
        column_names_gp = gpar(fontsize = 10),  # Adjust column text size
        row_names_gp = gpar(fontsize = 8)      # Adjust row text size
)

#########

median_best_vs_next_AUROC <- one_v_one_AUROCs_matrix_flat[one_v_one_AUROCs_matrix_flat > 0.5] %>% median(na.rm=TRUE)

dim(one_v_one_AUROCs_matrix_flat)
one_v_one_AUROCs_matrix_flat[1:5,1:6]

# # Vector of column names meeting the criterion
# cols_keep <- colnames(one_v_one_AUROCs_matrix_flat)[
#   apply(one_v_one_AUROCs_matrix_flat, 2, function(x) {
#     x <- x[!is.na(x) & x != 0]          # only non-NA, non-zero
#     sum(x > 0.6) >= 2                   # require at least 2 values above 0.6 - L5_ET
#   })
# ]
# length(cols_keep)
dim(one_v_one_AUROCs_matrix_flat)

vals <- as.vector(one_v_one_AUROCs_matrix_flat)
vals <- vals[!is.na(vals) & vals > 0.5]

ggplot(data.frame(AUROC = vals), aes(x = AUROC)) +
  geom_density(kernel = 'rect') + theme_bw() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  labs(title = "Density of AUROC values > 0.5", x = "AUROC", y = "Density")

ggplot(data.frame(AUROC = vals), aes(x = AUROC)) +
  geom_histogram(
    breaks = seq(0.5, 1, length.out = 31),  # 30 bins ending exactly at 1
    fill = "black",
    color = "white",
    linewidth = 0.5,
    closed = "right"
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Best vs next AUROC",
    y = "Frequency"
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.margin = margin(t = 5.5, r = 14, b = 5.5, l = 5.5)
  )

############################
############################
############################
############################

# Copy used only for ordering (replace NA with 0)
mat_for_order <- one_v_one_AUROCs_matrix_flat
mat_for_order[is.na(mat_for_order)] <- 0

# Row order from NA->0 matrix
row_order_zero_na <- rownames(mat_for_order)[hclust(dist(mat_for_order), method = "average")$order]

# Column order from NA->0 matrix
col_order_zero_na <- colnames(mat_for_order)[hclust(dist(t(mat_for_order)), method = "average")$order]

Heatmap(
  one_v_one_AUROCs_matrix_flat,    # plot original matrix (keep NAs as-is)
  row_order = row_order_zero_na,   # order based on NA->0 version
  column_order = col_order_zero_na,
  name = "AUROC",
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = col_anno,
  left_annotation = row_anno,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 8)
)

############################

get_row_nonNA_col_AUROC_order <- function(
    mat,
    center = 0.5,
    band = 0.02,
    power = 2,
    pos_weight = 3,    # prioritize AUROC > 0.5 (red)
    neg_weight = 1,    # still consider AUROC < 0.5 (blue)
    top_left_bias = 0.10,  # nudges red-heavy columns left
    sign_tie = 0.02        # small sign-aware tie-break (keeps red/blue distinction without separating too much)
) {
  stopifnot(is.matrix(mat))
  
  # =========================
  # 1) ROW ORDER: non-NA overlap first (Jaccard), then reorder leaves by red-priority AUROC score
  # =========================
  present_rows <- !is.na(mat)
  
  inter_r <- tcrossprod(present_rows * 1L)  # row intersections
  n_pres_r <- rowSums(present_rows)
  union_r <- outer(n_pres_r, n_pres_r, "+") - inter_r
  
  sim_r <- matrix(1, nrow(mat), nrow(mat), dimnames = list(rownames(mat), rownames(mat)))
  has_union_r <- union_r > 0
  sim_r[has_union_r] <- inter_r[has_union_r] / union_r[has_union_r]
  diag(sim_r) <- 1
  
  row_order <- rownames(mat)
  if (nrow(mat) > 1) {
    hc_rows <- hclust(as.dist(1 - sim_r), method = "average")
    
    # AUROC-based row score for leaf reordering (preserve cluster, only flip branches)
    z_raw <- mat - center
    row_red  <- rowSums(pmax(z_raw - band, 0)^power, na.rm = TRUE)
    row_blue <- rowSums(pmax(-z_raw - band, 0)^power, na.rm = TRUE)
    
    row_score <- pos_weight * row_red - neg_weight * row_blue
    row_score[!is.finite(row_score)] <- 0
    
    d_row <- as.dendrogram(hc_rows)
    d_row <- stats::reorder(d_row, wts = setNames(-row_score, rownames(mat)), agglo.FUN = mean)
    row_order <- rownames(mat)[order.dendrogram(d_row)]
  }
  
  # Reorder matrix by rows before computing column AUROC refinement
  m <- mat[row_order, , drop = FALSE]
  
  # =========================
  # 2) COLUMN ORDER: non-NA overlap first (Jaccard), then reorder leaves by AUROC key
  #    using |AUROC - 0.5| for closeness (so pos/neg can stay close), with red priority.
  # =========================
  present_cols <- !is.na(m)
  
  inter_c <- crossprod(present_cols * 1L)   # column intersections
  n_pres_c <- colSums(present_cols)
  union_c <- outer(n_pres_c, n_pres_c, "+") - inter_c
  
  sim_c <- matrix(1, ncol(m), ncol(m), dimnames = list(colnames(m), colnames(m)))
  has_union_c <- union_c > 0
  sim_c[has_union_c] <- inter_c[has_union_c] / union_c[has_union_c]
  diag(sim_c) <- 1
  
  col_order <- colnames(m)
  
  # AUROC-derived refinement metrics (computed on row-ordered matrix)
  na_mask <- is.na(m)
  m0 <- m
  m0[na_mask] <- 0
  
  z <- m0 - center
  z[na_mask] <- 0  # IMPORTANT: original NAs contribute no signal
  
  # Positive/red and negative/blue deviations
  zp <- pmax(z - band, 0)^power
  zn <- pmax(-z - band, 0)^power
  
  # fallback if everything is inside the dead band
  if (all((zp + zn) == 0)) {
    zp <- pmax(z, 0)^power
    zn <- pmax(-z, 0)^power
  }
  
  # Symmetric signal for "where should this column sit?" (keeps red/blue close when in same rows)
  w_sym <- zp + zn
  
  # Red-priority signal for tie-breaking / left-bias
  w_red_priority <- pos_weight * zp + neg_weight * zn
  
  eps <- 1e-12
  row_pos <- seq_len(nrow(m))
  
  col_w_sym <- colSums(w_sym)
  col_center_sym <- as.numeric(crossprod(row_pos, w_sym)) / pmax(col_w_sym, eps)
  col_center_sym[col_w_sym == 0] <- Inf
  
  col_red  <- colSums(zp)
  col_blue <- colSums(zn)
  col_balance <- (col_red - col_blue) / pmax(col_red + col_blue, eps)
  col_balance[!is.finite(col_balance)] <- 0
  
  col_peak_red <- apply(zp, 2, max)
  col_peak_red[!is.finite(col_peak_red)] <- 0
  
  # Top-left red bias (small nudge left for red-heavy columns)
  max_red <- suppressWarnings(max(col_red, na.rm = TRUE))
  if (!is.finite(max_red) || max_red <= 0) max_red <- 1
  col_red_norm <- col_red / max_red
  
  col_key <- col_center_sym -
    top_left_bias * nrow(m) * col_red_norm -
    sign_tie * col_balance
  
  # If no signal at all in a column, keep it at the end
  col_key[col_w_sym == 0] <- Inf
  
  if (ncol(m) > 1) {
    hc_cols <- hclust(as.dist(1 - sim_c), method = "average")
    d_col <- as.dendrogram(hc_cols)
    
    # Reorder leaves by AUROC key while preserving non-NA clusters
    d_col <- stats::reorder(
      d_col,
      wts = setNames(col_key, colnames(m)),
      agglo.FUN = mean
    )
    
    col_order <- colnames(m)[order.dendrogram(d_col)]
    
    # Optional tiny final tie-break among exact-equal keys inside preserved order:
    # (doesn't usually matter, but can stabilize left bias)
    idx <- match(col_order, colnames(m))
    ord_df <- data.frame(
      col = col_order,
      key = col_key[idx],
      red = col_red[idx],
      blue = col_blue[idx],
      peak_red = col_peak_red[idx],
      stringsAsFactors = FALSE
    )
    # preserve dendrogram order unless key is essentially tied
    # (stable sort by current order is implicit)
    col_order <- ord_df$col[order(
      round(ord_df$key, 6),
      -ord_df$red,
      ord_df$blue,
      -ord_df$peak_red,
      na.last = TRUE
    )]
  }
  
  list(
    row_order = row_order,
    col_order = col_order
  )
}

mat_for_order <- one_v_one_AUROCs_matrix_flat[filtered_names, , drop = FALSE]

ord <- get_row_nonNA_col_AUROC_order(
  mat_for_order,
  center = 0.2,
  band = 0.02,
  power = 2,
  pos_weight = 2,
  neg_weight = 1,
  top_left_bias = 0.012,
  sign_tie = 0.02
)

#col_keep_flag <- colnames(one_v_one_AUROCs_matrix_flat) %in% cols_keep

col_anno <- HeatmapAnnotation(
  study = getStudyId(colnames(one_v_one_AUROCs_matrix_flat)),
  col = list(
    study = study_color_map
  ),
  show_legend = FALSE, show_annotation_name = FALSE
)

row_anno <- rowAnnotation(
  study = getStudyId(rownames(one_v_one_AUROCs_matrix_flat)),
  col = list(study = study_color_map),
  show_annotation_name = FALSE
)

Heatmap(
  one_v_one_AUROCs_matrix_flat,   # plot original matrix (keep NAs as-is)
  row_order = ord$row_order,
  column_order = ord$col_order,
  name = "AUROC",
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = col_anno,
  row_labels = getCellType(rownames(one_v_one_AUROCs_matrix_flat)),
  left_annotation = row_anno,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 8)
)

################
###############
# annotated
#############

mat <- one_v_one_AUROCs_matrix_flat

row_order <- if (is.character(ord$row_order)) {
  match(ord$row_order, rownames(mat))
} else {
  ord$row_order
}

col_order <- if (is.character(ord$col_order)) {
  match(ord$col_order, colnames(mat))
} else {
  ord$col_order
}

displayed_row_labels <- getCellType(rownames(mat))
row_labels_ord <- displayed_row_labels[row_order]
mat_ord <- mat[row_order, col_order, drop = FALSE]

col_anno <- HeatmapAnnotation(
  study = getStudyId(colnames(mat)),
  col = list(study = study_color_map),
  show_legend = FALSE,
  show_annotation_name = FALSE
)

row_anno <- rowAnnotation(
  study = getStudyId(rownames(mat)),
  col = list(study = study_color_map),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat,
  row_order = row_order,
  column_order = col_order,
  name = "AUROC",
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = col_anno,
  row_labels = displayed_row_labels,
  left_annotation = row_anno,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 8)
)

ht <- draw(ht)

decorate_heatmap_body("AUROC", {
  nr <- nrow(mat_ord)
  nc <- ncol(mat_ord)
  
  row_runs <- rle(row_labels_ord)
  run_ends <- cumsum(row_runs$lengths)
  run_starts <- run_ends - row_runs$lengths + 1
  
  for (i in seq_along(row_runs$values)) {
    lab_text <- row_runs$values[i]
    
    if (is.na(lab_text) || lab_text == "" || row_runs$lengths[i] < 2) {
      next
    }
    
    r1 <- run_starts[i]
    r2 <- run_ends[i]
    
    block_mat <- mat_ord[r1:r2, , drop = FALSE]
    
    signal_cols <- which(
      colSums(!is.na(block_mat) & block_mat > 0.5) > 0
    )
    
    if (length(signal_cols) == 0) {
      next
    }
    
    c1 <- min(signal_cols)
    c2 <- max(signal_cols)
    
    x_left <- (c1 - 1) / nc
    x_right <- c2 / nc
    y_top <- 1 - (r1 - 1) / nr
    y_bottom <- 1 - r2 / nr
    
    grid.rect(
      x = unit((x_left + x_right) / 2, "npc"),
      y = unit((y_bottom + y_top) / 2, "npc"),
      width = unit(x_right - x_left, "npc"),
      height = unit(y_top - y_bottom, "npc"),
      gp = gpar(fill = NA, col = "black", lwd = 1)
    )
    
    x_text <- unit(min(x_right + 0.01, 0.98), "npc")
    y_text <- unit((y_bottom + y_top) / 2, "npc")
    
    label_gp <- gpar(cex = 0.9)
    
    label_grob <- textGrob(
      lab_text,
      x = x_text,
      y = y_text,
      just = "left",
      gp = label_gp
    )
    
    pad_x <- unit(1.1, "mm")
    pad_y <- unit(0.8, "mm")
    
    grid.roundrect(
      x = x_text + grobWidth(label_grob) / 2,
      y = y_text,
      width = grobWidth(label_grob) + 2 * pad_x,
      height = grobHeight(label_grob) + 2 * pad_y,
      r = unit(1.2, "mm"),
      just = "center",
      gp = gpar(fill = "white", col = "black", lwd = 0.8)
    )
    
    grid.draw(label_grob)
  }
})





############################
######################
##write out one-v-one - don't filter -reload
######################

one_v_one_AUROCs <- read_csv(paste0(base_folder, "aurocs_1v1.csv.gz"))
one_v_one_AUROCs %<>% rename(target = "...1")
one_v_one_AUROCs_matrix <- as.matrix(one_v_one_AUROCs %>% select(-target))
rownames(one_v_one_AUROCs_matrix) <- one_v_one_AUROCs %>% pull(target)
one_v_one_AUROCs_matrix <- one_v_one_AUROCs_matrix[sort(rownames(one_v_one_AUROCs_matrix)), sort(rownames(one_v_one_AUROCs_matrix))]

all_rows <- tibble()
for(target_type in target_types) {
  target_hits <- tibble(name = rownames(one_v_one_AUROCs_matrix), auc = one_v_one_AUROCs_matrix[, paste0(genesis_study_name, "|", target_type)]) %>% filter(auc > 0.5)  %>% arrange(name)
  target_hits %<>% mutate(reference_name = target_type)
  target_hits %<>% mutate(study_target = getStudyId(name))
  target_hits%<>% mutate(name = getCellType(name))
  target_hits %<>%
    pivot_wider(
      names_from = study_target,
      values_from = c(name, auc),
      id_cols = reference_name
    )
  all_rows %<>% bind_rows(target_hits)
}
all_rows %>% head(20)

all_rows %<>% mutate(
  Avg_AUROC = rowMeans(select(., starts_with("auc_")), na.rm = TRUE)
) %>%  arrange(Avg_AUROC)
#remove it's own 1-1
all_rows %<>% select(-!!paste0("name_", genesis_study_name), -!!paste0("auc_", genesis_study_name))
all_rows %>% mutate(Manual_agreement = "") %>% write_csv(paste0(output_folder, "one_v_one.top_hits.combined.csv"))


all_rows_combined <- all_rows %>%
  {
    # Extract unique dataset names from auc columns
    auc_cols <- grep("^auc_", names(.), value = TRUE)
    dataset_names <- sub("^auc_", "", auc_cols)
    
    # Create combined columns
    for (dataset in dataset_names) {
      auc_col <- paste0("auc_", dataset)
      name_col <- paste0("name_", dataset)
      
      # Create new combined column
      . <- mutate(., 
                  !!dataset := paste0(.[[name_col]], 
                                      " (", 
                                      format(round(.[[auc_col]], 2), nsmall = 2), 
                                      ")"))
    }
    .
  } %>%
  # Calculate average AUROC across all datasets
  mutate(
    Avg_AUROC = rowMeans(select(., starts_with("auc_")), na.rm = TRUE)
  ) %>%
  # Sort by average AUROC in descending order
  arrange(Avg_AUROC) %>%
  # Select only the reference_name, new combined columns, and Avg_AUROC
  select(reference_name, 
         any_of(other_studies),
         Avg_AUROC)

all_rows_combined %>% pull(Avg_AUROC) %>% median()
all_rows_combined %>% write_csv(paste0(output_folder, "one_v_one.top_hits.combined.slim.csv"))



all_rows_combined
filtered_out <- inner_join(top_hits_combined %>% filter(!(subtype %in% subtypes_passing_filter)),
                           all_rows_combined %>% rename_with(~ paste0("best_v_next_", .x), .cols = -1) %>% rename(subtype=reference_name))
#########################
#show filtered out:
print(filtered_out %>% select(-best_v_next_Avg_AUROC), n =30) 
################

######## output summary text 
summary_text <- paste(
  paste0(genesis_study_name, " versus GENA1-A3 cortical subclasses"),
  paste0(fmt_n(obs %>% filter(study_id == genesis_study_name) %>% nrow()), " cells"),
  paste0(fmt_n(nrow(var)), " HVGs"),
  paste0(obs %>% filter(study_id == genesis_study_name) %>% pull(cell_type) %>% unique() %>% length(),
         " supertype level 3 clusters from Donghoon"),
  paste0(length(subtypes_without_tophits), " without any AUROC > 0.8 hits: ", paste(subtypes_without_tophits, collapse=", ")),
  paste0(
    "Missing subclasses: ",
    if (length(missing_subclasses) == 0) {
      "all subclasses appear"
    } else {
      paste(missing_subclasses, collapse = ", ")
    }
  ),
  paste0(
    signif(proportion_subtypes_passing_filter * 100, 3),
    "% pass filter (",
    signif(proportion_cells_passing_filter * 100, 3),
    "% of cells pass filter)"
  ),
  paste0("Median AUROCs after filtering:"),
  paste0(top_hits_wide %>% filter(subtype %in% subtypes_passing_filter) %>% pull(Avg_AUROC) %>% median() %>% signif(digits =3), " for all versus all top hits"),
  paste0(median_best_vs_next_AUROC %>% signif(digits =3), " for best vs next"),
  paste0(basename(base_folder), "/", basename(output_folder)),
  sep = "\n"
)

cat(summary_text)

#write as csv:
summary_stats <- tibble(
  base_dir = basename(base_folder),
  comparison = paste0(genesis_study_name, " versus GENA1-A3 cortical subclasses"),
  n_cells = obs %>% filter(study_id == genesis_study_name) %>% nrow(),
  n_hvgs = nrow(var),
  n_lvl3_clusters = obs %>%
    filter(study_id == genesis_study_name) %>%
    pull(cell_type) %>%
    unique() %>%
    length(),
  n_nohit_subtypes = length(subtypes_without_tophits),
  nohit_subtypes = paste(subtypes_without_tophits, collapse = ", "),
  pct_subtypes_pass = proportion_subtypes_passing_filter * 100,
  pct_cells_pass = proportion_cells_passing_filter * 100,
  median_tophit_auroc = top_hits_wide %>%
    filter(subtype %in% subtypes_passing_filter) %>%
    pull(Avg_AUROC) %>%
    median(),
  median_best_next_auroc = median_best_vs_next_AUROC
)


write_csv(
  summary_stats,
  paste0(output_folder, "summary_stats_from_R.csv")
)


#show table of those that fail - all vs all wide?



###############
#graph
###############
one_v_one_AUROCs <- read_csv(paste0(base_folder, "aurocs_1v1.csv.gz"))
one_v_one_AUROCs %<>% rename(target = "...1")
dim(one_v_one_AUROCs)
one_v_one_AUROCs %<>% filter(target %in% keep_ids) %>% select(target, any_of(keep_ids))
#one_v_one_AUROCs %>% filter(target == "GEN_A4_PFC|3_2_3" ) %>% as.data.frame()
length(keep_ids)
dim(one_v_one_AUROCs)
print("Filtered out subtypes")

one_v_one_AUROCs_matrix <- as.matrix(one_v_one_AUROCs %>% select(-target))

rownames(one_v_one_AUROCs_matrix) <- one_v_one_AUROCs %>% pull(target)
one_v_one_AUROCs_matrix <- one_v_one_AUROCs_matrix[sort(rownames(one_v_one_AUROCs_matrix)), sort(rownames(one_v_one_AUROCs_matrix))]
dim(one_v_one_AUROCs_matrix)


cluster_graph = makeClusterGraph(one_v_one_AUROCs_matrix, low_threshold = 0.6)
#cluster_graph = makeClusterGraph(one_v_one_AUROCs_matrix)
plot(1)

plotClusterGraph(cluster_graph, size_factor=3) #, legend_cex = 0)


# 1. Extract study from node names
study <- sub("\\|.*", "", V(cluster_graph)$name)

# 2. Assign colors by study using gg_color
unique_studies <- unique(study)
gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors <- gg_color(length(unique_studies))
study_colors <- setNames(colors, unique_studies)
V(cluster_graph)$color <- study_colors[study]

# 3. Set labels: hide for genesis_study_name
vertex_labels <- V(cluster_graph)$name
vertex_labels[study == genesis_study_name] <- NA
vertex_labels[study != "GEN_A1"] <- NA
vertex_labels <- getCellType(vertex_labels)

# 4. Set node sizes: smaller for GEN_A2
default_size <- 5
small_size <- 3
vertex_sizes <- ifelse(study == genesis_study_name, small_size, default_size)

# 5. Plot
size_factor <- 1
label_cex <- 0.6 * size_factor
vertex_sizes <- vertex_sizes * size_factor
vertex_colors <- setNames(colors, unique_studies)

#layout <- layout_with_fr(cluster_graph, niter = 1000, area = vcount(graph)^2.5)
layout <- layout_nicely(cluster_graph) *4
plot(cluster_graph, layout = layout,
     vertex.color = V(cluster_graph)$color,
     vertex.label = vertex_labels,
     vertex.label.cex = label_cex,
     vertex.label.font = 2,
     vertex.label.color = "black", 
     vertex.size = vertex_sizes,
     vertex.frame.color = NA,
     edge.arrow.size = 0.1 * size_factor,
     edge.arrow.width = 0.5* size_factor)
plot(1)
graphics::legend("bottomright",
                 legend = names(vertex_colors),
                 pt.bg = vertex_colors,
                 pch = 21,
                 title = "Study")

############

g_tbl <- as_tbl_graph(cluster_graph) %>%               # convert once
  mutate(
    study        = sub("\\|.*", "", name),             # extract study
    label        = getCellType(name),                  # initial label
    label        = ifelse(study == genesis_study_name, NA, label),  # hide for GEN_A2
    label        = ifelse(study == "GEN_A5" | study == "GEN_A1", label, NA),  # hide for GEN_A2
    size         = ifelse(study == genesis_study_name, 2, 3)     # node size rule
  )

## ---------- 2. Colour palette identical to your gg_color() helper ----------
gg_color <- function(n){
  hues <- seq(15, 375, length.out = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
study_levels  <- sort(unique(g_tbl %>% pull(study)))
study_pal     <- setNames(gg_color(length(study_levels)), study_levels)


#TEMP
#TEMP
#TEMP
# study_map <- c(
#   "GEN_A1"     = "PsychAD",
#   "GEN_A2"     = "ROSMAP",
#   "GEN_A3_ctx" = "AMP-PD",
#   "GEN_A13"    = "Mathys 2023"
# )
# g_tbl <- g_tbl %>%
#   activate(nodes) %>%
#   mutate(
#     study = dplyr::coalesce(unname(study_map[study]), study)
#   )

# study_colors <- c("#F27EA5", "#81CEF3", "#FBB454", "#C1E5AD", "cyan")
# # Create named vector for study colors
# study_color_map <- c(
#   setNames("#C4E6B1", genesis_study_name),
#   setNames(study_colors[1:length(other_studies)], other_studies)
# )

#study_color_map <- study_color_map[match(names(study_map), names(study_color_map))]
#names(study_color_map) <- study_map[names(study_color_map)]
#end TEMP
#end TEMP
#end TEMP
#end TEMP
#end TEMP


set.seed(1)
layout_df <- create_layout(g_tbl, layout = "nicely")

# Combine layout with graph node data
layout_df$label <- layout_df$label  # already included, but explicitly added for clarity


ggraph(layout_df) +
  geom_edge_link(
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(1.5, "mm"),
    start_cap = circle(1.5, "mm"),
    colour = alpha("grey60", 0.6),
    width = 0.3
  ) +
  geom_node_point(
    aes(x = x, y = y, fill = study, size = size),
    shape = 21, stroke = ifelse(study == genesis_study_name, .25, .5)#, colour = ifelse(study == genesis_study_name, "grey60", NA)
  )  +
  ggrepel::geom_text_repel(
    data = layout_df[!is.na(layout_df$label), ],
    aes(x = x, y = y, label = label),
    size = 4.5,
    box.padding = 0.35, #change to move it farther
    point.padding = 0.3,
    segment.size = 0.2,
    max.overlaps = Inf,
    color = "black",
    bg.color = "white",
    bg.r = 0.15
  ) +
  scale_fill_manual(values = study_color_map, name = "Study") +
  scale_size_identity() +
  theme_graph(base_family = "sans") +
  #theme(legend.position = "bottom") + 
  guides(
    fill = guide_legend(
      override.aes = list(size = 6)  # increase circle size in legend
    )
  )

################
## use both edges of the best vs next
##################
cluster_graph <- makeClusterGraph(one_v_one_AUROCs_matrix) 

edge_threshold <- 0.6

g_tbl <- as_tbl_graph(cluster_graph) %>%
  activate(nodes) %>%
  mutate(
    study = sub("\\|.*", "", name),
    label = getCellType(name),
    label = ifelse(study == genesis_study_name, NA, label),
    label = ifelse(study %in% c("GEN_A1", "GEN_A5"), label, NA),
    size = ifelse(study == genesis_study_name, 2, 3)
  ) %>%
  activate(edges) %>%
  mutate(
    strong_edge = !is.na(weight) & weight > edge_threshold
  )

set.seed(1)
#layout_df <- create_layout(g_tbl, layout = "nicely")
layout_df <- create_layout(
  g_tbl,
  layout = "fr",
  weights = weight
)

ggraph(layout_df) +
  geom_edge_link(
    aes(filter = !strong_edge),
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(1.5, "mm"),
    start_cap = circle(1.5, "mm"),
    colour = alpha("grey90", 0.35),
    width = 0.25
  ) +
  geom_edge_link(
    aes(filter = strong_edge),
    arrow = arrow(length = unit(2, "mm"), type = "closed"),
    end_cap = circle(1.5, "mm"),
    start_cap = circle(1.5, "mm"),
    colour = alpha("grey60", 0.6),
    width = 0.3
  ) +
  geom_node_point(
    aes(fill = study, size = size),
    shape = 21,
    stroke = ifelse(layout_df$study == genesis_study_name, 0.25, 0.5)
  ) +
  ggrepel::geom_text_repel(
    data = layout_df[!is.na(layout_df$label), ],
    aes(x = x, y = y, label = label),
    size = 4.5,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.size = 0.2,
    max.overlaps = Inf,
    color = "black",
    bg.color = "white",
    bg.r = 0.15
  ) +
  scale_fill_manual(values = study_color_map, name = "Study") +
  scale_size_identity() +
  theme_graph(base_family = "sans") +
  guides(
    fill = guide_legend(
      override.aes = list(size = 6)
    )
  )


cat(summary_text)
