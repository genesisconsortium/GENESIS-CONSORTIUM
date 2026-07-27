library(data.table)
library(ggplot2)

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 02_qc_pearson_residual_files.R <residual_dir>")
}

RESIDUAL_DIR <- args[1]
OUTDIR <- file.path(RESIDUAL_DIR, "residual_qc")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(
  RESIDUAL_DIR,
  pattern = "_(class|subclass)_.*Pearson\\.tsv\\.gz$",
  full.names = TRUE
)


class_color_map <- c(
  EN = "#197EC0",
  IN = "#1A9993",
  Astro = "#D2AF81",
  Oligo = "#D5E4A2",
  OPC = "#FED439",
  Immune = "#C80813",
  Vascular = "#FD7446"
)

subclass_color_map <- c(
  "EN_L2-3_IT" = "#659EC7",
  "EN_L3-5_IT_1" = "#95B9C7",
  "EN_L3-5_IT_2" = "#6495ED",
  "EN_L3-5_IT_3" = "#79BAEC",
  "EN_L5-6_NP" = "#0020C2",
  "EN_L5_ET" = "#4863A0",
  "EN_L6B" = "#488AC7",
  "EN_L6_CT" = "#3BB9FF",
  "EN_L6_IT_1" = "#B4CFEC",
  "EN_L6_IT_2" = "#1589FF",
  "IN_CGE_KCNG1" = "#4E8975",
  "IN_CGE_LAMP5_RELN" = "#3B9C9C",
  "IN_CGE_VIP" = "#89C35C",
  "IN_LAMP5_LHX6" = "#7BCCB5",
  "IN_MGE_PVALB" = "#B2C248",
  "IN_MGE_PVALB_CHC" = "#008080",
  "IN_MGE_SST" = "#728C00",
  Astro = "#C19A6B",
  Oligo = "#ECE5B6",
  OPC = "#FFF380",
  Micro_PVM = "#F75D59",
  Adaptive_Immune = "#6F4E37",
  Endo = "#FFA62F",
  PC = "#E0B0FF",
  VLMC = "#B93B8F"
)


if (length(files) == 0) {
  stop("No Pearson residual files found in: ", RESIDUAL_DIR)
}

count_lines_gz <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))

  n <- 0
  repeat {
    x <- readLines(con, n = 100000)
    if (length(x) == 0) break
    n <- n + length(x)
  }

  n
}

get_level <- function(x) {
  if (grepl("_subclass_", basename(x))) return("subclass")
  if (grepl("_class_", basename(x))) return("class")
  NA
}

get_cell_type <- function(x) {
  b <- basename(x)

  if (grepl("_subclass_", b)) {
    b <- sub(".*_subclass_", "", b)
  } else {
    b <- sub(".*_class_", "", b)
  }

  b <- sub("_residualsPearson\\.tsv\\.gz$", "", b)
  b
}

qc <- data.frame()
resid_sample <- data.frame()

for (f in files) {
  message("QC: ", basename(f))

  this_level <- get_level(f)
  this_cell_type <- get_cell_type(f)

  dt_head <- fread(f, nrows = 1)
  n_samples <- ncol(dt_head) - 1
  n_genes <- count_lines_gz(f) - 1

  tmp <- fread(f, nrows = 300)
  vals <- as.numeric(unlist(tmp[, -1, with = FALSE]))
  vals <- vals[is.finite(vals)]

  if (length(vals) > 50000) {
    vals <- sample(vals, 50000)
  }

  resid_sample <- rbind(
    resid_sample,
    data.frame(
      value = vals,
      level = this_level,
      cell_type = this_cell_type
    )
  )

  qc <- rbind(
    qc,
    data.frame(
      file = basename(f),
      level = this_level,
      cell_type = this_cell_type,
      n_genes = n_genes,
      n_samples = n_samples,
      n_values_sampled = length(vals)
    )
  )
}

write.csv(
  qc,
  file.path(OUTDIR, "pearson_residual_file_qc.csv"),
  row.names = FALSE
)

p1 <- ggplot(qc, aes(x = level, y = n_genes, fill = level)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  theme_bw() +
  labs(
    title = "Genes per Pearson residual file",
    x = "",
    y = "Genes"
  ) +
  theme(legend.position = "none")

ggsave(file.path(OUTDIR, "genes_per_pearson_file.pdf"), p1, width = 6, height = 4)
ggsave(file.path(OUTDIR, "genes_per_pearson_file.png"), p1, width = 6, height = 4, dpi = 300)

p2 <- ggplot(qc, aes(x = level, y = n_samples, fill = level)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  theme_bw() +
  labs(
    title = "Samples per Pearson residual file",
    x = "",
    y = "Samples"
  ) +
  theme(legend.position = "none")

ggsave(file.path(OUTDIR, "samples_per_pearson_file.pdf"), p2, width = 6, height = 4)
ggsave(file.path(OUTDIR, "samples_per_pearson_file.png"), p2, width = 6, height = 4, dpi = 300)

class_resid <- resid_sample[resid_sample$level == "class", ]
subclass_resid <- resid_sample[resid_sample$level == "subclass", ]

if (nrow(class_resid) > 0) {

class_resid$cell_type <- factor(class_resid$cell_type, levels = names(class_color_map))

p_class <- ggplot(class_resid, aes(x = value, fill = cell_type)) +
  geom_histogram(bins = 45, color = NA) +
  facet_wrap(~ cell_type, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = class_color_map, drop = FALSE) +
  theme_bw() +
  labs(
    title = "Pearson residual distribution by class",
    x = "Pearson residual",
    y = "Sampled values"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 8)
  )

ggsave(file.path(OUTDIR, "pearson_residual_histogram_by_class.pdf"), p_class, width = 14, height = 3.5)
ggsave(file.path(OUTDIR, "pearson_residual_histogram_by_class.png"), p_class, width = 14, height = 3.5, dpi = 300)


}

if (nrow(subclass_resid) > 0) {

subclass_resid$cell_type <- factor(subclass_resid$cell_type, levels = names(subclass_color_map))

p_subclass <- ggplot(subclass_resid, aes(x = value, fill = cell_type)) +
  geom_histogram(bins = 45, color = NA) +
  facet_wrap(~ cell_type, nrow = 3, scales = "free_y") +
  scale_fill_manual(values = subclass_color_map, drop = FALSE) +
  theme_bw() +
  labs(
    title = "Pearson residual distribution by subclass",
    x = "Pearson residual",
    y = "Sampled values"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 7)
  )

ggsave(file.path(OUTDIR, "pearson_residual_histogram_by_subclass.pdf"), p_subclass, width = 16, height = 9)
ggsave(file.path(OUTDIR, "pearson_residual_histogram_by_subclass.png"), p_subclass, width = 16, height = 9, dpi = 300)

}

class_color_map <- c(
  EN = "#197EC0",
  IN = "#1A9993",
  Astro = "#D2AF81",
  Oligo = "#D5E4A2",
  OPC = "#FED439",
  Immune = "#C80813",
  Vascular = "#FD7446"
)

subclass_color_map <- c(
  "EN_L2-3_IT" = "#659EC7",
  "EN_L3-5_IT_1" = "#95B9C7",
  "EN_L3-5_IT_2" = "#6495ED",
  "EN_L3-5_IT_3" = "#79BAEC",
  "EN_L5-6_NP" = "#0020C2",
  "EN_L5_ET" = "#4863A0",
  "EN_L6B" = "#488AC7",
  "EN_L6_CT" = "#3BB9FF",
  "EN_L6_IT_1" = "#B4CFEC",
  "EN_L6_IT_2" = "#1589FF",
  "IN_CGE_KCNG1" = "#4E8975",
  "IN_CGE_LAMP5_RELN" = "#3B9C9C",
  "IN_CGE_VIP" = "#89C35C",
  "IN_LAMP5_LHX6" = "#7BCCB5",
  "IN_MGE_PVALB" = "#B2C248",
  "IN_MGE_PVALB_CHC" = "#008080",
  "IN_MGE_SST" = "#728C00",
  Astro = "#C19A6B",
  Oligo = "#ECE5B6",
  OPC = "#FFF380",
  Micro_PVM = "#F75D59",
  Adaptive_Immune = "#6F4E37",
  Endo = "#FFA62F",
  PC = "#E0B0FF",
  VLMC = "#B93B8F"
)


class_colors <- data.frame(
  level = "class",
  cell_type = names(class_color_map),
  color_hex = unname(class_color_map)
)

subclass_colors <- data.frame(
  level = "subclass",
  cell_type = names(subclass_color_map),
  color_hex = unname(subclass_color_map)
)

write.csv(class_colors, file.path(OUTDIR, "class_color_map.csv"), row.names = FALSE)
write.csv(subclass_colors, file.path(OUTDIR, "subclass_color_map.csv"), row.names = FALSE)

print(qc)

cat("\nSaved QC outputs to:\n", OUTDIR, "\n")
