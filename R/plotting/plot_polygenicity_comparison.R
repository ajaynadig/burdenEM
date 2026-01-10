#!/usr/bin/env Rscript

# Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(argparse))

# --- Argument Parsing ---
parser <- ArgumentParser(description='Compare polygenicity metrics across two tables.')
parser$add_argument('input_polygenicity_table1', type='character', help='Path to first polygenicity TSV file (source for effective_polygenicity).')
parser$add_argument('input_polygenicity_table2', type='character', help='Path to second polygenicity TSV file (source for effective_mutational_polygenicity).')
parser$add_argument('--output_dir_base', type='character', default='figures/polygenicity', help='Base directory for output plots.')
args <- parser$parse_args()

# --- Configuration ---
input_filename_base1 <- tools::file_path_sans_ext(basename(args$input_polygenicity_table1))
filename_parts1 <- unlist(strsplit(input_filename_base1, "\\."))
annotation1 <- filename_parts1[length(filename_parts1)]
output_dir <- file.path(args$output_dir_base)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Processing: ", args$input_polygenicity_table1, " vs ", args$input_polygenicity_table2)

# --- Load data ---
d1 <- read_tsv(args$input_polygenicity_table1, col_types = cols(.default = "c")) %>%
  mutate(across(ends_with("_polygenicity") | ends_with("_se"), as.numeric))

d2 <- read_tsv(args$input_polygenicity_table2, col_types = cols(.default = "c")) %>%
  mutate(across(ends_with("_polygenicity") | ends_with("_se"), as.numeric))

# Datasets present in both tables
datasets_to_plot <- intersect(unique(d1$dataset), unique(d2$dataset))

for (current_dataset_name in datasets_to_plot) {
  message(paste0("Processing dataset: ", current_dataset_name))

  eff <- d1 %>%
    filter(dataset == current_dataset_name) %>%
    transmute(
      abbreviation,
      dataset,
      value = effective_polygenicity,
      se = effective_polygenicity_se
    ) %>%
    filter(!is.na(value) & !is.na(se) & value > 0 & se <= value) %>%
    mutate(type = "Effective")

  mut <- d2 %>%
    filter(dataset == current_dataset_name) %>%
    transmute(
      abbreviation,
      dataset,
      value = effective_mutational_polygenicity,
      se = effective_mutational_polygenicity_se
    ) %>%
    filter(!is.na(value) & !is.na(se) & value > 0 & se <= value) %>%
    mutate(type = "Mutational")

  common_traits <- intersect(eff$abbreviation, mut$abbreviation)
  eff <- eff %>% filter(abbreviation %in% common_traits)
  mut <- mut %>% filter(abbreviation %in% common_traits)

  plot_data <- bind_rows(eff, mut) %>%
    mutate(
      log10_value = log10(value),
      log10_se = se / (value * log(10))
    )

  if (nrow(plot_data) == 0) {
    message(paste("Dataset:", current_dataset_name, "- no overlapping traits after filtering. Skipping."))
    next
  }

  # Order traits by effective polygenicity
  eff_order <- plot_data %>% filter(type == "Effective") %>% arrange(log10_value) %>% pull(abbreviation)
  plot_data <- plot_data %>% mutate(abbreviation = factor(abbreviation, levels = eff_order))

  # Axis limits and breaks based on all points in this dataset
  min_data_log10 <- min(plot_data$log10_value, na.rm = TRUE)
  max_data_log10 <- max(plot_data$log10_value, na.rm = TRUE)
  plot_limit_min_log10 <- floor(min(0, min_data_log10))
  plot_limit_max_log10 <- ceiling(max_data_log10)
  custom_breaks <- seq(from = plot_limit_min_log10, to = plot_limit_max_log10, by = 1)
  custom_labels <- 10^(custom_breaks)

  p <- ggplot(plot_data, aes(x = log10_value, y = abbreviation, color = type)) +
    geom_errorbarh(aes(xmin = log10_value - log10_se, xmax = log10_value + log10_se), height = 0.1, linewidth = 0.5) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c(Effective = "#1f77b4", Mutational = "#d62728"), name = NULL) +
    scale_x_continuous(breaks = custom_breaks, labels = custom_labels, limits = c(plot_limit_min_log10, plot_limit_max_log10)) +
    labs(x = "Polygenicity (s.e.)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", color = "gray80"),
      axis.text.y = element_text(size = rel(0.8)),
      axis.line.y = element_line(),
      axis.ticks.y = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.position = "right"
    )

  # Output
  dataset_sanitized <- stringr::str_replace_all(tolower(current_dataset_name), "[^a-z0-9_]+", "_")
  output_filename_loop <- sprintf("%s_vs_%s.%s.%s.pdf", input_filename_base1, tools::file_path_sans_ext(basename(args$input_polygenicity_table2)), dataset_sanitized, tolower(annotation1))
  output_filepath_loop <- file.path(output_dir, output_filename_loop)

  num_traits <- dplyr::n_distinct(plot_data$abbreviation)
  plot_height <- max(3, num_traits * 0.18 + 1.5)
  plot_width <- 6

  ggsave(output_filepath_loop, plot = p, width = plot_width, height = plot_height, units = "in")
  message("Saved plot to: ", output_filepath_loop)
}

message("All polygenicity comparison plots generated.")
