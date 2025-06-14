library(dplyr)
library(readr)
library(stringr)
source("plotting/plotting_functions.R")
suppressPackageStartupMessages(library(argparse))

# --- Argument Parsing ---
parser <- ArgumentParser(description = 'Plot heritability comparison figures from a TSV file.')
parser$add_argument('input_heritability_table', help = 'Path to the input TSV file with heritability estimates.')
parser$add_argument('--output_dir_base', default = "figures/heritability_comparison", help = 'Base directory for output figures (default: figures/heritability_comparison).')
parser$add_argument('--dataset1_name', default = "genebass", help = 'Name of the first dataset for comparison (default: genebass).')
args <- parser$parse_args()

# --- Configuration ---
input_filename_base <- tools::file_path_sans_ext(basename(args$input_heritability_table))
filename_parts <- unlist(strsplit(input_filename_base, "\\."))
annotation <- filename_parts[length(filename_parts)]
output_dir <- file.path(args$output_dir_base, annotation)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Color palette for emphasized traits
emphasized_traits <- get_emphasized_traits()

# Load data
heritability_data_raw <- read_tsv(args$input_heritability_table, col_types = cols(.default = "c")) %>%
  mutate(
    across(c(total_h2, total_h2_se, positive_h2, positive_h2_se, 
             negative_h2, negative_h2_se, feature1_h2, feature1_h2_se, 
             feature2_h2, feature2_h2_se, feature3_h2, feature3_h2_se, 
             feature4_h2, feature4_h2_se, feature5_h2, feature5_h2_se), 
           as.numeric)
  )

# Process Dataset 1 data (e.g., genebass)
dataset1_h2_data <- heritability_data_raw %>%
  filter(dataset == args$dataset1_name) %>%
  select(abbreviation, h2_observed_dataset1 = total_h2, h2_se_dataset1 = total_h2_se)

# Determine Dataset 2 list (all datasets in heritability_data_raw except dataset1_name)
aou_datasets_to_plot <- heritability_data_raw %>%
  filter(dataset != args$dataset1_name) %>%
  distinct(dataset) %>%
  pull(dataset)

input_filename_base <- tools::file_path_sans_ext(basename(args$input_heritability_table))
dataset1_sanitized <- str_replace_all(tolower(args$dataset1_name), "[^a-z0-9_]+", "_")

for (dataset2_name in aou_datasets_to_plot) {
  dataset2_h2_data <- heritability_data_raw %>%
    filter(dataset == dataset2_name) %>%
    select(abbreviation, h2_observed_dataset2 = total_h2, h2_se_dataset2 = total_h2_se)
  
  current_plot_data <- inner_join(dataset1_h2_data, dataset2_h2_data, by = "abbreviation")

  if (nrow(current_plot_data) == 0) {
    message(paste0("Skipping ", dataset2_name, " vs ", args$dataset1_name, " (", annotation, ") due to no common traits after filtering."))
    next
  }
  
  dataset2_sanitized <- str_replace_all(tolower(dataset2_name), "[^a-z0-9_]+", "_")
  
  plot_title <- paste0("Heritability Comparison: ", tools::toTitleCase(args$dataset1_name), " vs ", tools::toTitleCase(dataset2_name), " (", annotation, ")")
  xlab_label <- paste0(tools::toTitleCase(args$dataset1_name), " Heritability (h²)")
  ylab_label <- paste0(tools::toTitleCase(dataset2_name), " Heritability (h²)")
  
  p <- plot_with_emphasis(
    data = current_plot_data,
    x_var = "h2_observed_dataset1",
    y_var = "h2_observed_dataset2",
    x_lower_var = "h2_se_dataset1",
    x_upper_var = "h2_se_dataset1",
    y_lower_var = "h2_se_dataset2",
    y_upper_var = "h2_se_dataset2",
    color_var = "abbreviation", 
    label_var = "abbreviation",
    emphasize_traits = names(emphasized_traits),
    colors = emphasized_traits,
    title = plot_title,
    xlab = xlab_label,
    ylab = ylab_label,
    add_y_equals_x_line = TRUE,
    make_square_axes = TRUE,
    errorbar_scale = 1, 
    show_non_emphasized = TRUE,
    point_size = 1.5,
    emphasized_point_size = 3
  )
  
  output_filename <- sprintf("%s.%s.vs.%s.pdf", input_filename_base, dataset2_sanitized, dataset1_sanitized)
  output_filepath <- file.path(output_dir, output_filename)
  
  ggsave(output_filepath, p, width = 8, height = 7)
  message(paste0("Saved plot for ", dataset2_name, " vs ", args$dataset1_name, " (", annotation, ") to ", output_filepath))
}

cat("All heritability comparison plots generated.\n")
