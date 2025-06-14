#!/usr/bin/env Rscript

# Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(tidyr)) # For pivot_longer
suppressPackageStartupMessages(library(argparse))

# --- Argument Parsing ---
parser <- ArgumentParser(description='Plot distributions of polygenicity metrics.')
parser$add_argument('input_polygenicity_table', type='character', help='Path to the input polygenicity TSV file.')
parser$add_argument('--output_dir_base', type='character', default='figures/polygenicity_distributions', help='Base directory for output plots.')

args <- parser$parse_args()

# --- Configuration ---
input_filename_base <- tools::file_path_sans_ext(basename(args$input_polygenicity_table))
filename_parts <- unlist(strsplit(input_filename_base, "\\."))
annotation <- filename_parts[length(filename_parts)]
output_dir <- file.path(args$output_dir_base)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# --- Load and Prepare Data ---
polygenicity_data_raw <- read_tsv(args$input_polygenicity_table, col_types = cols(.default = "c")) %>%
  mutate(across(ends_with("_polygenicity") | ends_with("_se"), as.numeric))

# Identify metric columns (those ending in _polygenicity) and preserve their original order
# Exclude _se columns
original_metric_cols <- names(polygenicity_data_raw)[stringr::str_ends(names(polygenicity_data_raw), "_polygenicity") & !stringr::str_ends(names(polygenicity_data_raw), "_se")]

if (length(original_metric_cols) == 0) {
  stop("No metric columns ending with '_polygenicity' (and not '_se') found in the input file.")
}

# Get unique datasets to iterate over
datasets_to_plot <- unique(polygenicity_data_raw$dataset)

for (current_dataset_name in datasets_to_plot) {
  message(paste0("Processing dataset: ", current_dataset_name))

  current_dataset_data <- polygenicity_data_raw %>%
    filter(dataset == current_dataset_name)

# Cleaned metric names, maintaining original order for factor levels
  cleaned_metric_names_ordered <- original_metric_cols %>%
  stringr::str_remove("_polygenicity$") %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_to_title()

# Pivot to long format
  long_data_dataset <- current_dataset_data %>%
  tidyr::pivot_longer(
    cols = all_of(original_metric_cols),
    names_to = "metric_type_raw",
    values_to = "value"
  ) %>%
  dplyr::filter(!is.na(value), value > 0) # Filter out NA/non-positive values before log transform

if (nrow(long_data_dataset) == 0) {
    message(paste("Dataset:", current_dataset_name, "- No positive data available for plotting after filtering. Skipping plot for this dataset."))
    next # Skip to next dataset
  }

# Create metric_type factor with levels based on original column order (reversed for ggplot y-axis)
  long_data_dataset <- long_data_dataset %>%
  dplyr::mutate(
    metric_type_cleaned = stringr::str_remove(metric_type_raw, "_polygenicity$"),
    metric_type_cleaned = stringr::str_replace_all(metric_type_cleaned, "_", " "),
    metric_type_cleaned = stringr::str_to_title(metric_type_cleaned),
    metric_type = factor(metric_type_cleaned, levels = rev(cleaned_metric_names_ordered))
  )

# --- Plotting ---
# Determine breaks and limits for the x-axis (log10 scale)
  # Data 'value' is on original scale here.
  min_data_val <- min(long_data_dataset$value, na.rm = TRUE)
  max_data_val <- max(long_data_dataset$value, na.rm = TRUE)

# Define plot limits on the log10 scale, ensuring 1 (10^0) is included
plot_limit_min_log10 <- floor(log10(min(1, min_data_val)))
plot_limit_max_log10 <- ceiling(log10(max_data_val))

x_axis_breaks_powers <- seq(from = plot_limit_min_log10, to = plot_limit_max_log10, by = 1)
x_axis_labels <- 10^(x_axis_breaks_powers)

  p_dist <- ggplot2::ggplot(long_data_dataset, ggplot2::aes(x = value, y = metric_type)) +
  ggplot2::geom_boxplot(outlier.colour = "gray60", outlier.size = 0.5, outlier.alpha = 0.5) +
  ggplot2::scale_x_log10(
    breaks = x_axis_labels, # Breaks are on the original scale for scale_x_log10
    labels = x_axis_labels,
    limits = c(10^plot_limit_min_log10, 10^plot_limit_max_log10)
  ) +
  ggplot2::labs(
    title = paste0("Polygenicity Metric Distributions: ", tools::toTitleCase(current_dataset_name), " (", annotation, ")"),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.9)),
    # plot.title is NULL, so no need to format it
    panel.grid.major.x = ggplot2::element_line(linetype = "dashed", color = "gray80"),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank()
  )

# Define output filename
  dataset_sanitized <- stringr::str_replace_all(tolower(current_dataset_name), "[^a-z0-9_]+", "_")
  output_filename_loop <- sprintf("%s.%s.%s.pdf", input_filename_base, dataset_sanitized, tolower(annotation))
  output_filepath_loop <- file.path(output_dir, output_filename_loop)

  # Save plot
  num_metrics <- dplyr::n_distinct(long_data_dataset$metric_type)
  plot_height_dist <- max(3, num_metrics * 0.5 + 1.5)
  plot_width_dist <- 7

  ggplot2::ggsave(output_filepath_loop, plot = p_dist, width = plot_width_dist, height = plot_height_dist, units = "in")
  message("Saved distribution plot to: ", output_filepath_loop)
}

message("All polygenicity distribution plots generated.")
