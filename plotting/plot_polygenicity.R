#!/usr/bin/env Rscript

# Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(forcats)) # For reordering factors
suppressPackageStartupMessages(library(argparse))

# --- Argument Parsing ---
parser <- ArgumentParser(description='Plot polygenicity metrics.')
parser$add_argument('input_polygenicity_table', type='character', help='Path to the input polygenicity TSV file.')
parser$add_argument('--output_dir_base', type='character', default='figures/polygenicity', help='Base directory for output plots.')

args <- parser$parse_args()

# --- Configuration ---
input_filename_base <- tools::file_path_sans_ext(basename(args$input_polygenicity_table))
filename_parts <- unlist(strsplit(input_filename_base, "\\."))
annotation <- filename_parts[length(filename_parts)]
output_dir <- file.path(args$output_dir_base)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# --- Plotting ---
message("Processing: ", args$input_polygenicity_table)

# Load data
polygenicity_data_raw <- read_tsv(args$input_polygenicity_table, col_types = cols(.default = "c")) %>%
  mutate(across(ends_with("_polygenicity") | ends_with("_se"), as.numeric))


# Check for required columns
required_cols_polygenicity <- c("abbreviation", "dataset", "effective_polygenicity", "effective_polygenicity_se")
if (!all(required_cols_polygenicity %in% names(polygenicity_data_raw))) {
  stop(paste("File:", args$input_polygenicity_table, "- missing one or more required columns:", paste(required_cols_polygenicity, collapse=", ")))
}


# Get unique datasets to iterate over
datasets_to_plot <- unique(polygenicity_data_raw$dataset)

for (current_dataset_name in datasets_to_plot) {
  message(paste0("Processing dataset: ", current_dataset_name))

  # Filter data for the current dataset, join with pheno_info, calculate log10 values and SE
  plot_data_dataset <- polygenicity_data_raw %>%
    filter(dataset == current_dataset_name) %>%
  dplyr::filter(
    !is.na(effective_polygenicity) & 
    !is.na(effective_polygenicity_se) &
    effective_polygenicity_se <= effective_polygenicity &
    effective_polygenicity > 0 # Crucial for log and division
  ) %>%
  dplyr::mutate(
    log10_effective_polygenicity = log10(effective_polygenicity),
    # Delta method for SE of log10(X): SE(log10(X)) approx SE(X) / (X * ln(10))
    log10_effective_polygenicity_se = effective_polygenicity_se / (effective_polygenicity * log(10))
  )

if (nrow(plot_data_dataset) == 0) {
    message(paste("Dataset:", current_dataset_name, "in file:", args$input_polygenicity_table, "- no data left after filtering or all effective_polygenicity are zero/negative. Skipping plot for this dataset."))
    next # Skip to the next dataset
  }

# Order traits by log10_effective_polygenicity for plotting
  plot_data_dataset <- plot_data_dataset %>%
    dplyr::mutate(abbreviation = forcats::fct_reorder(abbreviation, log10_effective_polygenicity, .desc = FALSE))

# Determine breaks, labels, and limits for the x-axis
# Data (log10_effective_polygenicity) is already log10 transformed.
# We want to ensure '1' (which is 0 on this log10 scale) is visible.
min_data_log10 <- min(plot_data_dataset$log10_effective_polygenicity, na.rm = TRUE)
max_data_log10 <- max(plot_data_dataset$log10_effective_polygenicity, na.rm = TRUE)

# Define plot limits on the log10 scale
plot_limit_min_log10 <- floor(min(0, min_data_log10))
plot_limit_max_log10 <- ceiling(max_data_log10)

custom_breaks <- seq(from = plot_limit_min_log10, to = plot_limit_max_log10, by = 1)
custom_labels <- 10^(custom_breaks)
  
# Create plot
    p <- ggplot2::ggplot(plot_data_dataset, ggplot2::aes(x = log10_effective_polygenicity, y = abbreviation)) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = log10_effective_polygenicity - log10_effective_polygenicity_se, 
                 xmax = log10_effective_polygenicity + log10_effective_polygenicity_se),
    height = 0.1, 
    color = "gray40",
    linewidth = 0.5
  ) +
  ggplot2::geom_point(size = 1.5, color = "black") +
  ggplot2::scale_x_continuous(breaks = custom_breaks, labels = custom_labels, limits = c(plot_limit_min_log10, plot_limit_max_log10)) +
  ggplot2::labs(
    title = paste0("Effective Polygenicity: ", tools::toTitleCase(current_dataset_name), " (", annotation, ")"), 
    x = "Effective Polygenicity (SE)", 
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(linetype = "dashed", color = "gray80"),
    axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.8)),
    axis.line.y = ggplot2::element_line(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# Define output filename
  dataset_sanitized <- stringr::str_replace_all(tolower(current_dataset_name), "[^a-z0-9_]+", "_")
  output_filename_loop <- sprintf("%s.%s.%s.pdf", input_filename_base, dataset_sanitized, tolower(annotation))
  output_filepath_loop <- file.path(output_dir, output_filename_loop)

  # Save plot
    num_traits <- dplyr::n_distinct(plot_data_dataset$abbreviation)
  plot_height <- max(3, num_traits * 0.18 + 1.5)
  plot_width <- 6

  ggplot2::ggsave(output_filepath_loop, plot = p, width = plot_width, height = plot_height, units = "in")
  message("Saved plot to: ", output_filepath_loop)
}

message("All polygenicity plots generated.")