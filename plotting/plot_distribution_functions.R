#!/usr/bin/env Rscript

# --- Libraries ---
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(argparse))

# --- Argument Parsing ---
parser <- ArgumentParser(description='Plot heritability explained by top N genes for different traits and datasets.')
parser$add_argument('input_distribution_table', type='character', help='Path to the input TSV file with distribution data.')
parser$add_argument('--output_dir_base', type='character', default='figures/distributions', help='Base directory for output plots.')
parser$add_argument('--dataset', type='character', default=NULL, help='Optional: Specific dataset to plot. If NULL, all datasets are plotted.')

args <- parser$parse_args()

# --- Configuration ---
input_filename_base <- tools::file_path_sans_ext(basename(args$input_distribution_table))
filename_parts <- unlist(strsplit(input_filename_base, "\\."))
annotation <- filename_parts[length(filename_parts)] # Assumes annotation is the last part, e.g., pLoF

output_dir <- file.path(args$output_dir_base, annotation)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ---
message(paste("Loading data from:", args$input_distribution_table))
tryCatch({
  all_data_raw <- read_tsv(args$input_distribution_table, col_types = cols(.default = "c")) %>%
    mutate(
      H2_Percent = as.numeric(H2_Percent),
      Num_Genes = as.numeric(Num_Genes),
      xmin_genes = as.numeric(xmin_genes),
      xmax_genes = as.numeric(xmax_genes)
      # Num_Genes_SE is not explicitly used in the original plotting logic for lines/points
    )
}, error = function(e) {
  stop(paste("Error loading or processing input file:", args$input_distribution_table, "Error:", e$message))
})

if (nrow(all_data_raw) == 0) {
  stop("Input file is empty or no data could be loaded.")
}

# --- Determine Datasets to Plot ---
unique_datasets_in_file <- unique(all_data_raw$dataset)

datasets_to_plot <- if (!is.null(args$dataset)) {
  if (args$dataset %in% unique_datasets_in_file) {
    args$dataset
  } else {
    warning(paste("Specified dataset '", args$dataset, "' not found in the input file. Available datasets: ", paste(unique_datasets_in_file, collapse=", "), ". Plotting all datasets instead.", sep=""))
    unique_datasets_in_file
  }
} else {
  unique_datasets_in_file
}

if (length(datasets_to_plot) == 0) {
  message("No datasets to plot after filtering.")
  quit(status = 0)
}

# --- Plotting Loop (per dataset) ---
for (current_dataset_name in datasets_to_plot) {
  message(paste0("\nProcessing dataset: ", current_dataset_name))
  
  dataset_data <- all_data_raw %>%
    filter(dataset == current_dataset_name)
  
  if (nrow(dataset_data) == 0) {
    message(paste0("No data for dataset: ", current_dataset_name, " after filtering. Skipping."))
    next
  }
  
  dataset_sanitized <- stringr::str_replace_all(tolower(current_dataset_name), "[^a-z0-9_]+", "_")
  output_pdf_filename <- sprintf("%s.%s.%s.pdf", input_filename_base, dataset_sanitized, tolower(annotation))
  output_pdf_path <- file.path(output_dir, output_pdf_filename)
  
  tryCatch({
    pdf(output_pdf_path, width = 10, height = 8)
    message(paste0("Opened PDF device for dataset ", current_dataset_name, ": ", output_pdf_path))
  }, error = function(e) {
    warning(paste0("Error opening PDF device for ", current_dataset_name, ": ", output_pdf_path, ". Error: ", e$message, ". Skipping this dataset."))
    next
  })
  
  unique_abbreviations <- unique(dataset_data$abbreviation)
  if (all(is.na(unique_abbreviations))) unique_abbreviations <- "UnknownTrait" # Handle case where all abbreviations are NA
  
  for (current_trait_abbreviation in unique_abbreviations) {
    trait_plot_data <- dataset_data %>%
      filter(is.na(abbreviation) | abbreviation == current_trait_abbreviation) %>%
      arrange(H2_Percent) # Ensure data is ordered for line plot if necessary
    
    if (nrow(trait_plot_data) == 0 || all(is.na(trait_plot_data$Num_Genes)) || all(is.na(trait_plot_data$H2_Percent))) {
      message(paste0("No valid data for trait '", current_trait_abbreviation, "' in dataset '", current_dataset_name, "'. Skipping plot."))
      next
    }
    
    has_ribbon <- all(c("xmin_genes", "xmax_genes") %in% names(trait_plot_data)) &&
                  any(!is.na(trait_plot_data$xmin_genes) & !is.na(trait_plot_data$xmax_genes))
    
    num_total_genes <- max(trait_plot_data$Num_Genes, na.rm = TRUE)
    if (is.infinite(num_total_genes) || is.na(num_total_genes)) num_total_genes <- 1000 # Fallback
    
    x_axis_limits <- c(1, max(1, num_total_genes)) # Ensure max limit is at least 1
    manual_breaks <- c(1, 10, 100, 1000, 10000, 100000)
    actual_breaks <- manual_breaks[manual_breaks <= x_axis_limits[2] & manual_breaks >= x_axis_limits[1]]
    if (x_axis_limits[1] == 1 && !(1 %in% actual_breaks)) {
      actual_breaks <- c(1, actual_breaks)
    }
    actual_breaks <- sort(unique(actual_breaks))
    if(length(actual_breaks) == 0 && x_axis_limits[1] >= 1) actual_breaks <- x_axis_limits[1]
    if(length(actual_breaks) == 0) actual_breaks <- 1 # Ultimate fallback for breaks

    trait_display_name <- ifelse(is.na(current_trait_abbreviation) || current_trait_abbreviation == "UnknownTrait", "Overall", current_trait_abbreviation)
    plot_title_text <- paste0("Heritability Explained by Top N Genes - ", trait_display_name, " (", current_dataset_name, ", ", annotation, ")")
    
    p <- ggplot(trait_plot_data, aes(x = Num_Genes, y = H2_Percent)) +
      geom_line(color = "seagreen", na.rm = TRUE) +
      geom_point(color = "seagreen", na.rm = TRUE) +
      scale_x_log10(limits = x_axis_limits, 
                   breaks = actual_breaks, 
                   labels = label_comma(accuracy = 1)) + 
      scale_y_continuous(limits = c(0, 100), 
                        breaks = seq(0, 100, by = 10)) +
      labs(title = plot_title_text,
           x = "Number of Genes (log scale)",
           y = "% Heritability Explained") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
    if (has_ribbon) {
        ribbon_data <- trait_plot_data %>%
            filter(!is.na(xmin_genes) & !is.na(xmax_genes) & !is.na(H2_Percent))
        if(nrow(ribbon_data) > 0) {
            p <- p + geom_ribbon(data = ribbon_data, aes(xmin = xmin_genes, xmax = xmax_genes), 
                                fill = "skyblue", 
                                alpha = 0.5, 
                                na.rm = TRUE)
        } else {
            message(paste0("Not enough valid ribbon data for trait: ", current_trait_abbreviation, " in dataset: ", current_dataset_name))
        }
    } else {
      message(paste0("No ribbon data (xmin_genes, xmax_genes) for trait: ", current_trait_abbreviation, " in dataset: ", current_dataset_name))
    }
    
    print(p)
    message(paste0("  Generated plot for trait: ", current_trait_abbreviation))
  }
  
  dev.off()
  message(paste0("Plots for dataset ", current_dataset_name, " saved to: ", output_pdf_path))
}

message("\nAll distribution plots generated.")
