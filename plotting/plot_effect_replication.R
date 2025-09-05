#!/usr/bin/env Rscript

# --- Libraries ---
# suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(gridExtra))
# suppressPackageStartupMessages(library(ggrepel))
source('plotting/plotting_functions.R')
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))

# --- Argument Parsing ---
parser <- ArgumentParser(description = 'Plot effect size replication figures from a CSV file.')
parser$add_argument('input_file', help = 'Path to the input CSV file with effect replication data.')
args <- parser$parse_args()

# --- Configuration ---
# Save alongside input top directory: top/tables/... -> top/figures/...
input_filename_base <- tools::file_path_sans_ext(basename(args$input_file))
input_dir <- dirname(args$input_file)
top_dir <- dirname(input_dir)  # parent of 'tables'
figures_dir <- file.path(top_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
output_dir <- figures_dir


# Color palette for emphasized traits
emphasized_traits <- get_emphasized_traits()

# --- Load Data ---
raw_effects <- read_tsv(args$input_file, col_types = cols(.default = "c"))

effects <- raw_effects %>%
  mutate(
    trait = abbreviation,
    ancestry_group = stringr::str_extract(replication_dataset, "[^_]+$"),
    count = as.numeric(count),
    expected = as.numeric(expected_posterior_mean),
    expected_se = as.numeric(expected_posterior_mean_se),
    mean_aou = as.numeric(observed_effect_replication),
    mean_aou_se = as.numeric(observed_effect_replication_se),
    mean_ukbb = as.numeric(mean_effect_primary),
    mean_ukbb_se = as.numeric(mean_effect_primary_se)
  ) %>%
  select(
    trait, ancestry_group, bin, count, annotation, # Keep original and new identifiers
    expected, expected_se,      
    mean_aou, mean_aou_se,      
    mean_ukbb, mean_ukbb_se     
  )

# --- Main Script ---

# 2. Effect Size Plots
for (anc in unique(effects$ancestry_group)) {
  # Filter data for current ancestry
  eff_anc <- effects %>% filter(ancestry_group == anc)
  
  # Plot 1: Observed vs Posterior-Mean Effect Size
  p1 <- plot_with_emphasis(
    data = eff_anc,
    x_var = "expected",
    y_var = "mean_aou",
    x_lower_var = "expected_se",
    x_upper_var = "expected_se",
    y_lower_var = "mean_aou_se",
    y_upper_var = "mean_aou_se",
    color_var = "trait",
    label_var = "trait",
    emphasize_traits = names(emphasized_traits),
    colors = emphasized_traits,
    title = paste("Observed vs Posterior-Mean Effect:", toupper(anc)),
    xlab = "Posterior Mean Effect Size (95% CI)",
    ylab = "Observed Effect Size (95% CI)",
    make_square_axes = TRUE,
    add_y_equals_x_line = TRUE,
    show_non_emphasized = TRUE,
    errorbar_scale = 1.96,
    point_size = 0,
    emphasized_point_size = 1
  )
  
  # Plot 2: Observed AOU vs Observed UKBB Effect Size
  p2 <- plot_with_emphasis(
    data = eff_anc,
    x_var = "mean_ukbb",
    y_var = "mean_aou",
    x_lower_var = "mean_ukbb_se",
    x_upper_var = "mean_ukbb_se",
    y_lower_var = "mean_aou_se",
    y_upper_var = "mean_aou_se",
    color_var = "trait",
    label_var = "trait",
    emphasize_traits = names(emphasized_traits),
    colors = emphasized_traits,
    title = paste("AOU vs UKBB Effect Sizes:", toupper(anc)),
    xlab = "UKBB Effect Size (95% CI)",
    ylab = "AOU Effect Size (95% CI)",
    make_square_axes = TRUE,
    add_y_equals_x_line = TRUE,
    errorbar_scale = 1.96,
    point_size = 0,
    emphasized_point_size = 1
  )
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, ncol = 2, 
                              top = paste("Effect Size Replication:", toupper(anc)))
  
  # Save combined plot
  ggsave(
    filename = file.path(output_dir, paste0("effect_size_comparison_", anc, ".pdf")),
    plot = combined_plot,
    width = 16,
    height = 8
  )
}

# --- 4. Meta-analyzed Effect Size Plots ---

# Define ancestry colors for these specific meta-plots
# Ensure these names match the actual values in effects$ancestry_group
# (e.g., "EUR", "AFR", "EAS", "AMR", "SAS")
# User requested "eur", "afr", "eas" to be emphasized. Data shows "eur", "afr", "amr".
# Emphasizing "eur", "afr" as per request and data presence.
meta_ancestry_colors <- c(
  "eur" = "#1f77b4", # Blue
  "afr" = "#ff7f0e", # Orange
  "amr" = "#d62728", # Red
  "eas" = "#2ca02c", # Green (kept in case it appears)
  "sas" = "#9467bd", # Purple (kept in case it appears)
  "oth" = "grey50"   # Default for other groups if any
  # Add/modify as per actual ancestry groups in your data to ensure all are colored
)
emphasized_ancestries_for_meta <- c("eur", "afr", "amr") # Values from 'ancestry_group' to emphasize

# Aggregate data: 
# 1. Create a bin_rank within each trait & ancestry_group
# 2. Group by ancestry_group and bin_rank for meta-analysis
effects_with_rank <- effects %>%
  group_by(trait, ancestry_group) %>%
  # Assuming data is already ordered by bin as desired for ranking
  mutate(bin_rank = row_number()) %>%
  ungroup()

effects_meta_summary <- effects_with_rank %>%
  group_by(ancestry_group, bin_rank) %>%
  summarise(
    mean_expected = mean(expected, na.rm = TRUE),
    sd_expected = sd(expected, na.rm = TRUE),
    mean_mean_ukbb = mean(mean_ukbb, na.rm = TRUE),
    sd_mean_ukbb = sd(mean_ukbb, na.rm = TRUE),
    mean_mean_aou = mean(mean_aou, na.rm = TRUE),
    sd_mean_aou = sd(mean_aou, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(ancestry_group, bin_rank) # Ensure correct order for lines

# You might want to print unique ancestries to confirm color mapping:
# print(paste("Unique ancestries for meta effect plot:", paste(unique(effects_meta_summary$ancestry_group), collapse=", ")))

# Plot 1: Observed AoU vs. Expected (Meta-analyzed)
p_meta_eff_aou_vs_exp <- plot_with_emphasis(
  data = effects_meta_summary,
  x_var = "mean_expected",
  y_var = "mean_mean_aou",
  x_lower_var = "sd_expected", 
  x_upper_var = "sd_expected", 
  y_lower_var = "sd_mean_aou", 
  y_upper_var = "sd_mean_aou", 
  color_var = "ancestry_group",
  label_var = NULL, # Use legend for labels
  make_square_axes = TRUE,
  emphasize_traits = c("eur", "afr", "amr"), # Explicitly list for meta-plot
  colors = meta_ancestry_colors, 
  title = NULL,
  xlab = "Posterior-mean effect (UKBB) (s.e.)",
  ylab = "Mean effect (AoU) (s.e.)",
  point_size = 0.8, 
  emphasized_point_size = 1.2, 
  add_y_equals_x_line = TRUE,
  use_legend_for_labels = TRUE,
  legend_title_for_color = NULL
)


# Plot 2: Observed AoU vs. Observed UKBB (Meta-analyzed)
p_meta_eff_aou_vs_ukbb <- plot_with_emphasis(
  data = effects_meta_summary,
  x_var = "mean_mean_ukbb",
  y_var = "mean_mean_aou",
  x_lower_var = "sd_mean_ukbb",
  x_upper_var = "sd_mean_ukbb",
  y_lower_var = "sd_mean_aou",
  y_upper_var = "sd_mean_aou",
  color_var = "ancestry_group",
  label_var = NULL, # Labels will come from legend
  make_square_axes = TRUE,
  emphasize_traits = emphasized_ancestries_for_meta,
  colors = meta_ancestry_colors,
  title = NULL,
  xlab = "Mean effect (UKBB) (s.e.)",
  ylab = "Mean effect (AoU) (s.e.)",
  point_size = 0.8,
  emphasized_point_size = 1.2,
  add_y_equals_x_line = TRUE,
  use_legend_for_labels = TRUE,
  legend_title_for_color = NULL
)

# Combine the two meta-analyzed plots into a single figure and save once
combined_meta <- grid.arrange(
  p_meta_eff_aou_vs_exp,
  p_meta_eff_aou_vs_ukbb,
  ncol = 2
)

ggsave(
  filename = file.path(figures_dir, paste0(input_filename_base, ".pdf")),
  plot = combined_meta,
  width = 16,
  height = 8
)

cat("Saved plot to:", file.path(figures_dir, paste0(input_filename_base, ".pdf")), "\n")

