#!/usr/bin/env Rscript

# --- Libraries ---
# suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr)) # Added for str_extract
# suppressPackageStartupMessages(library(ggrepel))
source('plotting/plotting_functions.R')

# Helper function for binomial confidence intervals using the Clopper-Pearson exact method
binomial_ci <- function(x, n, conf_level = 0.95) {
    # Ensure x is not NA and n is not NA and n > 0
    if (is.na(x) || is.na(n) || n == 0) {
        return(data.frame(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    # Ensure x is an integer (number of successes)
    x_int <- round(x) # if x is a rate, x*n should be the number of successes
    if (x_int < 0) x_int <- 0
    if (x_int > n) x_int <- n

    alpha <- 1 - conf_level
    lower <- qbeta(alpha/2, x_int, n - x_int + 1)
    upper <- qbeta(1 - alpha/2, x_int + 1, n - x_int)
    # estimate from the function is x_int/n, which might differ slightly from original observed_rate due to rounding
    # we will keep the original observed_rate for consistency in plots, and use these CIs
    data.frame(estimate_binom = x_int/n, lower_ci_binom = lower, upper_ci_binom = upper)
}

# --- Argument Parsing ---
parser <- ArgumentParser(description = "Plot P-value replication rates from a results table.")
parser$add_argument("input_file", help = "Path to the input TSV file containing replication data.")
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
rates_raw <- read_tsv(args$input_file, col_types = cols(.default = "c")) # Read all as char first to handle mixed types if any, then convert

rates <- rates_raw %>%
  rename(
    trait = abbreviation,
    ancestry_group_full = replication_dataset,
    count_char = n_genes_tested_replication,
    expected_rate_char = expected_replication_rate,
    observed_rate_char = observed_replication_rate
  ) %>%
  mutate(
    ancestry_group = str_extract(ancestry_group_full, "[^_]+$"), # Extract part after last underscore
    count = as.integer(count_char),
    expected_rate = as.numeric(expected_rate_char),
    observed_rate = as.numeric(observed_rate_char)
  ) %>%
  # Calculate x for binomial_ci (number of observed successes)
  # If observed_rate or count is NA, x_binomial will be NA
  mutate(x_binomial = if_else(is.na(observed_rate) | is.na(count), NA_real_, round(observed_rate * count))) %>%
  rowwise() %>%
  mutate(ci_calcs = list(binomial_ci(x_binomial, count))) %>%
  ungroup() %>%
  tidyr::unnest_wider(ci_calcs, names_sep = "_") %>%
  # Rename to use the calculated CIs
  rename(lower_ci = ci_calcs_lower_ci_binom, upper_ci = ci_calcs_upper_ci_binom) %>%
  # Recalculate lengths for error bars
  mutate(
    y_lower_length = observed_rate - lower_ci,
    y_upper_length = upper_ci - observed_rate
  ) %>%
  # Select relevant columns, ensure correct types for plotting
  select(
    trait, ancestry_group, count, expected_rate, observed_rate, 
    lower_ci, upper_ci, y_lower_length, y_upper_length
  ) %>%
  # Filter out rows where essential plotting data might be NA after calculations
  filter(!is.na(expected_rate) & !is.na(observed_rate) & !is.na(lower_ci) & !is.na(upper_ci))

# --- Main Script ---

# 1. P-value Replication Plots
for (anc in unique(rates$ancestry_group)) {
  rates_anc <- rates %>% 
    filter(ancestry_group == anc) %>%
    group_by(trait) %>%
    slice(1) %>%
    ungroup()
  
  p3 <- plot_with_emphasis(
    data = rates_anc,
    x_var = "expected_rate",
    y_var = "observed_rate",
    y_lower_var = "y_lower_length",
    y_upper_var = "y_upper_length",
    color_var = "trait",
    label_var = "trait",
    emphasize_traits = names(emphasized_traits),
    colors = emphasized_traits,
    title = paste("P-value Replication Rate:", toupper(anc)),
    xlab = "Expected Replication Rate (95% CI)",
    ylab = "Observed Replication Rate (95% CI)",
    make_square_axes = TRUE,
    add_y_equals_x_line = TRUE,
    errorbar_scale = 1 # Since we are providing lengths directly
  )
  
  filename_base <- tools::file_path_sans_ext(basename(args$input_file))
  ggsave(
    filename = file.path(output_dir, paste0(filename_base, "_", anc, ".pdf")),
    plot = p3,
    width = 8,
    height = 8
  )
  cat("Saved plot to:", file.path(output_dir, paste0(filename_base, "_", anc, ".pdf")), "\n")
}

# 3. Meta-analyzed P-value Replication Rates Plot

# Calculate combined counts (sum of rate * count) and total significant genes
meta_rates_data_for_print <- rates %>%
  group_by(ancestry_group) %>%
  summarise(
    expected_sum_prod = sum(expected_rate * count, na.rm = TRUE),
    observed_sum_prod = sum(observed_rate * count, na.rm = TRUE),
    total_significant_genes = sum(count, na.rm = TRUE) # Assuming 'count' is the number of significant genes for that trait
  )

# Print total significant genes per ancestry group
cat("\nTotal significant genes per ancestry group (for meta-analyzed plot caption):\n")
for (i in 1:nrow(meta_rates_data_for_print)) {
  cat(paste0(meta_rates_data_for_print$ancestry_group[i], ": ", meta_rates_data_for_print$total_significant_genes[i], "\n"))
}
cat("\n")

meta_rates_summary <- meta_rates_data_for_print %>%
  select(ancestry_group, expected_sum_prod, observed_sum_prod) %>%
  tidyr::pivot_longer(
    cols = c(expected_sum_prod, observed_sum_prod),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(
    type = factor(type, levels = c("expected_sum_prod", "observed_sum_prod"),
                  labels = c("Expected count", "Observed count")) # Updated labels
  )

# Create the horizontal bar plot
p_meta_rates <- ggplot(meta_rates_summary, aes(x = value, y = ancestry_group, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.2f", value)), 
            position = position_dodge(width = 0.9), 
            hjust = -0.1, size = 3) + # Adjust hjust if labels overlap bars
  scale_fill_manual(values = c("Expected count" = "skyblue", "Observed count" = "salmon")) +
  labs(
    title = NULL, # No title
    x = NULL,     # No x-axis label
    y = NULL,     # No y-axis label
    fill = NULL   # No legend title
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_blank(), # Ensure title element is blank
    legend.position = "top",
    panel.grid.major = element_blank(), # No major grid lines
    panel.grid.minor = element_blank(), # No minor grid lines
    axis.ticks = element_blank(),       # No axis ticks
    axis.text.x = element_blank(),       # Remove x-axis text as counts are next to bars
    axis.text.y = element_text(size=10)  # Keep y-axis text for ancestry groups, or element_blank()
  ) +
  # Ensure x-axis starts at 0 and extends to accommodate labels
  coord_cartesian(xlim = c(0, max(meta_rates_summary$value, na.rm = TRUE) * 1.15)) 

# Save the meta-analyzed plot
ggsave(
  filename = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(args$input_file)), ".pdf")),
  plot = p_meta_rates,
  width = 10,
  height = 6
)
cat("Saved plot to:", file.path(output_dir, paste0(input_filename_base, ".pdf")), "\n")
