#!/usr/bin/env Rscript

# --- Configuration ---
convert_to_per_allele <- TRUE
annotation <- "pLoF"
output_dir <- "tables"
ancestry_groups <- c("afr", "amr", "eur")
# trait_subset <- c("height", "BMI", "low density lipoprotein")
trait_subset <- NULL
significance_label <- "nominal"
pval_threshold <- ifelse(significance_label == "nominal", 0.01, NULL) 

# --- Libraries ---
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# --- Source the replication script to access the function ---
source("luke/replication.R")

# --- Helper function for binomial confidence intervals ---
binomial_ci <- function(x, n, conf.level = 0.95) {
  alpha <- 1 - conf.level
  lower <- qbeta(alpha/2, x, n - x + 1)
  upper <- qbeta(1 - alpha/2, x + 1, n - x)
  data.frame(lower = lower, upper = upper)
}

# --- Create output directory ---
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- Load and prepare phenotype data ---
pheno_info <- read_csv("tables/aou_phenotype_info.csv", col_types = cols())

traits <- pheno_info %>%
  filter(!is.na(ukb_phenoname) & ancestry %in% ancestry_groups & trait_type == "continuous") %>%
  select(name = description, ukbb = ukb_phenoname, aou = phenoname, ancestry) %>%
  distinct()

if (!is.null(trait_subset)) {
  traits <- traits %>%
    filter(name %in% trait_subset)
}

# --- Initialize data frames to store results ---
rates_table <- data.frame()
effects_table <- data.frame()

# --- Process each trait ---
for (trait_info in split(traits, 1:nrow(traits))) {
  cat(paste0("\nProcessing trait: ", trait_info$name, " (", trait_info$ancestry, ")... "))
  
  # Run replication analysis with error handling
  result <- tryCatch({
    perform_replication(
      annotation = annotation,
      ukbb_pheno = trait_info$ukbb,
      aou_pheno = trait_info$aou,
      aou_ancestry = trait_info$ancestry,
      convert_to_per_allele = convert_to_per_allele,
      pval_threshold = pval_threshold
    )
  }, error = function(e) {
    cat("Error processing", trait_info$name, "(", trait_info$ancestry, "):", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Skip if there was an error
  if (is.null(result)) {
    cat("Skipping", trait_info$name, "(", trait_info$ancestry, ") due to error\n")
    next
  }
  
  # Extract p-value summary
  pval_summary <- result$pval_summary
  mean_summary <- result$mean_summary
  
  # Add to rates table
  if (!is.null(pval_summary)) {
    ci <- binomial_ci(
      round(pval_summary$observed_rate * pval_summary$n_genes),
      pval_summary$n_genes
    )
    
    rates_table <- rbind(rates_table, data.frame(
      trait = trait_info$name,
      ancestry_group = trait_info$ancestry,
      count = pval_summary$n_genes,
      expected_rate = pval_summary$expected_rate,
      observed_rate = pval_summary$observed_rate,
      lower_ci = ci$lower,
      upper_ci = ci$upper,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add to effects table
  if (!is.null(mean_summary)) {
    effects_data <- mean_summary %>%
      select(bin, count, expected, expected_se, mean_aou, mean_aou_se, mean_ukbb, mean_ukbb_se) %>%
      mutate(
        trait = trait_info$name,
        ancestry_group = trait_info$ancestry,
        .before = 1
      )
    
    effects_table <- rbind(effects_table, effects_data)
  }
  
  cat("done")
}


# --- Save tables ---
if (nrow(rates_table) > 0) {
  save_name <- ifelse(significance_label == "nominal", "replication_rates_nominal.csv", "replication_rates.csv")
  write.csv(rates_table, file.path(output_dir, save_name), row.names = FALSE)
  cat("\nSaved replication rates to:", file.path(output_dir, save_name), "\n")
}

if (nrow(effects_table) > 0) {
  write.csv(effects_table, file.path(output_dir, "effect_replication.csv"), row.names = FALSE)
  cat("Saved effect replication data to:", file.path(output_dir, "effect_replication.csv"), "\n")
}

if (nrow(rates_table) == 0 && nrow(effects_table) == 0) {
  cat("\nNo data was generated. Please check the input parameters and data.\n")
}
