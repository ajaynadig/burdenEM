# main.R
# Contains core functions for running burdenEM
# For command-line interface, use main_cli.R instead

# Paths used for default values
TOP_LEVEL_DATA_DIR <- "../../burdenEM_results/"
OUTPUT_DIR <- "./burdenEM_results"

# --- Libraries ---
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
options(vroom.show_problems = FALSE)

# --- Source Required Functions ---
source("R/io.R")
source("luke/intercept.R")
source("luke/variant_to_gene.R")
source("R/model.R")
source("R/likelihoods.R")
source("R/EM.R") # Source for EM_fit, bootstrap_EM, null_EM_rvas
source("luke/gene_features.R") # Source for get_bins
source("luke/burdenEM.R") # Source for fit_burdenem_model
source("luke/negative_binomial.R") # Source for estimate_overdispersion, aggregate_variants_to_gene_lists, gene_likelihood


#' Run BurdenEM RVAS Workflow
#'
#' Loads variant data, processes it to gene-level, initializes the BurdenEM model,
#' fits the model using EM, and optionally performs bootstrapping and null simulations.
#'
#' @param variant_dir Directory containing variant files.
#' @param data_name Dataset name prefix (e.g., 'genebass').
#' @param pheno Phenotype name (e.g., '50_NA').
#' @param annotation_to_process Character string for the single functional category to load.
#' @param intercept_frequency_bin_edges Numeric vector defining AF bin edges, only for intercept calculation.
#' @param feature_col_name Optional string: name of the column in the LD-corrected scores file to use for gene features (e.g., "lof_oe"). Default: NULL.
#' @param num_feature_bins Integer: number of quantile bins for the feature column. Default: 5.
#' @param num_positive_components: Number of positive components for BurdenEM; total number is 2*num_positive_components + 1.
#' @param burdenem_grid_size Grid size for likelihood calculation.
#' @param num_iter Number of iterations for EM algorithm.
#' @param verbose Logical, whether to print detailed messages.
#' @param frequency_range Numeric vector c(min, max) for AF range filter [min,max). Default: c(0, 0.001).
#' @param per_allele_effects Logical, whether to calculate per-allele effect sizes instead of per-gene. Default: FALSE.
#' @param correct_for_ld Logical, whether to apply LD correction to burden scores and gamma. Default: FALSE.
#' @param correct_genomewide Logical, whether to apply genome-wide burden correction. Default: FALSE.
#' @return A list object representing the fitted BurdenEM model.
#' @export
run_burdenEM_rvas <- function(
    variant_dir = file.path(TOP_LEVEL_DATA_DIR, "data", "genebass", "var_txt"),
    ld_corrected_scores_file = NULL,
    variant_file_pattern = NULL, # CLI provides this pattern for sumstats
    # data_name and pheno removed, now part of variant_file_pattern
    annotation_to_process='pLoF',
    intercept_frequency_bin_edges = c(0, 1e-5,1e-4, 1e-3),
    frequency_range = c(0, 0.001),
    feature_col_name = "lof.oe",
    num_feature_bins = 5,
    num_positive_components = 10,
    burdenem_grid_size = 100,
    num_iter = 10000,
    bootstrap = TRUE,
    output_file_prefix = NULL, # CLI provides the full path prefix for output files
    verbose = FALSE,
    per_allele_effects = FALSE,
    customize_components = FALSE,
    correct_for_ld = TRUE,
    correct_genomewide = FALSE
) {

    # --- 1. Load Variant-Level Data ---
    if(verbose) message("\n--- Loading Variant-Level Data ---")
    variant_data <- load_variant_files_with_category(
        variant_dir = variant_dir,
        variant_file_pattern = variant_file_pattern, # Use pattern from CLI
        # data_name and pheno arguments removed
        annotations_to_process = c(annotation_to_process),
        frequency_range = frequency_range
    )
    if(verbose) message(paste("Successfully loaded", nrow(variant_data), "variants."))

    # --- Detect Trait Type ---
    required_binary_cols <- c("AC_cases", "N")
    is_binary <- all(required_binary_cols %in% names(variant_data))
    trait_type <- if(is_binary) "binary" else "continuous"
    message(paste("Detected trait type:", trait_type))


    # --- 2. Load LD-Corrected Burden Scores ---
    if(!is.null(ld_corrected_scores_file)){ # Do not change this to 'if correct_for_ld'
      if(verbose) message("\n--- Loading LD-Corrected Burden Scores ---")
      lower_fmt <- format(frequency_range[1], nsmall=1, scientific=FALSE)
      upper_fmt <- format(frequency_range[2], nsmall=1, scientific=FALSE)
      dynamic_ld_scores_file <- ld_corrected_scores_file %>%
          gsub("<ANNOTATION>", annotation_to_process, ., fixed = TRUE) %>%
          gsub("<LOWER>", lower_fmt, ., fixed = TRUE) %>%
          gsub("<UPPER>", upper_fmt, ., fixed = TRUE)
      ld_corrected_scores_df <- data.table::fread(dynamic_ld_scores_file) %>%
          dplyr::mutate(gene = as.character(gene)) %>%
          dplyr::rename(burden_score_ld = burden_score)
      if(verbose) message(paste("Loaded", nrow(ld_corrected_scores_df),
          "LD-corrected scores from:", dynamic_ld_scores_file))
    }

    # --- 3. Process Variants to Gene-Level ---
    if (verbose) message(paste("\n--- Processing Variants to Gene-Level (", trait_type, ") ---"))
    if (trait_type == "binary") {
        if(verbose) message("Estimating overdispersion for binary trait...")
        variant_data <- estimate_overdispersion(variant_data, intercept_frequency_bin_edges = intercept_frequency_bin_edges, verbose = verbose)
        if(verbose) message("Aggregating variants to gene lists for binary trait...")
        gene_level_data <- aggregate_variants_to_gene_lists(variant_data, verbose = verbose)
        # Join feature_col_name from ld_corrected_scores_df if available for binary traits
        if (!is.null(ld_corrected_scores_df) && !is.null(feature_col_name) && feature_col_name %in% names(ld_corrected_scores_df)) {
            if (verbose) message(paste("Joining feature column '", feature_col_name, "' to gene_level_data for binary trait.", sep=""))
            # Ensure gene column exists in both for merging
            if (!("gene" %in% names(gene_level_data)) || !("gene" %in% names(ld_corrected_scores_df))){
                stop("Gene column missing in gene_level_data or ld_corrected_scores_df, cannot join features.")
            }
        } else if (!is.null(feature_col_name)) {
            if (verbose) warning(paste("Feature column '", feature_col_name, "' not found in ld_corrected_scores_df or ld_corrected_scores_df is NULL. Cannot join for binary trait.", sep=""))
        }
    } else { # Continuous trait
        gene_level_data <- process_variant_to_gene(
            variant_data = variant_data,
            frequency_bin_edges = intercept_frequency_bin_edges
        ) %>%
            dplyr::filter(!is.na(gene))
    }

    if (nrow(gene_level_data) == 0) {
        stop("No gene-level data remained after processing and filtering. Check input data and filters.")
    }
    message(paste("Aggregated to", nrow(gene_level_data), "gene-category rows."))

    # --- 4. Append LD-Corrected Scores ---
    if(!is.null(ld_corrected_scores_file)){ # Do not change this to 'if correct_for_ld'
      if(verbose) message("\n--- Appending LD-Corrected Burden Scores ---")
      join_on_column_name <- "gene"
      if (grepl("^ENSG", gene_level_data$gene[1])) {
          gene_level_data <- gene_level_data %>% dplyr::rename(gene_id = gene)
          join_on_column_name <- "gene_id"
      }
      gene_level_data <- gene_level_data %>%
          left_join(ld_corrected_scores_df, by = join_on_column_name)
    }

    # --- 5. Specify Effect Sizes (Continuous Only) & Prepare for Model ---
    if(trait_type == "continuous"){
        if(verbose) message("\n--- Specifying Gene Effect Sizes (Continuous) ---")
        # Pass per_allele_effects and correct_for_ld from function arguments
        message(paste("Applying per-allele effects:", per_allele_effects))
        gene_level_data <- specify_gene_effect_sizes(gene_level_data,
                                                   per_allele_effects = per_allele_effects,
                                                   correct_for_ld = correct_for_ld)
    } else { # Binary trait
        if(verbose) message("\n--- Skipping Gene Effect Size Specification (Binary traits use list-columns directly) ---")
    }

    # --- 6. Prepare Gene Features & Add to Gene Data ---
    # This step uses the output 'gene_level_data' from step 3/4, independent of trait type
    # It adds the 'features' column based on 'feature_col_name'
    if (!is.null(feature_col_name) && num_feature_bins > 1) {
        if (verbose) message(paste("\n--- Generating", num_feature_bins, "feature bins for:", feature_col_name, "---"))
        if (!(feature_col_name %in% names(gene_level_data))){
             stop(paste("Specified feature_col_name '", feature_col_name, "' not found in gene_level_data after aggregation.", sep=""))
        }
        # Ensure get_bins is available (e.g., source("luke/gene_features.R"))
        if (!exists("get_bins")) stop("Function 'get_bins' not found. Source luke/gene_features.R")
        gene_level_data <- get_bins(gene_level_data, column_name = feature_col_name, num_bins = num_feature_bins, verbose = verbose)
    } else {
        # Assign a single feature bin if no feature column is specified
        gene_level_data$features <- replicate(nrow(gene_level_data), c(1), simplify = FALSE)
    }

    # --- 7. Define Likelihood Function & Fit Model ---
    if(verbose) message("\n--- Defining Likelihood Function & Fitting Model ---")

    # Define likelihood_fn and h2_fn based on trait type
    if(trait_type == "continuous") {
        message("Using Normal likelihood for continuous trait.")
        likelihood_function <- function(row, beta_vec) {
            dnorm(row$effect_estimate, mean = beta_vec, sd = row$effect_se)
        }
        
        pval_function <- function(row) pnorm(row$effect_estimate / row$effect_se, lower.tail = FALSE)

        get_power_function <- function(pval_threshold, samplesize_ratio){
            power_function <- function(beta, row) {
                z_critical <- qnorm(pval_threshold, lower.tail = FALSE)
                mu <- abs(beta * sqrt(samplesize_ratio)) / (row$effect_se)
                pnorm(z_critical, mean = mu, lower.tail = FALSE)
            }
            return(power_function)
        }

        rel_samplesize_function <- function(row) 1 / row$effect_se^2
        
        if (per_allele_effects) {
            current_h2_function <- function(beta, row) row$burden_score * beta^2
        } else {
            current_h2_function <- function(beta, row) beta^2
        }
        current_drop_columns <- NULL
    } else { # Binary trait
        per_allele_effects <- TRUE # TODO

        message("Using Negative Binomial likelihood for binary trait.")
        likelihood_function <- gene_likelihood # From luke/negative_binomial.R
        
        pval_function <- function(row) ppois(row$CAC_cases, lambda = row$expected_CAC_cases, lower.tail = FALSE)
        
        get_power_function <- function(pval_threshold, samplesize_ratio){
            power_function <- function(beta, row) {
                mu <- samplesize_ratio * row$expected_CAC_cases
                AC_critical_upper <- qpois(pval_threshold, lambda = mu, lower.tail = FALSE)
                AC_critical_lower <- qpois(pval_threshold, lambda = mu, lower.tail = TRUE)
                result <- beta
                result[beta == 0] <- 0
                result[beta > 0] <- ppois(AC_critical_upper, lambda = mu * exp(beta[beta > 0]), lower.tail = FALSE)
                result[beta < 0] <- ppois(AC_critical_lower, lambda = mu * exp(beta[beta < 0]), lower.tail = TRUE)
                return(result)
            }
            return(power_function)
        }

        rel_samplesize_function <- function(row) row$expected_CAC_cases

        current_h2_function <- function(beta, row) {
            row$prevalence / (1-row$prevalence) * row$burden_score * (exp(beta) - 1)^2 
        }
        current_drop_columns <- c("AC_cases", "expected_count", "overdispersion")
    }

    # Fit the BurdenEM model using the generalized function
    burdenem_model <- fit_burdenem_model(
        gene_level_data = gene_level_data,
        likelihood_function = likelihood_function,
        h2_function = current_h2_function,
        burdenem_grid_size = burdenem_grid_size,
        num_positive_components = num_positive_components,
        num_iter = num_iter,
        per_allele_effects = per_allele_effects, 
        drop_columns = current_drop_columns, 
        verbose = verbose
    )
    burdenem_model$pval_function <- pval_function
    burdenem_model$get_power_function <- get_power_function
    burdenem_model$rel_samplesize_function <- rel_samplesize_function

    if(verbose) print(burdenem_model$delta)

    # --- 8. Save Model & Data ---
    if(verbose) message("\n--- Saving Model & Data ---")
    # output_prefix construction removed (info now in output_file_prefix from CLI)
    # Directory creation is now handled by the CLI before calling this function

    if (is.null(output_file_prefix)) {
        stop("output_file_prefix is NULL. Cannot save model or data.")
    }

    model_output_path <- paste0(output_file_prefix, ".rds")
    saveRDS(burdenem_model, file = model_output_path)
    if(verbose) message("Model saved successfully.")
    return(burdenem_model)
}