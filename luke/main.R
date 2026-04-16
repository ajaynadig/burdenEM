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
source("R/EM.R") # Source for EM_fit
source("luke/gene_features.R") # Source for get_bins
source("luke/burdenEM.R") # Source for fit_burdenem_model
source("luke/binary_trait_likelihood.R") 


#' Run BurdenEM RVAS Workflow
#'
#' Loads variant data, processes it to gene-level, initializes the BurdenEM model,
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
#' @param removed_variants_file Optional: file(s) listing variants to remove for LD pruning
#'   (binary traits only). Each file must have a \code{locus} (chr:pos) or \code{SNP}
#'   (chr:pos:ref:alt) column. Two forms:
#'   \itemize{
#'     \item Single path (character): applied to all variants regardless of cohort.
#'     \item Named list: per-cohort pruning for meta-analysis. Keys must match the
#'       \code{dataset} column in the variant data (e.g.,
#'       \code{list(aou_afr = "afr_removed.tsv", aou_eur = "eur_removed.tsv")}).
#'   }
#'   Default: NULL.
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
    burdenem_grid_size = 4,
    num_iter = 5000,
    output_file_prefix = NULL, # CLI provides the full path prefix for output files
    verbose = FALSE,
    per_allele_effects = FALSE,
    correct_for_ld = FALSE,
    removed_variants_file = NULL,
    binary_trait_model_type = "betabinom",
    optimizer = "EM",
    variant_data = NULL
) {

    # --- 1. Load Variant-Level Data ---
    if(verbose) message("\n--- Loading Variant-Level Data ---")
    if (!is.null(variant_data)) {
        # Use pre-loaded variant data (e.g., on-the-fly meta construction)
        variant_data <- variant_data %>% distinct()
        if(verbose) message(paste("Using pre-loaded variant data:", nrow(variant_data), "variants."))
    } else {
        variant_data <- load_variant_files_with_category(
            variant_dir = variant_dir,
            variant_file_pattern = variant_file_pattern, # Use pattern from CLI
            # data_name and pheno arguments removed
            annotations_to_process = c(annotation_to_process),
            frequency_range = frequency_range
        ) %>% distinct()
        if(verbose) message(paste("Successfully loaded", nrow(variant_data), "variants."))
    }

    if (nrow(variant_data) == 0) {
        warning(paste("No variant data available. Skipping this study."))
        return(NULL)
    }

    # --- Detect Trait Type ---
    required_binary_cols <- c("AC_cases", "N")
    is_binary <- all(required_binary_cols %in% names(variant_data))
    trait_type <- if(is_binary) "binary" else "continuous"
    message(paste("Detected trait type:", trait_type))

    # --- 2. LD-Prune Variants for Binary Traits ---
    # removed_variants_file: single path (applied to all variants) or named list
    # mapping dataset -> path (per-cohort pruning for meta-analysis)
    if (!is.null(removed_variants_file) && is_binary) {
        if(verbose) message("\n--- LD-Pruning Variants for Binary Trait ---")

        # Helper: read a removed-variants file and return locus vector
        .read_removed_loci <- function(filepath) {
            df <- data.table::fread(filepath)
            if ("locus" %in% names(df)) return(df$locus)
            if ("SNP" %in% names(df)) return(sub("^([^:]+:[^:]+):.*", "\\1", df$SNP))
            stop(paste("Cannot identify locus column in:", filepath, "— expected 'locus' or 'SNP'."))
        }

        # Construct .locus column for matching
        has_chr_pos <- all(c("CHR", "POS") %in% names(variant_data))
        has_locus <- "locus" %in% names(variant_data)
        if (!has_chr_pos && !has_locus) {
            stop("Cannot construct locus for variant filtering. Need CHR+POS or locus columns.")
        }
        variant_data <- variant_data %>%
            dplyr::mutate(.locus = if (has_chr_pos) paste0(CHR, ":", POS) else locus)

        original_count <- nrow(variant_data)

        if (is.list(removed_variants_file)) {
            # Named list: per-cohort pruning (e.g., list(aou_afr="path1", aou_eur="path2"))
            if (!("dataset" %in% names(variant_data))) {
                stop("removed_variants_file is a named list but variant_data has no 'dataset' column.")
            }
            for (ds in names(removed_variants_file)) {
                removed_loci <- .read_removed_loci(removed_variants_file[[ds]])
                variant_data <- variant_data %>%
                    dplyr::filter(!(dataset == ds & .locus %in% removed_loci))
            }
        } else {
            # Single file: applied to all variants
            removed_loci <- .read_removed_loci(removed_variants_file)
            variant_data <- variant_data %>%
                dplyr::filter(!(.locus %in% removed_loci))
        }

        variant_data <- variant_data %>% dplyr::select(-.locus)
        message(paste("LD pruning: removed", original_count - nrow(variant_data), "of", original_count, "variants"))
    }

    # --- 3. Load LD-Corrected Burden Scores (for features and continuous trait LD correction) ---
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

    # --- 4. Process Variants to Gene-Level ---
    if (verbose) message(paste("\n--- Processing Variants to Gene-Level (", trait_type, ") ---"))
    if (trait_type == "binary") {
        if (binary_trait_model_type == "betabinom") {
        if(verbose) message("Estimating overdispersion for binary trait using beta-binomial model...")
        variant_data <- estimate_overdispersion_binomial(variant_data, intercept_frequency_bin_edges = intercept_frequency_bin_edges, verbose = verbose)
        } else if (binary_trait_model_type == "nbinom") {
            if(verbose) message("Estimating overdispersion for binary trait using negative binomial model...")
            variant_data <- estimate_overdispersion_poisson(variant_data, intercept_frequency_bin_edges = intercept_frequency_bin_edges, verbose = verbose)    
        } else {
            variant_data$overdispersion <- 0
            if (verbose) message("Not estimating overdispersion for binary trait.")
        }
        
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

        moment_h2_est = sum(gene_level_data$gamma_per_sd^2 - gene_level_data$gene_intercept) 
        message(paste("Moment-based H2 estimate:", moment_h2_est))
    }

    if (nrow(gene_level_data) == 0) {
        stop("No gene-level data remained after processing and filtering. Check input data and filters.")
    }
    message(paste("Aggregated to", nrow(gene_level_data), "gene-category rows."))

    # --- 5. Append LD-Corrected Scores ---
    if(!is.null(ld_corrected_scores_file)){ # Do not change this to 'if correct_for_ld'
      if(verbose) message("\n--- Appending LD-Corrected Burden Scores ---")
      join_on_column_name <- "gene"
      if (grepl("^ENSG", gene_level_data$gene[1])) {
          gene_level_data <- gene_level_data %>% dplyr::rename(gene_id = gene)
          join_on_column_name <- "gene_id"
      }
      gene_level_data <- gene_level_data %>%
          left_join(ld_corrected_scores_df, by = join_on_column_name)

      # For binary traits: burden_score is computed in-pipeline from (pruned) variants
      # LD-corrected scores file is only used for feature columns (e.g., lof.oe)
      if (trait_type == "binary") {
          if(verbose) message("Binary trait: using in-pipeline burden_score (LD-corrected scores used for features only)")
      }
    }

    # --- 6. Specify Effect Sizes (Continuous Only) & Prepare for Model ---
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

    # --- 7. Prepare Gene Features & Add to Gene Data ---
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

    # --- 8. Validate Gene-Level Data for Binary Traits ---
    if (trait_type == "binary") {
        if(verbose) message("\n--- Validating Binary Trait Data ---")
        # Check that required list columns exist and are valid
        required_list_cols <- c("AC_cases", "expected_count", "overdispersion")
        for (col in required_list_cols) {
            if (!(col %in% names(gene_level_data))) {
                stop(paste("Required column", col, "not found in gene_level_data"))
            }
            # Check that it's actually a list column
            if (!is.list(gene_level_data[[col]])) {
                stop(paste("Column", col, "is not a list column as expected"))
            }
            # Check for NULL or empty lists
            n_invalid <- sum(sapply(gene_level_data[[col]], function(x) is.null(x) || length(x) == 0))
            if (n_invalid > 0) {
                warning(paste("Found", n_invalid, "genes with empty", col, "list. Filtering them out."))
                gene_level_data <- gene_level_data %>%
                    dplyr::filter(sapply(!!rlang::sym(col), function(x) !is.null(x) && length(x) > 0))
            }
        }
        # Check for required list columns (per-variant N and prevalence for meta-analysis h2)
        for (col in c("N", "prevalence")) {
            if (!(col %in% names(gene_level_data))) {
                stop(paste("Required column", col, "not found in gene_level_data"))
            }
            if (!is.list(gene_level_data[[col]])) {
                stop(paste("Column", col, "should be a list column"))
            }
        }
        # Check for required scalar columns
        for (col in c("burden_score")) {
            if (!(col %in% names(gene_level_data))) {
                stop(paste("Required column", col, "not found in gene_level_data"))
            }
            n_na <- sum(is.na(gene_level_data[[col]]))
            if (n_na > 0) {
                warning(paste("Found", n_na, "genes with NA", col, ". Filtering them out."))
                gene_level_data <- gene_level_data %>%
                    dplyr::filter(!is.na(!!rlang::sym(col)))
            }
        }
        if (nrow(gene_level_data) == 0) {
            stop("No valid genes remaining after validation")
        }
        message(paste("Validation passed:", nrow(gene_level_data), "genes ready for fitting"))
    }

    # --- 9. Define Likelihood Function & Fit Model ---
    if(verbose) message("\n--- Defining Likelihood Function & Fitting Model ---")

    # Define likelihood_fn and h2_fn based on trait type
    to_per_allele_effects <- if (per_allele_effects) {
        function(row, beta) beta
    } else {
        function(row, beta) {
            if (row[["burden_score"]] > 0) {
                beta / sqrt(row[["burden_score"]])
            } else {
                rep(0, length(beta))
            }
        }
    }

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
        likelihood_function <- get_gene_likelihood(model_type = binary_trait_model_type, to_per_allele_effects = to_per_allele_effects)
        pval_function <- function(row) ppois(row[["CAC_cases"]], lambda = row[["expected_CAC_cases"]], lower.tail = FALSE)

        get_power_function <- function(pval_threshold, samplesize_ratio){
            power_function <- function(beta, row) {
                beta_per_allele <- to_per_allele_effects(row, beta)
                mu <- samplesize_ratio * row[["expected_CAC_cases"]]
                AC_critical_upper <- qpois(pval_threshold, lambda = mu, lower.tail = FALSE)
                AC_critical_lower <- qpois(pval_threshold, lambda = mu, lower.tail = TRUE)
                result <- beta
                result[beta == 0] <- 0
                result[beta > 0] <- ppois(AC_critical_upper, lambda = mu * exp(beta_per_allele[beta > 0]), lower.tail = FALSE)
                result[beta < 0] <- ppois(AC_critical_lower, lambda = mu * exp(beta_per_allele[beta < 0]), lower.tail = TRUE)
                return(result)
            }
            return(power_function)
        }

        rel_samplesize_function <- function(row) row[["expected_CAC_cases"]]

        current_h2_function <- function(beta, row) {
            beta_per_allele <- to_per_allele_effects(row, beta)
            rate_ratio <- pmin(exp(beta_per_allele), 1/row[["mean_prevalence"]])
            row[["h2_numer_scale"]] * (rate_ratio - 1)^2 / row[["h2_denom"]]
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
        verbose = verbose,
        optimizer = optimizer
    )
    burdenem_model$pval_function <- pval_function
    burdenem_model$get_power_function <- get_power_function
    burdenem_model$rel_samplesize_function <- rel_samplesize_function
    burdenem_model$to_per_allele_effects <- to_per_allele_effects

    if(verbose) print(burdenem_model$delta)

    # --- 10. Save Model & Data ---
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