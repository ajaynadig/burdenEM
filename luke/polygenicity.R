#!/usr/bin/env Rscript

# --- Libraries ---
# Install future and furrr if not already installed
if (!requireNamespace("future", quietly = TRUE)) {
  message("Installing 'future' package...")
  install.packages("future", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("furrr", quietly = TRUE)) {
  message("Installing 'furrr' package...")
  install.packages("furrr", repos = "http://cran.us.r-project.org")
}
library(future)
library(furrr)
library(dplyr)
library(readr)
library(stringr)
library(purrr) # For bind_rows

# --- Source R scripts ---
# Note: model.R is sourced within process_single_study_for_polygenicity for furrr compatibility

# --- Helper function to process a single study for polygenicity --- 
process_single_study_for_polygenicity <- function(study_row, current_annotation, current_verbose) {
  source("R/model.R") # Source here for furrr workers
  if (current_verbose) {
    message(sprintf("Processing polygenicity for study: %s, dataset: %s, trait_type: %s, annotation: %s", 
                    study_row$abbreviation, study_row$dataset, study_row$trait_type, current_annotation))
  }

  model_filename_pattern <- study_row$model_filename
  actual_model_filename <- stringr::str_replace(model_filename_pattern, "<ANNOTATION>", current_annotation)

  model_fit <- tryCatch({
    readRDS(actual_model_filename)
  }, error = function(e) {
    if (current_verbose) message(sprintf("  Error loading model '%s': %s", actual_model_filename, e$message))
    return(NULL)
  })

  metric_values_list <- list()
  if (is.null(model_fit)) {
    for(metric_name in names(polygenicity_functions)) {
      metric_values_list[[metric_name]] <- NA_real_
      metric_values_list[[paste0(metric_name, "_se")]] <- NA_real_
    }
  } else {
    pi_i <- model_fit$pi_i 
    sigma_sq_i <- model_fit$sigma_sq_i

    for (metric_name in names(polygenicity_functions)) {
      polygenicity_fn <- polygenicity_functions[[metric_name]]
      calculated_metric <- polygenicity_fn(model_fit) 
      metric_values_list[[metric_name]] <- calculated_metric$mean 
      metric_values_list[[paste0(metric_name, "_se")]] <- calculated_metric$se
    }
  }
  
  # Combine selected study identifiers with the list of metric values
  selected_study_info <- as.data.frame(study_row)[, c("abbreviation", "dataset"), drop = FALSE]
  current_result_row <- dplyr::bind_cols(selected_study_info, as.data.frame(metric_values_list))
  return(current_result_row)
}

# --- Define polygenicity functions ---
effective_polygenicity <- function(model) {
  h2 <- sum(posterior_expectation2(model, model$h2_function))
  integrand <- function(x,row) model$h2_function(x,row)^2 / h2
  inside_term <- posterior_expectation_with_se(model, integrand)
  mean <- h2 / sum(inside_term$mean)
  se <- abs(mean) * sqrt(sum(inside_term$se^2)) / sum(inside_term$mean)
  return(list(mean = mean, se = se))
}


effective_mutational_polygenicity <- function(model) {

  # Mutational variance, as opposed to heritability, contributed by each gene
  mutational_h2_function <- function(x,row) {
    if (!"lof.mu" %in% names(row)) return(rep(NA, length(x)))
    result <- model$h2_function(x,row) * row$lof.mu / row$burden_score
    if(is.na(row$lof.mu)) 0*x else result
  }

  h2 <- sum(posterior_expectation2(model, mutational_h2_function))
  integrand <- function(x,row) mutational_h2_function(x,row)^2 / h2
  inside_term <- posterior_expectation_with_se(model, integrand)
  mean <- h2 / sum(inside_term$mean)
  se <- abs(mean) * sqrt(sum(inside_term$se^2)) / sum(inside_term$mean)
  return(list(mean = mean, se = se))
}

effective_effect_var_polygenicity <- function(model) {

  # Mutational variance, as opposed to heritability, contributed by each gene
  effect_var_function <- function(x,row) {
    model$h2_function(x,row) / row$burden_score
  }

  h2 <- sum(posterior_expectation2(model, effect_var_function))
  integrand <- function(x,row) effect_var_function(x,row)^2 / h2
  inside_term <- posterior_expectation_with_se(model, integrand)
  mean <- h2 / sum(inside_term$mean)
  se <- abs(mean) * sqrt(sum(inside_term$se^2)) / sum(inside_term$mean)
  return(list(mean = mean, se = se))
}

entropy_polygenicity <- function(model) {
  h2 <- sum(posterior_expectation2(model, model$h2_function))
  integrand <- function(x, row) {
    h2_gene <- model$h2_function(x, row)
    ifelse(h2_gene == 0, 0, -h2_gene * log(h2_gene) / h2)
  }
  inside_term <- posterior_expectation_with_se(model, integrand)
  mean <- h2 * exp(sum(inside_term$mean))
  se <- h2 * exp(sum(inside_term$mean)) * sqrt(sum(inside_term$se^2))
  return(list(mean = mean, se = se))
}

# Avoid exponent overflow; assumes that a gene cannot have heritability
# greater than exp(MAX_EXPONENT) * offset, where offset is defined as
# the maximum expected heritability of any gene.

softmax_polygenicity <- function(model) {
  MAX_EXPONENT <- 10 
  gene_h2 <- posterior_expectation2(model, model$h2_function)
  h2 <- sum(gene_h2)
  offset <- max(gene_h2)
  integrand <- function(x, row) {
    p_gene <- model$h2_function(x, row)
    exponent <- pmin(1/offset-1/p_gene, MAX_EXPONENT)
    ifelse(p_gene == 0, 0, p_gene * exp(exponent) / h2)
  }
  inside_term <- posterior_expectation_with_se(model, integrand)
  mean <- h2 * (1/offset - log(sum(inside_term$mean)))
  se <- h2 * sqrt(sum(inside_term$se^2)) / sum(inside_term$mean)
  
  return(list(mean = mean, se = se))
}

polygenicity_functions <- list(
  softmax_polygenicity = softmax_polygenicity,
  effective_polygenicity = effective_polygenicity,
  entropy_polygenicity = entropy_polygenicity,
  effective_mutational_polygenicity = effective_mutational_polygenicity,
  effective_effect_var_polygenicity = effective_effect_var_polygenicity
)

# --- Main function to be called by CLI --- 
calculate_polygenicity_metrics <- function(studies_df, annotation, verbose = FALSE, run_sequentially_param = FALSE) {
  if (!is.data.frame(studies_df) || nrow(studies_df) == 0) {
    if (verbose) message("Input studies_df is empty or not a dataframe. Returning empty results.")
    # Return an empty df with expected column structure if studies_df is empty
    # Ensure only abbreviation and dataset are kept from the original studies_df structure
    empty_df_template <- studies_df[0, c("abbreviation", "dataset"), drop = FALSE]
    # Create template for metric names including SE columns
    metric_col_names <- character(0)
    for (metric_name in names(polygenicity_functions)) {
      metric_col_names <- c(metric_col_names, metric_name, paste0(metric_name, "_se"))
    }
    metric_names_template <- stats::setNames(data.frame(matrix(ncol = length(metric_col_names), nrow = 0)), metric_col_names)
    return(dplyr::bind_cols(empty_df_template, metric_names_template))
  }

  if (!run_sequentially_param) {
    if (verbose) {
      message(sprintf("Setting up parallel execution for polygenicity calculation using up to %d workers.", future::availableCores()))
    }
    future::plan(future::multisession)
    final_results_df <- furrr::future_map_dfr(
      .x = 1:nrow(studies_df),
      .f = ~process_single_study_for_polygenicity(studies_df[.x, , drop = FALSE], annotation, verbose),
      .progress = verbose
    )
  } else {
    if (verbose) {
      message("Running polygenicity calculation in sequential mode.")
    }
    final_results_df <- purrr::map_dfr(
      .x = 1:nrow(studies_df),
      .f = ~process_single_study_for_polygenicity(studies_df[.x, , drop = FALSE], annotation, verbose)
    )
  }
  
  return(final_results_df)
}

# --- Helper function for standalone processing of a single trait ---
process_trait_for_polygenicity_metrics_standalone <- function(trait_name, phenotype_code, current_annotation_name, verbose_standalone = TRUE) {
  if (verbose_standalone) {
    cat(paste0("Processing Trait: ", trait_name, " (", phenotype_code, ") for annotation: ", current_annotation_name, "\n"))
  }

  model_file <- file.path("fitted_models",
                          paste0("burdenEM_fit_genebass_", phenotype_code, "_", current_annotation_name, ".rds"))

  if (!file.exists(model_file)) {
    if (verbose_standalone) {
      cat(paste0("  Model file not found: ", model_file, ". Skipping.\n"))
    }
    return(NULL) # Return NULL if model file not found
  }
  
  model <- tryCatch({
      readRDS(model_file)
  }, error = function(e) {
      if (verbose_standalone) {
          cat(paste0("  Error loading model file ", model_file, ": ", e$message, "\n"))
      }
      return(NULL) # Return NULL if model cannot be loaded
  })

  if (is.null(model)) {
      return(NULL)
  }
  
  current_trait_results <- list(
    trait = trait_name,
    phenotype_code = phenotype_code
    # annotation = current_annotation_name # Optional: add annotation to results
  )

  for (metric_name in names(polygenicity_functions)) {
    metric_calculator_func <- polygenicity_functions[[metric_name]]
    if (verbose_standalone) {
      cat(paste0("  Calculating metric: ", metric_name, "\n"))
    }

    metric_output <- tryCatch({
      metric_calculator_func(model) # Expected to return list(mean = SCALAR, se = SCALAR)
    }, error = function(e) {
      if (verbose_standalone) {
        cat(paste0("    ERROR calculating ", metric_name, " for ", trait_name, ": ", e$message, "\n"))
      }
      list(mean = NA_real_, se = NA_real_) # Return NA on error
    })

    current_trait_results[[paste0(metric_name, "_mean")]] <- metric_output$mean
    current_trait_results[[paste0(metric_name, "_se")]] <- metric_output$se
  }
  return(current_trait_results)
}

# --- Standalone script execution logic (if file is run directly) ---
if (sys.nframe() == 0) { # Check if the script is being run directly (not sourced)
  cat("Running polygenicity.R as a standalone script...\n")
  
  # --- Load and prepare phenotype data ---
  pheno_info_raw <- readr::read_csv("tables/aou_phenotype_info.csv", 
    col_types = readr::cols(),  
    show_col_types = FALSE)
  
  traits_df <- pheno_info_raw %>%
    dplyr::filter(trait_type == "continuous", !is.na(ukb_phenoname)) %>%
    dplyr::select(name = description, ukbb_pheno = ukb_phenoname) %>%
    dplyr::distinct()
  
  # --- Common parameters ---
  annotation_name <- "pLoF" # Default annotation for standalone run
  output_file <- paste0("tables/polygenicity_metrics_standalone_", annotation_name, ".csv")
  
  if (nrow(traits_df) > 0) {
      all_trait_metrics_list <- purrr::map(1:nrow(traits_df), function(i) {
          process_trait_for_polygenicity_metrics_standalone(
              trait_name = traits_df$name[i],
              phenotype_code = traits_df$ukbb_pheno[i],
              current_annotation_name = annotation_name,
              verbose_standalone = TRUE 
          )
      })
  
      # Filter out NULL results (e.g., from missing model files or load errors)
      all_trait_metrics_list <- purrr::compact(all_trait_metrics_list)
  
      if (length(all_trait_metrics_list) > 0) {
        final_df <- dplyr::bind_rows(all_trait_metrics_list)
  
        # Define desired column order
        fixed_cols <- c("trait", "phenotype_code")
        metric_cols_mean <- paste0(names(polygenicity_functions), "_mean")
        metric_cols_se <- paste0(names(polygenicity_functions), "_se")
        
        all_expected_cols <- c(fixed_cols, metric_cols_mean, metric_cols_se)
        for(col_name in all_expected_cols){
            if(!col_name %in% names(final_df)){
                final_df[[col_name]] <- NA_real_ # Ensure numeric NA
            }
        }
        
        final_df <- final_df[, all_expected_cols] 
  
        readr::write_csv(final_df, output_file)
        cat(paste0("Standalone polygenicity metrics saved to: ", output_file, "\n"))
      } else {
        cat("No polygenicity metrics could be calculated for any traits in standalone mode.\n")
      }
  } else {
      cat("No traits selected for polygenicity calculation in standalone mode.\n")
  }
  
  cat("Standalone script part of polygenicity.R finished.\n")
}

