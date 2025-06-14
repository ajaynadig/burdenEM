# --- Define polygenicity functions ---
effective_polygenicity <- function(genes_df) {
  h2_gene <- genes_df$gene_h2
  h2_gene <- h2_gene[!is.na(h2_gene) & h2_gene != 0]
  if (length(h2_gene) == 0) {
    return(0)
  }
  poly <- sum(h2_gene)^2 / sum(h2_gene^2)
}

effective_effect_var_polygenicity <- function(genes_df) {
  h2_gene <- genes_df$effectSize^2
  h2_gene <- h2_gene[!is.na(h2_gene) & h2_gene != 0]
  if (length(h2_gene) == 0) {
    return(0)
  }
  poly <- sum(h2_gene)^2 / sum(h2_gene^2)
}

entropy_polygenicity <- function(genes_df) {
  h2_gene <- genes_df$gene_h2
  h2_gene <- h2_gene[!is.na(h2_gene) & h2_gene != 0]
  if (length(h2_gene) == 0) {
    return(0)
  }

  h2 <- sum(h2_gene)
  inside_term <- 0-h2_gene * log(h2_gene) / h2
  poly <- h2 * exp(sum(inside_term, na.rm=TRUE))
}

# Avoid exponent overflow; assumes that a gene cannot have heritability
# greater than exp(MAX_EXPONENT) * offset, where offset is defined as
# the maximum expected heritability of any gene.

softmax_polygenicity <- function(genes_df) {
  h2_gene <- genes_df$gene_h2
  MAX_EXPONENT <- 10 
  h2 <- sum(h2_gene)
  offset <- max(h2_gene)
  integrand <- function(x) {
    exponent <- pmin(1/offset-1/x, MAX_EXPONENT)
    ifelse(x == 0, 0, x * exp(exponent) / h2)
  }
  inside_term <- integrand(h2_gene)
  poly <- h2 * (1/offset - log(sum(inside_term)))
}

polygenicity_functions <- list(
  softmax_polygenicity = softmax_polygenicity,
  effective_polygenicity = effective_polygenicity,
  entropy_polygenicity = entropy_polygenicity,
  effective_effect_var_polygenicity = effective_effect_var_polygenicity
)

calculate_true_polygenicity <- function(genes_df) {
  results <- list()
  for (metric_name in names(polygenicity_functions)) {
    metric_calculator_func <- polygenicity_functions[[metric_name]]
    metric_output <- metric_calculator_func(genes_df)
    results[[metric_name]] <- metric_output
  }
  return(results)
}

#' Meta-analyze true and estimated polygenicity results
#'
#' Calculates summary statistics on a log10 scale for true and estimated polygenicity.
#' For true results, it calculates the mean of log10-transformed values.
#' For estimated results, it calculates the mean of log10-transformed values and
#' the root-mean-square (RMS) of log10-transformed standard errors.
#' The formula for SE transformation is: se(log10(X)) = se(X) / (X * ln(10)).
#'
#' @param true_results_df Data frame of true polygenicity results, with 'abbreviation'.
#' @param estimated_results_df Data frame of estimated polygenicity results, with 'abbreviation',
#'        metric columns, and corresponding metric_se columns.
#' @return A data frame with meta-analyzed results, joined by 'abbreviation'.
meta_analyze_polygenicity <- function(true_results_df, estimated_results_df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for meta_analyze_polygenicity.")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("tibble package is required for meta_analyze_polygenicity.")
  }

  # Summarize true results: mean of log10 for all numeric columns
  true_summary <- true_results_df %>%
    dplyr::group_by(abbreviation) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = where(is.numeric),
        .fns = ~mean(log10(ifelse(.x > 0, .x, NA_real_)), na.rm = TRUE),
        .names = "{.col}_mean_log10"
      ),
      .groups = "drop_last"
    ) %>%
    dplyr::rename_with(~ paste0(., ".true"), -abbreviation)

  # Summarize estimated results
  estimated_summary <- estimated_results_df %>%
    dplyr::group_by(abbreviation) %>%
    dplyr::group_modify(~ {
      .x_data <- .x
      .x_numeric_cols <- .x_data %>% dplyr::select(where(is.numeric))
      
      current_value_cols <- names(.x_numeric_cols)[!grepl("_se$", names(.x_numeric_cols))]
      
      summary_row_list <- list()
      for (col_name in current_value_cols) {
        vals <- .x_numeric_cols[[col_name]]
        valid_val_indices <- vals > 0 & !is.na(vals)
        
        if(any(valid_val_indices)){
          summary_row_list[[paste0(col_name, "_mean_log10")]] <- mean(log10(vals[valid_val_indices]), na.rm = TRUE)
        } else {
          summary_row_list[[paste0(col_name, "_mean_log10")]] <- NA_real_
        }
        
        se_col_name <- paste0(col_name, "_se")
        if (se_col_name %in% names(.x_numeric_cols)) {
          ses <- .x_numeric_cols[[se_col_name]]
          # Indices where both value and SE are valid for transformation
          valid_transform_indices <- vals > 0 & !is.na(vals) & !is.na(ses)
          
          if(any(valid_transform_indices)){
            # Use only valid pairs for transformation
            vals_for_transform <- vals[valid_transform_indices]
            ses_for_transform <- ses[valid_transform_indices]
            
            log10_transformed_se_values <- ses_for_transform / (vals_for_transform * log(10)) # log(10) is ln(10) in R
            summary_row_list[[paste0(col_name, "_rms_log10_se")]] <- sqrt(mean(log10_transformed_se_values^2, na.rm = TRUE))
          } else {
            summary_row_list[[paste0(col_name, "_rms_log10_se")]] <- NA_real_
          }
        }
      }
      tibble::as_tibble(summary_row_list)
    }, .keep = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::rename_with(~ paste0(., ".estimated"), -abbreviation)

  # Merge the two summaries
  merged_results <- dplyr::full_join(true_summary, estimated_summary, by = "abbreviation")

  # Ensure original order of abbreviations from the 'true' set is preserved
  if (nrow(true_results_df) > 0 && "abbreviation" %in% names(true_results_df)) {
    ordered_abbreviations <- unique(as.character(true_results_df$abbreviation))
    merged_results$abbreviation <- as.character(merged_results$abbreviation)
    
    merged_results <- merged_results[match(ordered_abbreviations, merged_results$abbreviation), ]
    # Remove any all-NA rows that might have been introduced by match if an abbreviation was missing
    merged_results <- merged_results[rowSums(is.na(merged_results)) < (ncol(merged_results) -1), ] # -1 for abbreviation col
  }
  
  return(merged_results)
}
