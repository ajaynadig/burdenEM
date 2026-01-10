#' Calculate Delta Values
#'
#' Computes delta values for variants based on per-SD effect sizes and burden weights.
#'
#' @param beta_per_sd Numeric vector of per-standard-deviation effect sizes.
#' @param burden_weight Numeric vector of burden weights.
#' @return Numeric vector of delta values.
#' @export
calculate_delta_values <- function(beta_per_sd, burden_weight) {
  if (!is.numeric(beta_per_sd) || !is.numeric(burden_weight) || length(beta_per_sd) != length(burden_weight)) {
    stop("Inputs 'beta_per_sd' and 'burden_weight' must be numeric vectors of the same length.")
  }
  if (length(beta_per_sd) == 0) return(numeric(0))

  w_dot_w <- sum(burden_weight^2, na.rm = TRUE)
  if (w_dot_w == 0) return(rep(NA_real_, length(beta_per_sd)))

  w_dot_beta <- sum(burden_weight * beta_per_sd, na.rm = TRUE)
  gamma <- w_dot_beta / sqrt(w_dot_w)

  numerator <- beta_per_sd - gamma * burden_weight / sqrt(w_dot_w)
  denominator <- sqrt(1 - burden_weight^2 / w_dot_w)

  delta <- ifelse(denominator == 0, NA_real_, numerator / denominator)
  return(delta)
}

# --- Prepare Data for BurdenEM --- 

#' Prepare Gene-Level Data for BurdenEM Model Input
#'
#' Selects and renames columns from the gene-level data frame to match the 
#' expected input format for the BurdenEM model functions.
#'
#' @param gene_level_df A data frame containing gene-level aggregated data. 
#'                      Must include 'gamma_per_sd' and 'gene_intercept' columns.
#' @return A data frame with columns 'effect_estimate' (from 'gamma_per_sd') and 
#'         'effect_se' (from 'gene_intercept'), along with all other columns from the input.
#' @importFrom dplyr %>% mutate select any_of
#' @export
specify_gene_effect_sizes <- function(gene_level_df, per_allele_effects=FALSE, correct_for_ld=FALSE) {
    required_cols <- c("gamma_per_sd", "gene_intercept", "burden_score")
    if (!all(required_cols %in% names(gene_level_df))) {
        stop("Input dataframe must contain columns: ", paste(required_cols, collapse = ", "))
    }

    if (correct_for_ld) {
        gene_level_df <- gene_level_df %>% 
              dplyr::mutate(
                  ld_correction_factor = ifelse(is.na(burden_score_ld) | burden_score_ld == 0, 1, burden_score_ld / burden_score_no_ld),
                  gamma_per_sd = gamma_per_sd / sqrt(ld_correction_factor), 
                  burden_score = burden_score * ld_correction_factor
              ) 
    }
    message("Mean LD correction factor: ", mean(gene_level_df$ld_correction_factor, na.rm = TRUE))


    # Create the effect estimate and standard error columns
    if (per_allele_effects) {
          gene_level_df <- gene_level_df %>% 
              dplyr::mutate(
                  effect_estimate = gamma_per_sd / sqrt(burden_score), 
                  effect_se = sqrt(gene_intercept / burden_score)
              ) 
    }
    else{
          gene_level_df <- gene_level_df %>% 
              dplyr::mutate(
                  effect_estimate = gamma_per_sd, 
                  effect_se = sqrt(gene_intercept)
              )
        }

    return(gene_level_df)
}

calculate_variant_intercept <- function(variant_data, frequency_bin_edges) {
  # --- Input Validation ---
  # Basic checks
  if (!is.data.frame(variant_data)) stop("'variant_data' must be a data frame.")
  if (!is.numeric(frequency_bin_edges) || length(frequency_bin_edges) < 2) stop("'frequency_bin_edges' must be a numeric vector with at least two elements.")
  if (any(diff(frequency_bin_edges) <= 0)) stop("'frequency_bin_edges' must be strictly increasing.")

  # Check required columns
  required_cols <- c("beta_per_sd", "AF", "functional_category") # Using beta_per_sd now
  missing_cols <- setdiff(required_cols, names(variant_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in 'variant_data':", paste(missing_cols, collapse = ", ")))
  }

  # --- Data Preparation ---
  # Assign frequency bins
  # Include lowest value, make right side open except for the last bin
  freq_labels <- paste0("[", head(frequency_bin_edges, -1), ",", tail(frequency_bin_edges, -1), ")")
  # Ensure the last label is closed: [bin_n-1, bin_n]
  freq_labels[length(freq_labels)] <- gsub("\\)$", "\\]", freq_labels[length(freq_labels)])

  variant_data <- variant_data %>%
    dplyr::mutate(
      frequency_bin = cut(AF, breaks = frequency_bin_edges, labels = freq_labels, right = FALSE, include.lowest = TRUE)
    )

  # --- Calculate Intercepts per Bin ---
  # Group and calculate weighted averages
  intercept_calc <- variant_data %>%
    dplyr::group_by(gene, functional_category, frequency_bin) %>%
    dplyr::summarise(
      gamma_per_allele = sum(beta_per_sd * sqrt(variant_variance)) / sum(variant_variance),
      beta_centered = beta_per_sd - gamma_per_allele * sqrt(variant_variance),
      n_variants = dplyr::n(),
      variance = mean(beta_centered^2 / (1 - variant_variance / sum(variant_variance))),
      variance = ifelse(n_variants > 1, variance, 0), # degrees of freedom are n_variants-1
      .groups = 'drop'
    )
  
  # Aggregate genes
  intercept_summary <- intercept_calc %>%
    dplyr::group_by(functional_category, frequency_bin) %>%
    dplyr::summarise(
      intercept = sum(variance * (n_variants - 1)) / sum(n_variants - 1),
      .groups = 'drop'
    )%>%
    dplyr::select(functional_category, frequency_bin, intercept) %>%
    dplyr::arrange(functional_category, frequency_bin)

  # Check for any NA intercepts which might indicate issues
  if(any(is.na(intercept_summary$intercept))) {
      warning("Some intercepts resulted in NA. This might be due to bins with zero weight sum or calculation issues.")
  }

  return(intercept_summary)
}