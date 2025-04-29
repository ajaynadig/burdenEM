#' Aggregate Variant Data to Gene Level with Integrated Processing
#'
#' Processes variant-level summary statistics, calculates necessary metrics
#' ('beta_per_sd', 'gamma_per_sd', 'burden_score', 'gene_intercept'),
#' calculates and joins intercepts, and filters the final gene-level data.
#'
#' @param variant_data Dataframe with required columns: 'gene', 'AF' (Allele Frequency), 'beta', and 'variant_variance'.
#' @param frequency_bin_edges Numeric vector defining the edges for allele frequency bins.
#' @return A dataframe aggregated per gene, with columns: 'gene', 'gamma_per_sd', 'gene_intercept',
#'         'burden_score', 'n_variants'. Assumes input is for a single functional category.
#'         Returns an empty dataframe if processing fails.
#' @importFrom dplyr %>% group_by summarize mutate filter n left_join select rename all_of arrange distinct case_when if_else
#' @importFrom tidyr separate
#' @importFrom stats setNames complete.cases sd
process_variant_to_gene <- function(variant_data, frequency_bin_edges) {

  # --- Input Validation ---
  required_cols <- c("gene", "AF", "beta", "variant_variance") # functional_category implicitly handled
  if (!is.data.frame(variant_data) || !all(required_cols %in% names(variant_data))) {
    stop(paste("'variant_data' must be a dataframe containing columns:", paste(required_cols, collapse=", ")))
  }
  if (!is.numeric(frequency_bin_edges) || length(frequency_bin_edges) < 2) {
    stop("'frequency_bin_edges' must be a numeric vector with at least two values.")
  }

  # --- 1. Calculate beta_per_sd ---
  message("Calculating 'beta_per_sd'...")
  variant_data <- variant_data %>%
      mutate(beta_per_sd = beta * sqrt(variant_variance))

  # --- 2. Calculate Variant Intercepts ---
  intercept_summary <- calculate_variant_intercept(variant_data, frequency_bin_edges)

  # --- Add Frequency Bins to variant_data (needed for join) ---
  # Include lowest value, make right side open except for the last bin
  freq_labels <- paste0("[", head(frequency_bin_edges, -1), ",", tail(frequency_bin_edges, -1), ")")
  # Ensure the last label is closed: [bin_n-1, bin_n]
  freq_labels[length(freq_labels)] <- gsub("\\)$", "\\]", freq_labels[length(freq_labels)])

  variant_data <- variant_data %>%
    dplyr::mutate(
      frequency_bin = cut(AF, breaks = frequency_bin_edges, labels = freq_labels, right = FALSE, include.lowest = TRUE)
    ) %>%
    dplyr::filter(!is.na(frequency_bin)) # Remove variants that fall outside defined bins

  # Join intercepts back to variant data
  variant_data <- variant_data %>%
      dplyr::left_join(intercept_summary %>% select(functional_category, frequency_bin, intercept),
                       by = c("functional_category", "frequency_bin"))

  # --- 3. Group by Gene and Summarize ---
  message("Aggregating results per gene...")
  gene_level_summary <- variant_data %>%
    dplyr::group_by(gene, functional_category) %>%
    dplyr::summarize(
      gamma_per_sd = sum(beta_per_sd * sqrt(variant_variance)),
      burden_score = sum(variant_variance),
      gene_intercept = sum(intercept * variant_variance),
      n_variants = n(),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      gamma_per_sd = if_else(burden_score > 0, gamma_per_sd / sqrt(burden_score), 0),
      gene_intercept = if_else(burden_score > 0, gene_intercept / burden_score, 0)
    )

  # Select final columns and ungroup
  gene_level_summary <- gene_level_summary %>%
    dplyr::select(
      gene,
      functional_category,
      gamma_per_sd,
      gene_intercept,
      burden_score,
      n_variants
    ) %>%
    dplyr::ungroup()


  message("Finished processing variant data to gene level.")
  return(gene_level_summary)
}
