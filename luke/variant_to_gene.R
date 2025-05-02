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
    dplyr::filter(!is.na(frequency_bin)
    ) %>%
      dplyr::left_join(intercept_summary %>% select(functional_category, frequency_bin, intercept), 
                       by = c("functional_category", "frequency_bin"))

  # --- 3. Group by Gene and Summarize --- 
  message("Aggregating results per gene...")
  gene_level_summary <- variant_data %>%
    group_by(gene, functional_category) %>% 
    summarize(
      gamma_per_sd = sum(beta_per_sd * sqrt(variant_variance)),
      burden_score = sum(variant_variance), 
      gene_intercept = sum(intercept * variant_variance),
      n_variants = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      gamma_per_sd = if_else(burden_score > 0, gamma_per_sd / sqrt(burden_score), 0),
      gene_intercept = if_else(burden_score > 0, gene_intercept / burden_score, 0)
    )

  # Select final columns and ungroup
  gene_level_summary <- gene_level_summary %>% 
    select(
      gene,
      functional_category,
      gamma_per_sd,
      gene_intercept,
      burden_score,
      n_variants
    ) %>% 
    ungroup()


  message("Finished processing variant data to gene level.")
  return(gene_level_summary)
}

#' Process variant-level data (including pre-calculated likelihoods) to gene-level aggregates.
#'
#' This function aggregates variant-level statistics and pre-calculated log-likelihoods
#' to the gene x functional category level. It assumes the input data frame
#' contains columns for 'gene', 'functional_category', 'variant_variance',
#' and additional columns representing the log-likelihood for each grid point.
#'
#' @param variant_data_with_ll A data frame containing variant information and
#'   log-likelihoods. Must include columns 'gene', 'functional_category',
#'   'variant_variance', and columns for log-likelihoods (e.g., "LL_1", "LL_2", ...).
#'
#' @return A list containing two elements:
#'   \describe{
#'     \item{gene_level_stats}{A data frame summarized at the gene x category
#'       level with columns: 'gene', 'functional_category', 'burden_score'
#'       (sum of variant variances), 'n_variants'.}
#'     \item{gene_level_likelihoods}{A matrix containing the summed log-likelihoods
#'       for each gene x category combination. Rows are named 'gene:functional_category',
#'       columns correspond to the input likelihood columns.}
#'   }
#' @importFrom dplyr group_by summarize n ungroup select any_of across all_of ends_with contains relocate
#' @importFrom tidyr unite pivot_longer pivot_wider
#' @export
process_variant_to_gene_binary <- function(variant_data_with_ll) {

    # --- Input Checks ---
    # Identify expected non-likelihood columns needed for stats summary
    stat_cols <- c("gene", "functional_category", "variant_variance") # Add others if needed

    # Identify likelihood columns specifically by prefix
    likelihood_col_names <- names(variant_data_with_ll)[startsWith(names(variant_data_with_ll), "LL_")]

    # Check that required stat columns are present
    missing_cols <- setdiff(stat_cols, names(variant_data_with_ll))
    if (length(missing_cols) > 0) {
        stop(paste("Missing required columns in variant_data_with_ll:", paste(missing_cols, collapse=", ")))
    }
    if (length(likelihood_col_names) == 0) {
        stop("No columns identified as likelihood columns in variant_data_with_ll.")
    }
    # Check if likelihood columns are numeric
     likelihood_cols_are_numeric <- all(sapply(variant_data_with_ll[, likelihood_col_names, drop = FALSE], is.numeric))
     if(!likelihood_cols_are_numeric){
        stop("Identified likelihood columns must be numeric.")
     }


    # --- Aggregate to Gene x Category Level ---
    gene_summary <- variant_data_with_ll %>%
        dplyr::group_by(gene, functional_category) %>%
        dplyr::summarize(
            # Sum variant variances to get a simple burden score for stats output
            burden_score = sum(variant_variance, na.rm = TRUE),
            # Count number of variants per gene/category
            n_variants = dplyr::n(),
            # Sum the log-likelihoods across all variants in the group for each grid point (column)
            dplyr::across(dplyr::all_of(likelihood_col_names), ~ sum(.x, na.rm = TRUE)),
            .groups = 'drop'
        )

    # --- Separate Stats and Likelihoods ---
    # Gene-level statistics
    gene_level_stats <- gene_summary %>%
        dplyr::select(gene, functional_category, burden_score, n_variants) # Add more stats cols if summarized above

    # Gene-level likelihoods matrix
    # Create unique row names (gene:category)
    gene_summary <- gene_summary %>%
        tidyr::unite("gene_category", gene, functional_category, sep = ":", remove = FALSE)

    # Extract likelihood columns into a matrix
    # Ensure columns are ordered correctly if they came in like LL_1, LL_10, LL_2
    # If names are like "LL_neg2.0", "LL_neg1.8", etc. sorting works. If "LL_1", "LL_2", "LL_10", need numeric sort.
    # Assuming column names allow correct alphabetic sorting or are already ordered.
    gene_level_likelihoods <- gene_summary %>%
        dplyr::select(gene_category, dplyr::all_of(likelihood_col_names)) %>%
        tibble::column_to_rownames("gene_category") %>%
        as.matrix()

    # --- Return Results ---
    return(list(
        gene_level_stats = gene_level_stats,
        gene_level_likelihoods = gene_level_likelihoods
    ))
}
