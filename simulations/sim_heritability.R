


total_h2 <- function(genes_df) {
  h2_gene <- genes_df$gene_h2
  h2 <- sum(h2_gene)
}

positive_h2 <- function(genes_df) {
  h2_gene <- genes_df$gene_h2[genes_df$effectSize > 0]
  h2 <- sum(h2_gene)
}

negative_h2 <- function(genes_df) {
  h2_gene <- genes_df$gene_h2[genes_df$effectSize < 0]
  h2 <- sum(h2_gene)
}

heritability_functions <- list(
  total_h2 = total_h2,
  positive_h2 = positive_h2,
  negative_h2 = negative_h2
)

calculate_true_heritability <- function(genes_df) {
  results <- list()
  for (metric_name in names(heritability_functions)) {
    metric_calculator_func <- heritability_functions[[metric_name]]
    metric_output <- metric_calculator_func(genes_df)
    results[[metric_name]] <- metric_output
  }
  return(results)
}

#' Meta-analyze true and estimated heritability results
#'
#' This function takes two data frames, one with true heritability values and one with
#' estimated values, summarizes them across replicates, and merges them.
#'
#' @param true_results_df A data frame of true heritability results. Must contain
#'   an 'abbreviation' column and numeric metric columns.
#' @param estimated_results_df A data frame of estimated heritability results. Must
#'   contain an 'abbreviation' column and numeric metric columns.
#' @return A data frame with meta-analyzed results, containing mean and sd for
#'   each metric from both true and estimated inputs, merged by abbreviation.
meta_analyze_heritability <- function(true_results_df, estimated_results_df) {
  
  # Add replicate IDs to join on. Assumes rows are in the same replicate order for each abbreviation.
  true_df_with_id <- true_results_df %>%
    dplyr::group_by(abbreviation) %>%
    dplyr::mutate(replicate_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  estimated_df_with_id <- estimated_results_df %>%
    dplyr::group_by(abbreviation) %>%
    dplyr::mutate(replicate_id = dplyr::row_number()) %>%
    dplyr::ungroup()

  # To avoid name clashes, suffix the columns of the true dataframe
  true_metric_cols <- names(true_results_df)[sapply(true_results_df, is.numeric)]
  true_df_renamed <- true_df_with_id %>%
    dplyr::rename_with(~ paste0(., ".true"), .cols = all_of(true_metric_cols))
  
  # Now join with the estimated dataframe
  merged_replicates_df <- dplyr::full_join(
    estimated_df_with_id, 
    true_df_renamed, 
    by = c("abbreviation", "replicate_id")
  )
  
  # Calculate residuals for each metric
  estimated_metric_names <- estimated_results_df %>% 
    dplyr::select(where(is.numeric) & !ends_with(c("_se", ".se"))) %>% 
    names()
  
  for (metric in estimated_metric_names) {
    true_col_name <- paste0(metric, ".true")
    residual_col_name <- paste0(metric, "_residual")
    if (true_col_name %in% names(merged_replicates_df)) {
      merged_replicates_df <- merged_replicates_df %>%
        dplyr::mutate(!!residual_col_name := .data[[metric]] - .data[[true_col_name]])
    }
  }
  
  # Find the columns to summarize
  estimated_se_names <- estimated_results_df %>% dplyr::select(ends_with(c("_se", ".se"))) %>% names()
  residual_names <- names(merged_replicates_df)[endsWith(names(merged_replicates_df), "_residual")]

  summary_df <- merged_replicates_df %>%
    dplyr::group_by(abbreviation) %>%
    dplyr::summarise(
      # Mean and SD of estimated metrics
      dplyr::across(
        .cols = all_of(estimated_metric_names),
        .fns = list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}.estimated"
      ),
      # Mean and SD of true metrics
      dplyr::across(
        .cols = all_of(paste0(true_metric_cols, ".true")),
        .fns = list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
        .names = "{sub('\\\\.true$', '', .col)}_{.fn}.true"
      ),
      # Mean and RMS of estimated SEs
      dplyr::across(
        .cols = all_of(estimated_se_names),
        .fns = list(mean = ~mean(.x, na.rm = TRUE), rms = ~sqrt(mean(.x^2, na.rm = TRUE))),
        .names = "{.col}_{.fn}.estimated"
      ),
      # SD of residuals
      dplyr::across(
        .cols = all_of(residual_names),
        .fns = ~sd(.x, na.rm = TRUE),
        .names = "{.col}_sd"
      ),
      .groups = "drop"
    )
  
  # Ensure original order of abbreviations from the 'true' set is preserved
  if (nrow(true_results_df) > 0 && "abbreviation" %in% names(true_results_df)) {
    ordered_abbreviations <- unique(as.character(true_results_df$abbreviation))
    summary_df$abbreviation <- as.character(summary_df$abbreviation)
    
    summary_df <- summary_df[match(ordered_abbreviations, summary_df$abbreviation), ]
    # Remove any all-NA rows that might have been introduced by match if an abbreviation was missing
    summary_df <- summary_df[rowSums(is.na(summary_df)) < (ncol(summary_df) -1), ] # -1 for abbreviation col
  }
  
  return(summary_df)
}
