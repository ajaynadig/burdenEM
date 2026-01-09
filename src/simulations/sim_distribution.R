# sim_distribution.R: Functions for calculating and meta-analyzing true gene distribution metrics.

#' Calculate true distribution of genes contributing to heritability
#'
#' Determines the number of genes required to explain specified percentages of total true heritability.
#'
#' @param genes_df A data frame containing gene-level data, including a `gene_h2` column 
#'   representing the true heritability contributed by each gene.
#' @return A data frame with columns `variance` (the target heritability percentage) and 
#'   `needed_genes` (the number of top-contributing genes required to reach that H2 percentage).
calculate_true_distribution <- function(genes_df, gene_proportions) {
  if (!"gene_h2" %in% names(genes_df)) {
    stop("genes_df must contain a 'gene_h2' column.")
  }

  if (nrow(genes_df) == 0) {
    stop("genes_df must not be empty.")
  }

  total_num_genes <- nrow(genes_df)
  gene_counts <- floor(total_num_genes * gene_proportions)

  # Sort genes by their contribution to h2 in descending order
  sorted_gene_h2 <- sort(genes_df$gene_h2, decreasing = TRUE)
  cum_h2 <- c(0, cumsum(sorted_gene_h2))
  cum_h2 <- cum_h2 / cum_h2[length(cum_h2)]
  results_df <- data.frame(
    needed_genes = gene_counts,
    variance = cum_h2[1 + gene_counts],
    pergene_variance = cum_h2[1 + pmax(gene_counts, 1)] - cum_h2[pmax(gene_counts, 1)]
  )

  return(results_df)
}



#' Meta-analyze true and estimated distribution results
#'
#' This function takes two data frames, one with true distribution values and one with
#' estimated values, summarizes them across replicates, and merges them.
#'
#' @param true_results_df A data frame of true distribution results. Must contain
#'   'abbreviation', 'variance', and 'needed_genes' columns.
#' @param estimated_results_df A data frame of estimated distribution results. Must
#'   contain 'abbreviation', 'variance', and 'needed_genes' columns.
#' @return A data frame with meta-analyzed results, containing mean and sd for
#'   'needed_genes' from both true and estimated inputs, merged by 'abbreviation' and 'variance'.
meta_analyze_distribution <- function(true_results_df, estimated_results_df) {
  
  # Corrected summarize_results using a clearer rename approach
  summarize_results <- function(df, suffix) {
    if (nrow(df) == 0) {
      stop("Data frame is empty, cannot summarize results.")
    }
    summary_df <- df %>% 
      dplyr::group_by(abbreviation, needed_genes, .drop = FALSE) %>% 
      dplyr::summarise(
        !!paste0("variance_mean", suffix) := mean(variance, na.rm = TRUE),
        !!paste0("variance_sd", suffix) := sd(variance, na.rm = TRUE),
        .groups = "drop"
      )
    
    return(summary_df)
  }

  summarize_se <- function(df, suffix) {
    if (nrow(df) == 0) {
      stop("Data frame is empty, cannot summarize results.")
    }
    df %>% 
      dplyr::group_by(abbreviation, needed_genes, .drop = FALSE) %>% 
      dplyr::summarise(
        !!paste0("se_needed_genes_mean", suffix) := mean(needed_genes_se, na.rm = TRUE),
        !!paste0("se_variance_mean", suffix) := mean(variance_se, na.rm = TRUE),
        .groups = "drop"
      )
  }

  true_summary <- summarize_results(true_results_df, ".true")
  estimated_summary <- summarize_results(estimated_results_df, ".estimated")
  # se_summary <- summarize_se(estimated_results_df, ".estimated")
  
  merged_results <- dplyr::full_join(true_summary, estimated_summary, by = c("abbreviation", "needed_genes"))
  # merged_results <- dplyr::full_join(merged_results, se_summary, by = c("abbreviation", "variance"))
  
  # if (nrow(true_results_df) > 0 && "abbreviation" %in% names(true_results_df) && "variance" %in% names(true_results_df)){
  #   ordered_keys_df <- unique(true_results_df[, c("abbreviation", "variance")])
  #   if (nrow(ordered_keys_df) > 0 && all(c("abbreviation", "variance") %in% names(merged_results))) {
  #       merged_results <- dplyr::left_join(ordered_keys_df, merged_results, by = c("abbreviation", "variance"))
  #   }
  # }
  
  return(merged_results)
}
