# sim_distribution.R: Functions for calculating and meta-analyzing true gene distribution metrics.

#' Calculate true distribution of genes contributing to heritability
#'
#' Determines the number of genes required to explain specified percentages of total true heritability.
#'
#' @param genes_df A data frame containing gene-level data, including a `gene_h2` column 
#'   representing the true heritability contributed by each gene.
#' @return A data frame with columns `H2_Percent` (the target heritability percentage) and 
#'   `Num_Genes` (the number of top-contributing genes required to reach that H2 percentage).
calculate_true_distribution <- function(genes_df) {
  if (!"gene_h2" %in% names(genes_df)) {
    stop("genes_df must contain a 'gene_h2' column.")
  }

  # Filter for genes with positive heritability contribution
  genes_df_filtered <- genes_df[!is.na(genes_df$gene_h2) & genes_df$gene_h2 > 0, ]

  if (nrow(genes_df_filtered) == 0) {
    return(data.frame(H2_Percent = numeric(0), Num_Genes = numeric(0)))
  }

  total_true_h2 <- sum(genes_df_filtered$gene_h2, na.rm = TRUE)

  # Sort genes by their contribution to h2 in descending order
  sorted_genes_df <- genes_df_filtered[order(-genes_df_filtered$gene_h2), ]

  test_heritability_proportions <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999)
  
  results_list <- list()

  if (total_true_h2 == 0) {
      for (hp_g in test_heritability_proportions) {
          results_list[[length(results_list) + 1]] <- data.frame(
              H2_Percent = hp_g * 100,
              Num_Genes = 0 
          )
      }
      return(dplyr::bind_rows(results_list))
  }

  for (hp_g in test_heritability_proportions) {
    target_h2_value <- hp_g * total_true_h2
    cumulative_h2 <- 0
    genes_counted <- 0
    
    if (target_h2_value == 0) {
        results_list[[length(results_list) + 1]] <- data.frame(
            H2_Percent = hp_g * 100,
            Num_Genes = 0
        )
        next
    }

    for (i in 1:nrow(sorted_genes_df)) {
      cumulative_h2 <- cumulative_h2 + sorted_genes_df$gene_h2[i]
      genes_counted <- genes_counted + 1
      if (cumulative_h2 >= target_h2_value) {
        break
      }
    }
    
    results_list[[length(results_list) + 1]] <- data.frame(
      H2_Percent = hp_g * 100,
      Num_Genes = genes_counted
    )
  }
  
  return(dplyr::bind_rows(results_list))
}

#' Meta-analyze true and estimated distribution results
#'
#' This function takes two data frames, one with true distribution values and one with
#' estimated values, summarizes them across replicates, and merges them.
#'
#' @param true_results_df A data frame of true distribution results. Must contain
#'   'abbreviation', 'H2_Percent', and 'Num_Genes' columns.
#' @param estimated_results_df A data frame of estimated distribution results. Must
#'   contain 'abbreviation', 'H2_Percent', and 'Num_Genes' columns.
#' @return A data frame with meta-analyzed results, containing mean and sd for
#'   'Num_Genes' from both true and estimated inputs, merged by 'abbreviation' and 'H2_Percent'.
meta_analyze_distribution <- function(true_results_df, estimated_results_df) {
  
  summarize_results <- function(df, suffix) {
    if (nrow(df) == 0) {
        return(data.frame(abbreviation = character(0), H2_Percent = numeric(0), Num_Genes_mean = numeric(0), Num_Genes_sd = numeric(0)))
    }
    df %>% 
      dplyr::group_by(abbreviation, H2_Percent, .drop = FALSE) %>% 
      dplyr::summarise(
        Num_Genes_mean = mean(Num_Genes, na.rm = TRUE),
        Num_Genes_sd = sd(Num_Genes, na.rm = TRUE),
        .groups = "drop" # Ungroup completely after summarise
      ) %>% 
      dplyr::rename_with(~ paste0("Num_Genes", suffix, if(. == "Num_Genes_mean") "_mean" else if (.=="Num_Genes_sd") "_sd" else sub("Num_Genes_?", "", .)), starts_with("Num_Genes")) %>% # More robust renaming
      dplyr::rename_with(~ paste0(sub("Num_Genes", "", .x), suffix), .cols = starts_with("Num_Genes_")) # simpler rename: Num_Genes_mean -> Num_Genes_mean.true
  }

  # Corrected summarize_results using a clearer rename approach
  summarize_results_corrected <- function(df, suffix) {
    if (nrow(df) == 0) {
      return(data.frame(abbreviation = character(0), H2_Percent = numeric(0)))
    }
    summary_df <- df %>% 
      dplyr::group_by(abbreviation, H2_Percent, .drop = FALSE) %>% 
      dplyr::summarise(
        mean_val = mean(Num_Genes, na.rm = TRUE),
        sd_val = sd(Num_Genes, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Dynamically create new column names
    mean_col_name <- paste0("Num_Genes_mean", suffix)
    sd_col_name <- paste0("Num_Genes_sd", suffix)
    
    names(summary_df)[names(summary_df) == "mean_val"] <- mean_col_name
    names(summary_df)[names(summary_df) == "sd_val"] <- sd_col_name
    
    return(summary_df)
  }

  true_summary <- summarize_results_corrected(true_results_df, ".true")
  estimated_summary <- summarize_results_corrected(estimated_results_df, ".estimated")
  
  merged_results <- dplyr::full_join(true_summary, estimated_summary, by = c("abbreviation", "H2_Percent"))
  
  if (nrow(true_results_df) > 0 && "abbreviation" %in% names(true_results_df) && "H2_Percent" %in% names(true_results_df)){
    ordered_keys_df <- unique(true_results_df[, c("abbreviation", "H2_Percent")])
    if (nrow(ordered_keys_df) > 0 && all(c("abbreviation", "H2_Percent") %in% names(merged_results))) {
        merged_results <- dplyr::left_join(ordered_keys_df, merged_results, by = c("abbreviation", "H2_Percent"))
    }
  }
  
  return(merged_results)
}
