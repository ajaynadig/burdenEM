#' Create One-Hot Encoded Quantile Bins for a Gene Feature (Handling NAs)
#'
#' Calculates quantile bins for a specified column in a gene dataframe and
#' returns a dataframe with one-hot encoded columns representing bin membership.
#' NA values are assigned to a separate "NA" bin (6th bin) to keep the feature
#' matrix binary/one-hot, which is required for EM numerical stability.
#'
#' @param gene_df A dataframe containing gene-level data. Must contain the column
#'   specified by `column_name`.
#' @param column_name A character string specifying the name of the column in
#'   `gene_df` to be binned. Default is "lof_oe".
#' @param num_bins An integer specifying the number of quantile bins to create.
#'   Default is 5.
#' @param verbose Logical, whether to print messages about NA handling. Default TRUE.
#'
#' @return The input dataframe `gene_df` with an added list-column named `features`.
#'   Each element of the `features` column is a 1-row matrix containing the named
#'   one-hot encoded bin membership.
#' @export
#' @importFrom stats quantile model.matrix as.formula
#' @importFrom dplyr %>% mutate select any_of filter pull arrange rename_with tibble as_tibble starts_with case_when row_number all_of
#' @importFrom rlang := !! sym
#' @importFrom stringr str_replace
get_bins <- function(gene_df, column_name = "lof_oe", num_bins = 5, verbose = TRUE) {
  stopifnot(
    is.data.frame(gene_df),
    is.character(column_name), length(column_name) == 1,
    column_name %in% names(gene_df),
    is.numeric(num_bins), num_bins > 0, floor(num_bins) == num_bins,
    is.numeric(gene_df[[column_name]])
  )

  col_values <- gene_df[[column_name]]
  na_mask <- is.na(col_values)
  if (verbose && any(na_mask)) {
    message("Found ", sum(na_mask), " NA values in '", column_name,
            "'. Assigning to separate NA bin.")
  }

  # Calculate quantiles and cut into bins
  quantile_breaks <- unique(stats::quantile(col_values, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE))
  num_actual_bins <- length(quantile_breaks) - 1
  bin_labels <- paste0(column_name, "_bin_", seq_len(num_actual_bins))
  na_bin_label <- paste0(column_name, "_bin_NA")

  if (verbose) {
    message("Quantile breaks for '", column_name, "': ", paste(quantile_breaks, collapse = ", "))
  }
  bins <- cut(col_values, breaks = quantile_breaks,
              labels = bin_labels, include.lowest = TRUE)

  # One-hot encode non-NA bin values
  non_na_df <- data.frame(bin = bins[!na_mask])
  one_hot <- model.matrix(~ bin - 1, data = non_na_df)
  colnames(one_hot) <- gsub("^bin", "", colnames(one_hot))

  # All labels including NA bin
  all_labels <- c(bin_labels, na_bin_label)
  num_total_bins <- length(all_labels)

  # Construct features column: one-hot with NA bin as extra column
  features <- vector("list", length(bins))
  non_na_counter <- 0
  for (i in seq_along(bins)) {
    vec <- rep(0, num_total_bins)
    names(vec) <- all_labels
    if (is.na(bins[i])) {
      vec[na_bin_label] <- 1
    } else {
      non_na_counter <- non_na_counter + 1
      matched_bin <- as.character(bins[i])
      vec[matched_bin] <- 1
    }
    features[[i]] <- matrix(vec, nrow = 1, dimnames = list(NULL, all_labels))
  }

  if (verbose) {
    message(sprintf("Feature bins: %d quantile bins + 1 NA bin (%d genes in NA bin)",
                    num_actual_bins, sum(na_mask)))
  }

  gene_df$features <- features
  return(gene_df)
}
