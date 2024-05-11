#Functions for input/output and data preprocessing

process_data_trio <- function(input_data,
                              features) {

  if (!is.null(input_data$case_rate)) {
    input_data$expected_count = 2 * input_data$N * input_data$case_rate
  }

  if (any(is.na(input_data))) {
    stop("NAs present in input data, please check")
  } else if ( !is.null(features) & any(is.na(features))) {
    stop("NAs present in features, please check")
  }

  if (!is.null(features) & !(all(rownames(input_data) == rownames(features)))) {
    stop("features rownames do not match input data rownames, please check")
  }


return(input_data)
}
