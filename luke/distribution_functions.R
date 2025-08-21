# luke/distribution_functions.R
source("R/model.R")

#' PMF of the effect size distribution over genes evaluated at grid points
#'
#' @param model A model object
#' @return A numeric vector representing the PMF.
gene_PMF_from_prior <- function(model) {
  features <- do.call(rbind, model$df$features)
  gene_component_contributions <- features %*% model$delta
  weights <- colMeans(gene_component_contributions)
  pmf_matrix <- weights %*% model$components
  return(as.numeric(pmf_matrix))
}

#' PMF of the effect size distribution over genes evaluated at grid points
#'
#' @param model A model object
#' @return A numeric vector representing the PMF.
gene_PMF_from_posterior <- function(model) {
  pmf_matrix <- grid_posterior(model)
  return(colMeans(pmf_matrix))
}

gene_PMF <- gene_PMF_from_posterior


#' PMF of the effect size distribution over components of variance evaluated at grid points
#' 
#' @param model A model object
#' @return A numeric vector representing the PMF.
variance_PMF <- function(model) {
  gene_pmf_values <- gene_PMF(model)
  var_contributions <- (model$grid^2) * gene_pmf_values # TODO: generalize to models with other h2 functions
  total_variance <- sum(var_contributions)
  
  if (abs(total_variance) < 1e-12) { # Check if total_second_moment is effectively zero
    warning("Total second moment (sum of grid_i^2 * PMF_i) is close to zero. Returning a PMF of zeros.")
    return(rep(0, length(model$grid)))
  }
  
  return(var_contributions / total_variance)
}

#' CDF of the gene effect size distribution evaluated at grid points
#'
#' @param model A model object.
#' @return A numeric vector representing the CDF.
#' @export
gene_CDF <- function(model) {
  pmf <- gene_PMF(model)
  if (is.null(pmf)) {
    return(NULL)
  }
  cdf <- cumsum(pmf)
  # Ensure CDF ends at 1 (it should due to PMF summing to 1, but good practice for floating points)
  # cdf[length(cdf)] <- 1.0 
  return(cdf)
}

#' CDF of the variance distribution evaluated at grid points
#'
#'  @param model A model object.
#' @return A numeric vector representing the CDF.
#' @export
variance_CDF <- function(model) {
  pmf <- variance_PMF(model)
  if (is.null(pmf)) {
    return(NULL)
  }
  cdf <- cumsum(pmf)
  return(cdf)
}

#' Creates a piecewise linear function interpolating between points on a grid.
#' 
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @return A piecewise linear function mapping f(x[i]) -> y[i].
#' @export
as_piecewise_linear_smooth <- function(x, y) {
  n <- length(y)
  y_smoothed <- y
  for (i in 1:5) {
    y_smoothed[2:(n-1)] <- 0.5 * (y_smoothed[1:(n-2)] + y_smoothed[3:n])
  }
  return(as_piecewise_linear(x, y_smoothed))
  # smoothed <- smooth.spline(x, y)
  # as_piecewise_linear(smoothed$x, smoothed$y)
}

#' Creates a piecewise linear function interpolating between points on a grid.
#' 
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @return A piecewise linear function mapping f(x[i]) -> y[i].
#' @export
as_piecewise_linear <- function(x, y) {
  approxfun(x, y, method = "linear", rule = 1, ties = "ordered", f = 1)
}

#' Returns the function g(t) = -f(t)
#' 
#' @param f A function.
#' @return A function g(t) = -f(t).
#' @export
get_negated <- function(f) {
  function(t) -f(t)
}

# Aggregates weights for consecutive duplicate values. Assumes `values` are already sorted.
aggregate_weights_by_value <- function(values, weights) {
  r <- rle(values)
  unique_values <- r$values
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  summed_weights <- mapply(function(a, b) sum(weights[a:b]), starts, ends)
  list(values = unique_values, weights = summed_weights)
}

#' Returns the function g(t) = f(-t)
#' 
#' @param f A function.
#' @return A function g(t) = f(-t).
#' @export
get_reversed <- function(f) {
  function(t) f(-t)
}

#' The CDF of a transformed random variable, approximated as a piecewise linear function.
#' Also can be used to get a derivate of the CDF, passing the derivative of the PMF.
#' 
#' @param grid_values A numeric vector of x-coordinates.
#' @param grid_pmf The PMF evaluated at the grid values.
#' @param f A transformation to be applied to the grid values.
#' @param right_tail A logical value indicating whether to return the right tailed CDF.
#' @return The CDF of the transformed random variable.
#' @export
get_piecewise_linear_cdf <- function(grid_values, grid_pmf, f=function(t) t, right_tail=FALSE) {
  fn <- if (right_tail) get_negated(f) else f
  f_values <- fn(grid_values)
  argsort <- order(f_values)
  f_values <- f_values[argsort]
  grid_pmf <- grid_pmf[argsort]
  agg <- aggregate_weights_by_value(f_values, grid_pmf)
  result_fn <- as_piecewise_linear(agg$values, cumsum(agg$weights))
  result <- if (right_tail) get_reversed(result_fn) else result_fn
  return(result)
}

#' The quantile function of a transformed random variable, approximated as a piecewise linear function
#' 
#' @param grid_values A numeric vector of x-coordinates.
#' @param grid_pmf The PMF evaluated at the grid values.
#' @param f A transformation to be applied to the grid values.
#' @param right_tail A logical value indicating whether to return the right tailed QF.
#' @return The quantile function of the transformed random variable.
#' @export
get_piecewise_linear_qf <- function(grid_values, grid_pmf, f=function(t) t, right_tail=FALSE) {
  fn <- if (right_tail) get_negated(f) else f
  f_values <- fn(grid_values)
  argsort <- order(f_values)
  f_values <- f_values[argsort]
  grid_pmf <- grid_pmf[argsort]
  agg <- aggregate_weights_by_value(f_values, grid_pmf)
  result_fn <- as_piecewise_linear(cumsum(agg$weights), agg$values)
  result <- if (right_tail) get_negated(result_fn) else result_fn
  return(result)
}

#' The CDF of beta^2 over genes as a piecewise linear function.
#' 
#' @param model A model object.
#' @param right_tail A logical value indicating whether to return the right tailed CDF.
#' @return The CDF of beta^2 over genes.
#' @export
get_gene_cdf_betasq <- function(model, right_tail=FALSE) {
  pmf <- gene_PMF(model)
  return(get_piecewise_linear_cdf(grid_values = model$grid_effects, grid_pmf = pmf, f = function(t) t^2, right_tail=right_tail))
}

#' The CDF of beta^2 over components of variance as a piecewise linear function.
#' 
#' @param model A model object.
#' @param right_tail A logical value indicating whether to return the right tailed CDF.
#' @return The CDF of beta^2 over components of variance.
#' @export
get_variance_cdf_betasq <- function(model, right_tail=FALSE) {
  pmf <- variance_PMF(model)
  return(get_piecewise_linear_cdf(
              grid_values = model$grid_effects, 
              grid_pmf = pmf, 
              f = function(t) t^2, 
              right_tail=right_tail)
              )
}

#' The quantile function of beta^2 over genes as a piecewise linear function.
#' 
#' @param model A model object.
#' @param right_tail A logical value indicating whether to return the right tailed QF.
#' @return The quantile function of the distribution of beta^2 over genes.
#' @export
get_gene_qf_betasq <- function(model, right_tail=FALSE) {
  pmf <- gene_PMF(model)
  return(get_piecewise_linear_qf(
            grid_values = model$grid_effects, 
            grid_pmf = pmf, 
            f = function(t) t^2, 
            right_tail=right_tail)
            )
}

#' The quantile function of beta^2 over components of variance as a piecewise linear function.
#' 
#' @param model A model object.
#' @param right_tail A logical value indicating whether to return the right tailed QF.
#' @return The quantile function of the distribution of beta^2 over components of variance.
#' @export
get_variance_qf_betasq <- function(model, right_tail=FALSE) {
  pmf <- variance_PMF(model)
  return(get_piecewise_linear_qf(grid_values = model$grid_effects, grid_pmf = pmf, f = function(t) t^2, right_tail=right_tail))
}


# ---- Derivatives of distribution functions with respect to model mixture weights ----

#' Derivative of gene_PMF with respect to the mixture weight of component i, evaluated at the grid
#' 
#' This is the i-th column of `model$components`. It is constant w.r.t. the mixture weights.
#'
#' @param model A model object.
#' @param component_index The index `i` of the mixture component.
#' @return A numeric vector for the derivative.
#' @export
get_dG_dw <- function(model, component_index, right_tail=FALSE) {
  
  dG_dw <- get_piecewise_linear_cdf(
    grid_values = model$grid_effects, 
    grid_pmf = as.numeric(model$components[component_index, ]), 
    f = function(t) t^2, 
    right_tail = right_tail
  )
}

#' Derivative of variance_PDF with respect to the mixture weight of component i
#'
#' Uses the quotient rule to account for variance normalization.
#'
#' @param component_index The index `i` of the mixture component.
#' @param right_tail A logical value indicating whether to return the right tailed CDF.
#' @return A numeric vector for the derivative.
#' @export
get_dV_dw <- function(model, component_index, right_tail=FALSE) {
  grid_sq <- model$grid_effects^2

  # PMF = numerator/denominator
  numerator <- grid_sq * gene_PMF(model)
  denominator <- sum(numerator)
  d_numerator_dw <- grid_sq * as.numeric(model$components[component_index, ])
  d_denominator_dw <- sum(d_numerator_dw)

  # d/dw PMF
  d_quotient_dw <- (d_numerator_dw * denominator - numerator * d_denominator_dw) / (denominator^2)
  
  return(get_piecewise_linear_cdf(model$grid_effects, d_quotient_dw, f = function(t) t^2, right_tail=right_tail))

}


# ---- Derivatives of distribution functions with respect to their arguments ----

#' Creates a piecewise constant function interpolating between points on a grid.
#' 
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @return A piecewise constant function mapping f(x[i]) -> y[i].
#' @export
as_piecewise_constant <- function(x, y) {
  return(approxfun(x, y, method = "constant", 
                   rule = 1, f = 1, ties = "ordered"))
}

#' The PDF of a transformed random variable, approximated as a piecewise constant function
#' 
#' @param grid_values A numeric vector of x-coordinates.
#' @param grid_pmf The PMF evaluated at the grid values.
#' @param f A transformation to be applied to the grid values.
#' @return The PMF of the transformed random variable.
#' @export
get_piecewise_linear_pdf <- function(grid_values, grid_pmf, f=function(t) t) {
  cdf <- get_piecewise_linear_cdf(grid_values, grid_pmf, f=f)
  epsilon <- 1e-12
  pdf <- function(x) {
    (cdf(x + epsilon) - cdf(x)) / epsilon
  }
  return(pdf)
}

#' The PDF of beta^2 over genes, approximated as a piecewise constant function
#'
#' @param model A model object.
#' @return The PDF of beta^2 over genes.
#' @export
get_gene_pdf_betasq <- function(model) {
  pmf <- gene_PMF(model)
  return(get_piecewise_linear_pdf(grid_values = model$grid_effects, grid_pmf = pmf, f = function(t) t^2))
}

#' The PDF of beta^2 over components of variance, approximated as a piecewise constant function
#'
#' @param model A model object.
#' @return The PDF of beta^2 over components of variance.
#' @export
get_variance_pdf_betasq <- function(model) {
  pmf <- variance_PMF(model)
  return(get_piecewise_linear_pdf(grid_values = model$grid_effects, grid_pmf = pmf, f = function(t) t^2))
}

#' The derivative of the quantile function of a transformed random variable, 
#' approximated as a piecewise constant function
#' 
#' @param grid_values A numeric vector of x-coordinates.
#' @param grid_pmf The PMF evaluated at the grid values.
#' @param f A transformation to be applied to the grid values.
#' @return The PMF of the transformed random variable.
#' @export
get_piecewise_linear_dqf <- function(grid_values, grid_pmf, f=function(t) t) {
  pdf <- get_piecewise_linear_pdf(grid_values = grid_values, grid_pmf = grid_pmf, f = f)
  qf <- get_piecewise_linear_qf(grid_values = grid_values, grid_pmf = grid_pmf, f = f, right_tail = FALSE)

  # d/dx f^-1(x) = 1 / f'(f^-1(x))
  dqf <- function(q_values) {
    quantiles_at_q <- qf(q_values)
    pdf_values_at_quantiles <- pdf(quantiles_at_q)
    return(1 / pdf_values_at_quantiles)
  }
  
  return(dqf)
}

#' The derivative of the quantile function of beta^2 over genes, 
#' approximated as a piecewise constant function
#'
#' @param model A model object.
#' @return The derivative of the quantile function.
#' @export
get_gene_dqf_betasq <- function(model) {
  pmf <- gene_PMF(model)
  return(get_piecewise_linear_dqf(grid_values = model$grid_effects, grid_pmf = pmf, f = function(t) t^2))
}

#' The derivative of the quantile function of beta^2 over components of variance, 
#' approximated as a piecewise constant function
#'
#' @param model A model object.
#' @return The derivative of the quantile function.
#' @export
get_variance_dqf_betasq <- function(model) { 
  pmf <- variance_PMF(model)
  return(get_piecewise_linear_dqf(grid_values = model$grid_effects, grid_pmf = pmf, f = function(t) t^2))
}

#' Calculate the derivative of the variance quantile function Vinv(v) 
#' with respect to a mixture weight w_i.
#' 
#' Formula: (∂Vinv/∂w_i)(v) = [ (∂V/∂w_i)(Vinv(v)) ] / p_V(Vinv(v))
#'
#' @param model The model object.
#' @param component_index The index of the mixture weight.
#' @return A function that takes a quantile `v` and returns (∂Vinv/∂w_i)(v).
get_dVinv_dw <- function(model, component_index) {
  # Right-tailed CDF derivative: cumulative sum from the upper tail
  dV_dw <- get_dV_dw(model, component_index, right_tail = TRUE)
  Vinv <- get_variance_qf_betasq(model, right_tail = TRUE)
  variance_pdf <- get_variance_pdf_betasq(model) # dV/dx

  # The derivative function (∂Vinv/∂w_i)(v)
  dVinv_dw <- function(v) {
    q <- Vinv(v)
    numerator <- dV_dw(q)
    denominator <- variance_pdf(q)
    return(numerator / denominator)
  }
  
  return(dVinv_dw)
}

# ---- Functions for the fraction of genes needed to explain a given fraction of h2 ----

#' A function for the fraction of genes needed to explain a given fraction of h2.
#'
#' @param model A fitted model object.
#' @return A function F(v) that takes a proportion of heritability and returns
#'         the fraction of genes needed to explain that proportion of h2.
#' @export
get_needed_genes_fn <- function(model) {
  V_inv_fn <- get_variance_qf_betasq(model, right_tail=TRUE)
  G_fn <- get_gene_cdf_betasq(model, right_tail=TRUE)

  F <- function(v) {
    actual_variance_amount <- V_inv_fn(v)
    gene_cdf_proportion <- G_fn(actual_variance_amount)
    return(gene_cdf_proportion)
  }
  return(F)
}

#' The derivative of needed_genes_fn with respect to the mixture weights
#'
#' @param model A fitted model object.
#' @param component_index The index of the component to differentiate with respect to.
#' @return A function F(v) that takes a proportion of heritability and returns
#'         the fraction of genes needed to explain that proportion of h2.
#' @export
get_needed_genes_fn_derivative <- function(model, component_index) {
  
  # (∂G/∂w)(x) 
  dG_dw <- get_dG_dw(model, component_index, right_tail=TRUE)
  
  gene_pdf <- get_gene_pdf_betasq(model)
  Vinv <- get_variance_qf_betasq(model, right_tail=TRUE)
  
  # (∂Vinv/∂w)(v)
  dVinv_dw <- get_dVinv_dw(model, component_index)

  # dF/dw|v = (∂G/∂w)(Vinv(v)) - gene_pdf(Vinv(v)) * (∂Vinv/∂w)(v)
  dF_dw_fn <- function(v) {
    q <- Vinv(v)
    return(dG_dw(q) - gene_pdf(q) * dVinv_dw(v))
  }
  
  return(dF_dw_fn)
}

get_cov_w_matmat <- function(model) {
  feature_weights <- colMeans(do.call(rbind, model$df$features))

  matvec <- function(x) {
    result <- matrix(0, nrow = nrow(model$components) - 1, ncol = ncol(x))
    for (i in 1:length(feature_weights)) {
      offset <- 1e-6 * diag(diag(model$information[[i]]))#mean(diag(model$information[[i]])) * diag(nrow(model$information[[1]]))
      # print(offset)
      # print(model$information[[i]])
      product <- solve(model$information[[i]] + offset, x)
      result <- result + feature_weights[i] * product
      # print('x:')
      # print(x)
      # print('product:')
      # print(product)

    }
    return(result)
  }
}

#' Calculate quantiles of genes and variance with standard errors
#'
#' @param model A model object, which must contain:
#'   - `delta`: Matrix of mixture weights (features x components).
#'   - `null_index`: Index of the null component.
#'   - `information`: A list of information matrices, one for each feature category.
#'                    Each matrix corresponds to the non-null components.
#'   And be compatible with `get_needed_genes_fn` and `get_needed_genes_fn_derivative`.
#'
#' @return A function which takes a vector of variance quantiles and 
#' returns a dataframe with the following columns:
#'   - `variance`: The input vector
#'   - `needed_genes`: The fraction of genes needed to explain V[i] of variance
#'   - `se`: The standard errors of the quantiles
#'   - `dVdG`: The nth largest contribution to variance, where n=needed_genes[i]*num_genes;
#'            i.e., the variance quantile function evaluated at V[i]
#'   - variance_se: Standard error of the fraction of variance explained by contributions >= dVdG[i]
#'   - needed_genes_se: Standard error of the fraction of genes with contributions >= dVdG[i]
#' @export
get_quantiles_with_se_fn <- function(model) {

  quantiles_with_se_fn <- function(v) {
    num_genes <- nrow(model$df)
    needed_genes_fn <- get_needed_genes_fn(model)
    needed_genes_var_fn <- get_needed_genes_var_fn(model)
    Vinv_fn <- get_variance_qf_betasq(model, right_tail=TRUE)
    df <- data.frame(
      variance = v,
      needed_genes = needed_genes_fn(v),
      dVdG = Vinv_fn(v),
      needed_genes_se = sqrt(needed_genes_var_fn(v)) / num_genes
    )

    total_variance <- sum(posterior_expectation2(model, function_to_integrate = function(x) {x^2}))
    for (i in 1:length(v)) {
      q <- Vinv_fn(v[i])
      f <- function(x) {x^2 > q}
      result <- posterior_expectation_with_se(model, function_to_integrate = f)
      df$genes_cdf_se[i] <- sqrt(sum(result$se^2))  / num_genes

      f <- function(x, row) {model$h2_function(x, row) * (x^2 > q)}
      result <- posterior_expectation_with_se(model, function_to_integrate = f)
      df$variance_se[i] <- sqrt(sum(result$se^2))  / total_variance
    }
    return(df)
    
  }
  return(quantiles_with_se_fn)
  }

#' Calculate Point Estimates and Standard Errors for F(v)
#'
#' Computes F(v) (the fraction of genes needed to explain a proportion v of heritability)
#' and its standard error using the delta method. The standard error is derived from
#' the derivatives of F(v) with respect to mixture weights and the information matrix
#' from the model.
#'
#' @param model A model object, which must contain:
#'   - `delta`: Matrix of mixture weights (features x components).
#'   - `null_index`: Index of the null component.
#'   - `information`: A list of information matrices, one for each feature category.
#'                    Each matrix corresponds to the non-null components.
#'   And be compatible with `get_needed_genes_fn` and `get_needed_genes_fn_derivative`.
#' @param v_values A numeric vector of proportions of heritability (0 to 1) at which
#'   to evaluate F(v) and its standard error.
#'
#' @return A list with two named elements:
#'   - `point_estimates`: A numeric vector of F(v) values corresponding to `v_values`.
#'   - `standard_errors`: A numeric vector of standard errors for F(v) corresponding to `v_values`.
#'
#' @export
get_needed_genes_var_fn <- function(model) {

  needed_genes_var_fn <- function(v) {
    cov_w_matmat <- get_cov_w_matmat(model)
    Z <- matrix(NA, nrow = length(v), ncol = nrow(model$components))

    for (j in 1:nrow(model$components)) {
      dF_dw_j_fn <- get_needed_genes_fn_derivative(model, component_index = j)
      Z[, j] <- dF_dw_j_fn(v)
    }
    Z <- Z[, -model$null_index] - Z[, model$null_index] %*% matrix(1, nrow = 1, ncol = nrow(model$components)-1)
    result <- rowSums(Z * t(cov_w_matmat(t(Z))))
    print(result)
    return(result)
  }
  return(needed_genes_var_fn)
  }

#' Calculate the cumulative variance function, defined as the fraction of variance explained
#' by the top fraction p of genes.
#'
#' @param model A model object.
#' @return A function that takes a fraction of genes and returns the cumulative variance explained by the top p genes.
#' @export
get_cumulative_variance_fn <- function(model) {
  gene_upper_quantile_fn <- get_gene_qf_betasq(model, right_tail=TRUE)
  gene_upper_cdf <- get_gene_cdf_betasq(model, right_tail=TRUE)
  variance_upper_cdf <- get_variance_cdf_betasq(model, right_tail=TRUE)

  cumulative_variance_fn <- function(p) {
    v <- gene_upper_quantile_fn(p)
    V <- variance_upper_cdf(v)
    return(V)
  }
}