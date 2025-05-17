// File: EM_fit.cpp

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List EM_fit_cpp(List model,
                int max_iter,
                double tol = 1e-6,
                bool return_likelihood = true,
                bool return_numiter = true) {

  // 1) Extract everything from the R list
  arma::mat X  = as<arma::mat>(model["features"]);               // N × P
  arma::mat CL = as<arma::mat>(model["conditional_likelihood"]); // N × K
  arma::mat δ  = as<arma::mat>(model["delta"]);                  // P × K

  int N = X.n_rows;

  // 2) Pre-compute (X'X)^{-1} X'
  arma::mat OLS_den    = arma::inv( X.t() * X );   // P × P
  arma::mat OLS_t_feat = OLS_den * X.t();          // P × N

  // 3) Prepare storage for log-likelihood
  arma::vec ll(max_iter, arma::fill::none);
  ll.fill(arma::datum::nan);

  double ll_change = tol + 1.0;

  // --- first iteration ---
  arma::mat W = X * δ;                // N × K   : weights
  arma::mat post = W % CL;            // N × K   : elementwise multiply
  arma::vec rowsum = arma::sum(post,1);       // N × 1
  ll(0) = arma::accu(arma::log(rowsum));      // scalar
  post.each_col() /= rowsum;                 // normalize each row
  δ = OLS_t_feat * post;                     // P × K

  int iter = 1;

  // --- subsequent iterations ---
  while (iter < max_iter && ll_change > tol) {
    W    = X * δ;
    post = W % CL;
    rowsum = arma::sum(post,1);

    double ll_new = arma::accu(arma::log(rowsum));
    ll(iter) = ll_new;
    ll_change = std::abs(ll_new - ll(iter-1)) / std::abs(ll(iter-1));

    post.each_col() /= rowsum;
    δ = OLS_t_feat * post;

    ++iter;
  }

  // ————————————————
  // print actual iteration count
  if (return_numiter) {
  Rcpp::Rcout << "EM_fit_cpp: converged in " << iter << " iterations." << std::endl;
  }
  // ————————————————

  // 4) Pack results back into the model list
  model["delta"] = δ;
  if (return_likelihood) {
    // convert arma::vec → R numeric vector
    NumericVector Rll(max_iter);
    for (int i = 0; i < max_iter; ++i) Rll[i] = ll(i);
    model["ll"] = Rll;
  }
  return model;
}
