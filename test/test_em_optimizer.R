set.seed(123)

n_genes <- 100
n_components <- 5
n_features <- 3

cdl <- matrix(runif(n_genes * n_components), nrow = n_genes, ncol = n_components)
cdl <- cdl / rowSums(cdl)

features <- matrix(0, nrow = n_genes, ncol = n_features)
for (i in 1:n_genes) {
  features[i, sample(1:n_features, 1)] <- 1
}

delta_init <- matrix(runif(n_features * n_components), nrow = n_features, ncol = n_components)
delta_init <- delta_init / rowSums(delta_init)

source("R/EM.R")

result1 <- EM_optimize(cdl, features, delta_init, max_iter = 100)
result2 <- EM_optimize(cdl, features, delta_init, max_iter = 100)

cat("Delta difference:", max(abs(result1$delta - result2$delta)), "\n")
cat("Log-likelihood 1:", result1$ll, "\n")
cat("Log-likelihood 2:", result2$ll, "\n")
cat("LL difference:", abs(result1$ll - result2$ll), "\n")

if (max(abs(result1$delta - result2$delta)) < 1e-10 && abs(result1$ll - result2$ll) < 1e-10) {
  cat("Test PASSED: Both runs give identical results\n")
} else {
  cat("Test FAILED: Results differ\n")
}
