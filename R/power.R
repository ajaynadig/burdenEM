# Power test
NTPR <- function(bb, pp, nn, alpha) {
  # NTPR computes the not-by-chance true positive rate for a mixture of uniform
  # distributions with parameters [0,b] for b in bb, or [b,0] if b<0, and
  # mixture weights p in pp, at sample sizes nn and significance thresholds
  # alpha. It also computes h2gwas, the fraction of h2 that will be
  # explained by significant SNPs/genes at that significance threshold, and
  # h2power, the expected fraction of SNPs/genes that will be significant.

  # Ensure input dimensions are correct
  stopifnot(is.vector(nn))
  stopifnot(is.vector(alpha))
  stopifnot(is.vector(bb))
  stopifnot(is.vector(pp))

  chisq_threshold <- qchisq(1 - alpha, df = 1)

  # Power function
  powerfn <- function(x, a, n) {
    pnorm(sqrt(n) * x - sqrt(a)) + pnorm(-sqrt(n) * x - sqrt(a))
  }

  # Not-by-chance true positive rate function
  ntprfn <- function(x, a, n) {
    (pnorm(sqrt(n) * abs(x) - sqrt(a)) - pnorm(-sqrt(n) * abs(x) - sqrt(a))) / powerfn(x, a, n)
  }

  # Initialize total power, NTPR, h2gwas, positive and negative power components
  total_power <- matrix(0, nrow = length(alpha), ncol = length(nn))
  total_ntpr <- matrix(0, nrow = length(alpha), ncol = length(nn))
  total_h2gwas <- matrix(0, nrow = length(alpha), ncol = length(nn))
  total_power_pos <- matrix(0, nrow = length(alpha), ncol = length(nn))
  total_power_neg <- matrix(0, nrow = length(alpha), ncol = length(nn))

  for (ii in seq_along(pp)) {
    xx <- seq(0, 1, by = 0.001) * bb[ii]

    pow <- array(NA, dim = c(length(alpha), length(nn), length(xx)))
    for (i in seq_along(alpha)) {
      for (j in seq_along(nn)) {
        pow[i, j, ] <- powerfn(xx, chisq_threshold[i], nn[j])
      }
    }
    cpt_power <- apply(pow, c(1, 2), mean)

    ntpr <- array(NA, dim = c(length(alpha), length(nn), length(xx)))
    for (i in seq_along(alpha)) {
      for (j in seq_along(nn)) {
        ntpr[i, j, ] <- ntprfn(xx, chisq_threshold[i], nn[j])
      }
    }
    cpt_ntpr <- apply(pow * ntpr, c(1, 2), mean)


    h2gwas <- xx^2
    h2gwas_array <- array(h2gwas, dim = c(length(alpha), length(nn), length(xx))) # Ensure dimensions match
    print(dim(pow))
    print(dim(h2gwas_array))

    # Now perform element-wise multiplication after ensuring conformity
    cpt_h2gwas <- apply(pow * h2gwas_array, c(1, 2), mean)

    total_power <- total_power + pp[ii] * cpt_power
    total_ntpr <- total_ntpr + pp[ii] * cpt_ntpr
    total_h2gwas <- total_h2gwas + pp[ii] * cpt_h2gwas

    if (bb[ii] > 0) {
      total_power_pos <- total_power_pos + pp[ii] * cpt_power
    } else if (bb[ii] < 0) {
      total_power_neg <- total_power_neg + pp[ii] * cpt_power
    }
  }

  total_ntpr <- total_ntpr / total_power
  total_h2gwas <- total_h2gwas / sum(pp * bb^2 / 3)

  return(list(total_ntpr = total_ntpr, total_h2gwas = total_h2gwas,
              total_power = total_power, total_power_pos = total_power_pos,
              total_power_neg = total_power_neg))
}
