### Threshold functions
usethis::use_package('CompQuadForm')
### Find thresholds for quadratic testing
thresholdWrapperQuadratic <- function(y,
                                      K,
                                      cov.mat,
                                      first.thres,
                                      second.thres = NULL,
                                      t1 = NULL,
                                      t2 = NULL,
                                      alpha) {
  if (is.null(t1)) {
    t1 <- getProbThresholdQuadratic(quant = first.thres, cov.mat = cov.mat, K = K)
  }
  if (is.null(second.thres)) {
    first.thres  <- getThresholdQuadratic(prob = 1 - t1, cov.mat = cov.mat, K = K)
  }
  if (is.null(t2)) {
    t2 <- alpha^2 * t1
  }
  second.thres <- getThresholdQuadratic(prob = 1 - t2, cov.mat = cov.mat, K = K)
  return(list('first.threshold'  = first.thres,
              'second.threshold' = second.thres,
              't1' = t1,
              't2' = t2))
}

### Find threshold for linear testing
thresholdWrapperLinear <- function(y,
                                   a.vec,
                                   cov.mat,
                                   first.threshold.lower,
                                   first.threshold.upper,
                                   second.thres = NULL,
                                   t1.lower = NULL,
                                   t1.upper = NULL,
                                   t2.lower = NULL,
                                   t2.upper = NULL,
                                   alpha) {
  sd.test <- sqrt(t(a.vec) %*% cov.mat %*% a.vec)
  if (is.null(t1)) {
    t1.lower <- pnorm(first.threshold.lower, mean = 0, sd = sd.test)
    t1.upper <- pnorm(first.threshold.upper, mean = 0, sd = sd.test, lower.tail = FALSE)
  }
  if (is.null(second.thres)) {
    first.threshold.lower <- qnorm(prob = 1 - t1.lower, mean = 0, sd = sd.test)
    first.threshold.upper <- qnorm(prob = 1 - t1.upper, mean = 0, sd = sd.test)
  }
  if (is.null(t2)) {
    t2.lower <- (alpha / 2)^2 * t1.lower
    t2.upper <- (alpha / 2)^2 * t1.upper
  }
  second.threshold.lower <- qnorm(prob = 1 - t2.lower, mean = 0, sd = sd.test)
  second.threshold.upper <- qnorm(prob = 1 - t2.upper, mean = 0, sd = sd.test)
  return(list('first.threshold.lower'  = first.threshold.lower,
              'first.threshold.upper' = second.thres,
              't1.lower' = t1.lower,
              't1.upper' = t1.upper,
              't2.lower' = t2.lower,
              't2.upper' = t2.upper))
}
### Finds the qunailte distribution of weighted chi square, where lambda are the eigen values
### of the M matrix in X'MX
findQuantileWeightedChiLiu <- function(prob, lambda, jump.size = 2) {
  ### Finding limits for uniroot
  lower.search <- 0
  upper.search <- 1
  prob.upper.search  <- (1 - CompQuadForm::liu(upper.search, lambda))
  ## We want the upper bound to have probability larger than p
  ## assures us that upper.bound > q
  while (prob.upper.search < prob) {
    lower.search       <- upper.search ## Assures us that lower.bound < q
    upper.search       <- upper.search * jump.size ## Increase upper search bound
    prob.upper.search  <- (1 - CompQuadForm::liu(upper.search, lambda)) ## Update probability value
  }
  return(uniroot(function(q) (1 - CompQuadForm::liu(q, lambda)) - prob,
                 lower = lower.search,
                 upper = upper.search)$root)
}


### Find the eigen values of test matrix

#' Find Delta of test matrix
#'
#' Function that finds the delta of the test matrix
#'
#' @param cov.mat The covariance matrix of \code{y}.
#' @param K The test matrix \eqn{K} used in the aggregate test
#'
#' @details Finding the weighted chi-square eigen values based on the following
#' \eqn{y ~ N(\mu, \Sigma)}, and we have the following quadratic form
#' \eqn{y' K y = y' \Sigma^{-1/2} \Sigma^{1/2} K \Sigma^{1/2} \Sigma^{-1/2} y}.
#' The functions find the eigen-values of \eqn{\Sigma^{1/2} K \Sigma^{1/2}, where $\Sigma^{1/2}}
#' is found using cholesky decomposition.

getTestMatEigen <- function(cov.mat, K) {
  root.cov.mat <- chol(cov.mat)
  return(svd(root.cov.mat %*% K %*% t(root.cov.mat))$d)
}

### Finds threshold of quadratic test given a probability
getThresholdQuadratic <- function(prob, cov.mat, K, jump.size = 2) {
  lambda <- getTestMatEigen(cov.mat, K)
  return(findQuantileWeightedChiLiu(prob, lambda = lambda, jump.size = jump.size))
}

### Finds the probability given a threshold
getProbThresholdQuadratic <- function(quant, cov.mat, K, jump.size = 2) {
  lambda <- getTestMatEigen(cov.mat, K)
  return(1 - CompQuadForm::liu(quant, lambda = lambda))
}


