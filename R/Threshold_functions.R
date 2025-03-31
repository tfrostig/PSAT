### Threshold functions
usethis::use_package('CompQuadForm')
### Find thresholds for quadratic testing
thresholdWrapperQuadratic <- function(test.val,
                                      K,
                                      cov.mat,
                                      first.threshold,
                                      second.threshold = NULL,
                                      t1 = NULL,
                                      alpha,
                                      calc.t2) {
  if (is.null(t1)) {
    t1 <- 1 - getProbThresholdQuadratic(quant = first.threshold, cov.mat = cov.mat, K = K)
  }
  if (is.null(first.threshold)) {
    first.threshold  <- getThresholdQuadratic(prob = 1 - t1, cov.mat = cov.mat, K = K)
  }
  if (calc.t2) {
    t2 <- alpha^2 * t1
    second.threshold      <- getThresholdQuadratic(prob = 1 - t2, cov.mat = cov.mat, K = K)
  }
  if (!calc.t2) {
    second.threshold <- NULL
    t2   <- NULL
  }
  pass.second.threshold <- test.val > second.threshold
  return(list('first.threshold'  = first.threshold,
              'second.threshold' = second.threshold,
              't1' = t1,
              't2' = t2,
              'new.alpha' = alpha - t2 / t1,
              'pass.second.threshold' = pass.second.threshold))
}

### Find threshold for linear testing
thresholdWrapperLinear <- function(test.val,
                                   a.vec,
                                   cov.mat,
                                   first.threshold.lower,
                                   first.threshold.upper,
                                   t1.lower = NULL,
                                   t2.upper = NULL,
                                   alpha,
                                   calc.t2) {
  sd.test  <- sqrt(t(a.vec) %*% cov.mat %*% a.vec)
  if (is.null(t1.lower) & is.null(t2.upper)) {
    t1.lower <- pnorm(first.threshold.lower, mean = 0, sd = sd.test)
    t1.upper <- pnorm(first.threshold.upper, mean = 0, sd = sd.test, lower.tail = FALSE)
  }
  if (is.null(first.threshold.lower) & is.null(first.threshold.lower)) {
    first.threshold.lower <- qnorm(p = t1.lower, mean = 0, sd = sd.test)
    first.threshold.upper <- qnorm(p = 1 - t1.upper, mean = 0, sd = sd.test)
  }
  ### Finding the second threshold parameters
  if (calc.t2) {
    t2.lower <- (alpha)^2 * t1.lower
    t2.upper <- (alpha)^2 * t1.upper
    second.threshold.lower <- qnorm(p = t2.lower, mean = 0, sd = sd.test)
    second.threshold.upper <- qnorm(p = 1 - t2.upper, mean = 0, sd = sd.test)
  }
  if (!calc.t2) {
    t2.lower <- NULL
    t2.upper <- NULL
    second.threshold.lower <- NULL
    second.threshold.upper <- NULL
  }
  pass.second.threshold  <- test.val < second.threshold.lower | test.val > second.threshold.upper
  return(list('first.threshold.lower' = first.threshold.lower,
              'first.threshold.upper' = first.threshold.upper,
              'second.threshold.lower' = second.threshold.lower,
              'second.threshold.upper' = second.threshold.upper,
              't1.lower' = t1.lower,
              't1.upper' = t1.upper,
              't2.lower' = t2.lower,
              't2.upper' = t2.upper,
              'new.alpha' = (t2.lower + t2.upper) / (t1.lower + t1.upper),
              'pass.second.threshold' = pass.second.threshold))
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

#' @title Find quantile for weighted Chi-square distribution
#'
#' @description \code{getThresholdQuadratic} is used to find a quantile specified of the following test statsitic distribution
#' \deqn{y' K y \; y \sim N(y, \Sigma)}
#' @param prob probability
#' @param cov.mat The covariance matrix of \code{y}.
#' @param K The test matrix \eqn{K} used in the aggregate test. Defaults to \code{solve(cov.mat)}, WALD test.
#' @param jump.size parameter for searching the quantile (how fast the increase the searchs)
#' @export
getThresholdQuadratic <- function(prob, cov.mat, K, jump.size = 2) {
  lambda <- getTestMatEigen(cov.mat, K)
  return(findQuantileWeightedChiLiu(prob, lambda = lambda, jump.size = jump.size))
}

### Finds the probability given a threshold
getProbThresholdQuadratic <- function(quant, cov.mat, K, jump.size = 2) {
  lambda <- getTestMatEigen(cov.mat, K)
  return(1 - CompQuadForm::liu(quant, lambda = lambda))
}

getQudraticLam <- function(test.mat, cov.mat, justlam = TRUE) {
  c_ <- chol(cov.mat)
  tempmat <- c_ %*% test.mat %*% t(c_)

  eig <- svd(tempmat)
  vec <- eig$v
  P <- t(vec)
  lam <- eig$d
  if(justlam) {
    return(lam)
  } else {
    deltamat <- P %*% solve(t(c_))
    return(list(lam = lam, deltamat = deltamat))
  }
}
