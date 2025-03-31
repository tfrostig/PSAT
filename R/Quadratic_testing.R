usethis::use_package('devtools')
usethis::use_package('DEoptim')
usethis::use_package('tmg')
usethis::use_package('truncnorm')
usethis::use_package('MASS')
# Quadratic ---------------------------------------------------------------
#' @title Conduct Post Selection After Qudratic Aggregation Testing inference
#'
#' @description \code{PSATQuadratic} is used to estimate (point & interval) and inferon normal means model that was selected based
#' on a single quadratic aggregate test of the form:
#' \deqn{y' K y > c > 0}
#'
#' @param y The observed normal vector.
#' @param cov.mat The covariance matrix of \code{y}.
#' @param K The test matrix \eqn{K} used in the aggregate test. Defaults to \code{solve(cov.mat)}, WALD test.
#' @param threshold The threshold \eqn{c > 0} used in the aggregate test.
#' @param contrast An optional matrix of contrasts to be tested: must have number of columns
#' identical to the length of \code{y}. If left as \code{NULL}, the coordinates of \code{y} will be tested by default.
#' @param alpha Indicates the type I error level in which to find rejection areas and confidence intervals.
#' @param alternative A string.vector indicated the for which contrast what type of test to conduct
#' \code{'two.sided'}, \code{'less'} or \code{'greater'}. Should have the same length as number of contrasts.
#' Defaults to two.sided tests.
#' @param test.type Type of tests (UMPU and Symmetric), only relevant for two-sided testing.
#' @param pval.type  Type of p-values to calculate for each contrast.
#' Options available are \code{'global', 'hybrid', 'polyhedral', 'naive'}.
#' @param ci.type Type of confidence intervals to create for each contrast.
#'  Options available are \code{'global', 'hybrid', 'polyhedral', 'switch', 'naive'}.
#' @param est.type Type of point estimation, options available are \code{'moment', 'naive'}.
#' @param verbose Whether to report on the progress of the computation.
#' @param zscore.bound Number of the furthest bound before ignoring truncation and addressing the distribution as normal.
#' Defaults to the \code{pval.type}.
#' @param samp.para List of parameters for the tmg truncated normal sampler.
#' Contains two parameters \code{samp.size} and \code{burn.in}, for more details see \code{tmg::rtmg}.
#' @param optim.control allows the user to set some characteristics of the Differential Evolution optimization algorithm implemented in DEoptim.
#' For more details see code{\link[DEoptim]{DEoptim.control}}.
#'
#' @details The function is used to perform inference for normal mean vectors
#' that were selected based on a single quadratic aggregate test. To be exact, suppose
#' that \eqn{y ~ N(\mu,\Sigma)} and that we are interested in estimating \eqn{\mu}
#' only if we can determine that \eqn{\mu\neq 0} using an aggregate test of the form:
#' \deqn{y' K y > c > 0} for some predetermined constant \eqn{c}.
#'
#' \code{PSATQuadratic} offers several options for computing p-values. The "global-null"
#' method relies on comparing the magnitude of \eqn{y} to samples from the truncated
#' global-null distribution. This method is powerful when \eqn{\mu} is sparse and its
#' non-zero coordinates are not very large. The `polyhedral` method is exact when the
#' observed data is approximately normal and is quite robust to model misspecification.
#' It tends to be more powerful than the `global-null` method when the magnitude of
#' \eqn{\mu} is large. The "hybrid" method combines the strengths of the "global-null"
#' and `polyhedral` methods, possessing good power regardless of the sparsity or
#' magnitude of \eqn{\mu}. However it is less robust to the misspecification of the distribution
#' of \eqn{y} than the "polyhedral" method. The confidence interval methods are similar to the p-values ones.
#' For more details see \href{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12318}{Post-selection estimation and testing following aggregate association tests}.
#' @examples
#' \dontrun{
#' y         <- c(15, rnorm(p - 1))
#' cov.mat <- diag(p)
#' PSATQuadratic(y, cov.mat) ### Runs on default values
#' }
#'
#' \dontrun{
#' y         <- c(10, rnorm(p - 1))
#' contrast  <- rbind(rep(1, p), c(1, 1, 1, rep(0, p - 3)))
#' PSATQuadratic(y = y, cov.mat = cov.mat,
#'               contrast = contrast, alternative =  c('less', 'greater')) ## Test different alternatives
#' }
#' @export
PSATQuadratic <- function(y,
                          cov.mat,
                          alpha = 0.05,
                          K = solve(cov.mat),
                          threshold = qchisq(1 - alpha, nrow(cov.mat)),
                          contrast = NULL,
                          alternative = NULL,
                          test.type = c('symmetric', 'UMPU'),
                          pval.type = c('switch', 'global', 'hybrid', 'polyhedral', 'naive'),
                          ci.type   = pval.type,
                          est.type  = c('naive', 'moment'),
                          verbose   = FALSE,
                          zscore.bound = 5,
                          samp.para = list(samp.size = 10000,
                                           burn.in   = 200),
                          optim.control = DEoptim::DEoptim.control(),
                          optim.method = 'Nelder-Mead') {
  ### Initalizing paramters
  p         <- length(y)
  test.val  <- drop(t(y) %*% K %*% y)
  optim.control$trace <- verbose
  pval.type           <- unique(c(ci.type, pval.type))
  # Tests  ------------------------------------------------------------------
  if (test.val < threshold) {
    stop('Test static below threshold, y is not tested. Returns data.frame of NA')
  }
  ### Eigen decomposition
  eigen.decomp <- eigen(cov.mat)
  if (any(eigen.decomp$values <= 0)) {
    stop('Covariance matrix is not invertiable')
  }
  if (is.null(contrast)) {
    contrast  <- diag(p)
  }
  if (is.null(alternative)) {
    alternative <- rep('two.sided', nrow(contrast))
  }
  if (!any(unique(alternative) %in%  c('two.sided', 'less', 'greater'))) {
    stop('Testing sided are ill defined. \n
         alternative accepts only the following: two.sided, less and greater')
  }
  if (nrow(contrast) != length(alternative)) {
    stop('Must specifiy which side to test for each contrast')
  }
  ### Thresholding
  thres.list <- thresholdWrapperQuadratic(test.val = test.val,
                                          K = K,
                                          cov.mat = cov.mat,
                                          first.threshold = threshold,
                                          second.threshold = NULL, ## Should add to parameters?
                                          alpha = alpha,
                                          calc.t2 = 'switch' %in% ci.type)
  ### Finding the contrast values
  contrast.y     <- drop(contrast %*% y)
  var.contrast.y <- diag(contrast %*% cov.mat %*% t(contrast))
  ### Preparing data according to requested types
  trunc.para.mat <- apply(contrast, 1, findTruncParametersQuardratic,
                          cov.mat = cov.mat, K = K, threshold = threshold, y = y)
  trunc.para.mat <- cbind(t(trunc.para.mat), 'alternative' = sideToNum(alternative))
  if (any(c(pval.type, ci.type) %in% c('hybrid', 'global.null', 'switch'))) {
    samp.null         <- sampleGNQuadratic(y = y,
                                           null.mu = rep(0, p),
                                           sigma   = cov.mat,
                                           threshold = threshold,
                                           K.mat = K,
                                           num.samp = samp.para$samp.size,
                                           burn.in  = samp.para$burn.in)
  }
  # Pvalues -----------------------------------------------------------------
  pval.df <- pvalueWrapper(truncation.parameters = trunc.para.mat,
                           sample.global.null    = samp.null,
                           contrast              = contrast,
                           alpha                 = alpha,
                           pval.type             = pval.type,
                           test.type             = test.type,
                           optim.control         = optim.control,
                           zscore.bound          = zscore.bound)
  # Rejection area  ---------------------------------------------------------
  # reject.df <- rejectionAreaWrapper(truncation.parameters = trunc.para.mat,
  #                                   sample.global.null    = samp.null,
  #                                   contrast.y            = contrast.y,
  #                                   var.contrast.y        = var.contrast.y,
  #                                   contrast              = contrast,
  #                                   alpha                 = alpha,
  #                                   test.type             = test.type,
  #                                   rejection.area.type   = rej.type,
  #                                   optim.control         = optim.control)
  # CIs ---------------------------------------------------------------------
  ci.df <- CIwrapper(truncation.parameters = trunc.para.mat,
                     sample.global.null    = samp.null,
                     contrast              = contrast,
                     alpha                 = alpha,
                     test.type             = test.type,
                     ci.type               = ci.type,
                     thres.list            = thres.list,
                     zscore.bound          = zscore.bound)

  # Point estimation  -------------------------------------------------------
  est.df <- estimationWrapper(truncation.parameters = trunc.para.mat,
                              y                     = y,
                              zscore.bound          = zscore.bound,
                              test.mat              = K,
                              cov.mat               = cov.mat,
                              threshold             = threshold,
                              optim.method          = optim.method,
                              est.type              = est.type)
  return(list('Pvalues' = pval.df,
              ## 'Rejection' = reject.df,
              'CI' = ci.df,
              'Point.estimation'  = est.df,
              'Threshold.details' = thres.list))
}

### Finding the parameters of the truncated normal
findTruncParametersQuardratic <- function(y, cov.mat, K, threshold, eta) {
  p         <- length(y)
  I.p       <- diag(p)
  etaSigeta <- drop(t(eta) %*% cov.mat %*% eta)
  c.const   <- (1 / etaSigeta) * (cov.mat %*% eta)
  ## The effect orthogonal to the linear contrast
  W         <- (I.p - c.const %*% t(eta)) %*% y
  ### Matrix multiplication for later
  WKC     <- drop(t(W) %*% K %*% c.const)
  CKC     <- drop(t(c.const) %*% K %*% c.const)
  WKW     <- drop(t(W) %*% K %*% W)
  ### Finding delta
  delta       <- 4 * (WKC^2 - CKC * (WKW - threshold))
  trunc.bound <- findABQuad(WKC, CKC, delta)
  return(c('x' = t(eta) %*% y, 'var' = etaSigeta, 'bounds' = trunc.bound))
}


### Find A(W) and B(W)
findABQuad <- function(WKC, CKC, delta) {
  if (delta >= 0) {
    A.W <- (-2 * WKC + sqrt(delta)) / (2 * CKC)
    B.W <- (-2 * WKC - sqrt(delta)) / (2 * CKC)
  } else {
    A.W <-  Inf
    B.W <- -Inf
  }
  return(c('lower' = ifelse(is.nan(B.W), -Inf, B.W), ### NaN happens if CKC 0, assuming K SPD then -inf
           'upper' = ifelse(is.nan(A.W), Inf, A.W)))
}


### Sampling from global null
sampleGNQuadratic <- function(y, null.mu, sigma, threshold, K.mat,
                              num.samp = samp.para$samp.size, burn.in = 100) {
  p               <- length(null.mu)
  ### For some reason constraint.list needs to be a list of list, otherwise rtmg returns error
  constraint.list <- list(a = list(A = K.mat, B = rep(0, p), C = -threshold))
  omega           <- solve(sigma)
  return(tmg::rtmg(n = num.samp,
                   M = omega,
                   r = drop(omega %*% null.mu),
                   initial = y,
                   q = constraint.list,
                   burn.in = burn.in))
}


### Dealing with sides
numToSide <- function(x) {
  x[x == 0] <- 'two.sided'
  x[x == 1] <- 'greater'
  x[x == 2] <- 'less'
  return(as.character(x))
}

sideToNum <- function(x) {
  x[x == 'two.sided'] <- 0
  x[x == 'less']      <- 1
  x[x == 'greater']   <- 2
  if (any(is.na(x))) {
    stop('A test was not correctly defined')
  }
  return(as.numeric(x))
}


# Rejection area  ---------------------------------------------------------
### Finding rejection area
polyReject <- function(trunc.parameters, test.type = 'symmetric',
                       optim.control, alpha = 0.05,
                       zscore.bound) {
  x     <- trunc.parameters[1]
  sd    <- sqrt(trunc.parameters[2])
  lower <- trunc.parameters[3]
  upper <- trunc.parameters[4]
  ### Condition
  minimal.zscore <- (sign(lower) == sign(upper)) * min(abs(c(lower, upper))) / sd
  ### No selection effect
  if (lower == -Inf | upper == Inf | minimal.zscore > zscore.bound) {
    ### Rejection area
    rej.area.symm <- c('lower' = qnorm(alpha / 2, mean = 0, sd = sd),
                       'upper' = qnorm(1 - alpha / 2, mean = 0, sd = sd))
    return(rej.area.symm)
  }
  ### Selection effect, symmetric p-values
  if (test.type == 'symmetric') {
    rej.area.trunc.symm  <- c(findTruncQuantile(c(alpha / 2, 1 - alpha / 2),
                                                lower.trunc = lower.upper,
                                                upper.trunc = upper.trunc,
                                                mean = 0, sd = sd))
    return(rej.area.trunc.symm)
  }
  ### Selection effect, UMPU p-values
  if (test.type == 'UMPU') {
    rej.area.umpu <- findUMPU(test.stat = NA,
                              alpha = alpha,
                              lower.bound = lower,
                              upper.bound = upper,
                              mean = 0,
                              sd = sd,
                              optim.type = 'rejection',
                              optim.control = optim.control)
    rej.area.umpu <- rej.area.umpu$optim$bestmem
    return(rej.area.umpu)
  }
}

## Rejection area for GN test (Section 3.2 PSAT)
GNReject <- function(sample.mat, alpha) {
  apply(sample.mat, 2, function(x) quantile(x, probs = c('lower' = alpha / 2,
                                                         'upper' = 1 - alpha / 2)))
}

