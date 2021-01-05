#' @title Conduct Post Selection After qudratic Aggregation Testing inference
#'
#' Inference for Normal Means after Aggregate Testing
#'
#' @description \code{PSATLinear} is used to estimate a normal means model that was selected based
#' on a single linear aggregate test of the form:
#' \deqn{a'y > u or a'y < l,}
#'
#' @param y The observed normal vector.
#' @param cov.mat The covariance matrix of \code{y}.
#' @param a.vec The test vector \eqn{a} of size \code{length(y)} used in the aggregate test, defaults to sum of all the test statistics.
#' @param threshold.lower The threshold \eqn{l < 0} used in the aggregate test.
#' @param threshold.upper The threshold  \eqn{u > 0} use in the aggregate test.
#'  To specify a one-sided hypothesis testing, define the appropriate boundary as -Inf or Inf.
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
#' Defaults to the \code{pval.type}.
#' @param zscore.bound Number of the furthest bound before ignoring truncation and addressing the distribution as normal.
#' Defaults to the \code{pval.type}.
#' @param samp.para List of parameters for the tmg truncated normal sampler.
#' Contains two parameters \code{samp.size} and \code{burn.in}, for more details see \code{tmg::rtmg}.
#' @param optim.control allows the user to set some characteristics of the Differential Evolution optimization algorithm implemented in DEoptim.
#' For more details see code{\link[DEoptim]{DEoptim.control}}.
#' @details \code{PSATLinear}  The function is used to perform inference for normal mean vectors
#' that were selected based on a single linear aggregate test. To be exact, suppose
#' that \eqn{y ~ N(\mu,\Sigma)} and that we are interested in estimating \eqn{\mu}
#' only if we can determine that \eqn{\mu != 0} using an aggregate test of the form:
#' \eqn{a'y <l} or \eqn{a'y > u} for some predetermined constants \eqn{a, l, u}.
#'
#' The \code{threshold} parameter specifies the constants \eqn{l<u} which are used
#' to threshold the aggregate test. If only a single number is provided, then the threshold
#' will be set according to test_direction:
#' \itemize{
#' \item lower: a'y < threshold
#' \item upper: a'y > threshold
#' \item two-sided a'y < -threshold, or a'y > threshold}
#' For more details see \href{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12318}{Post-selection estimation and testing following aggregate association tests}.
#' @examples
#' \dontrun{
#' y         <- c(15, rnorm(p - 1))
#' cov.mat <- diag(p)
#' PSATLinear(y, cov.mat,
#' threshold.lower = qnorm(0.025, 0, sqrt(p)),
#' threshold.upper = qnorm(0.975, 0, sqrt(p)))
#' }
#'
#' @export
PSATLinear <- function(y,
                       cov.mat,
                       alpha = 0.05,
                       a.vec = rep(1, ncol(cov.mat)),
                       threshold.lower = -Inf,
                       threshold.upper = Inf,
                       contrast = NULL,
                       alternative = NULL,
                       test.type = c('symmetric', 'UMPU'),
                       pval.type = c('switch', 'global', 'hybrid', 'polyhedral', 'naive'),
                       ci.type   = pval.type,
                       est.type  = c('moment', 'naive'),
                       verbose   = FALSE,
                       zscore.bound = 5,
                       samp.para = list(samp.size = 10000,
                                        burn.in   = 200),
                       optim.control = DEoptim::DEoptim.control()) {
  ### Initalizing paramters
  p         <- length(y)
  pval.df   <- data.frame(matrix(NA, nrow = p, ncol = 0))
  reject.df <- data.frame(matrix(NA, nrow = p, ncol = 0))
  ci.df     <- data.frame(matrix(NA, nrow = p, ncol = 0))
  test.val  <- drop(t(a.vec) %*% y)
  test.var  <- drop(t(a.vec) %*% cov.mat %*% a.vec)
  optim.control$trace <- verbose
  pval.type <- unique(c(ci.type, pval.type))
  # Tests  ------------------------------------------------------------------
  if (is.infinite(threshold.upper) & is.infinite(threshold.lower)) {
    warning('No threshold specified, reverting to rejection area of 1 - alpha')
    threshold.upper <- qnorm(1 - alpha / 2, sd = sqrt(test.var))
    threshold.lower <- qnorm(alpha / 2, sd = sqrt(test.var))
  }
  if (test.val < threshold.upper & test.val > threshold.lower) {
    stop('Test static below threshold, y is not tested.')
  }
  ### Eigen decomposition
  eigen.decomp <- eigen(cov.mat)
  if (any(eigen.decomp$values <= 0)) {
    stop('Covariance matrix is not invertiable')
  }
  if (is.null(contrast)) {
    contrast <- diag(p)
  }
  if (is.null(alternative)) {
    alternative <- rep('two.sided', nrow(contrast))
  }
  if (nrow(contrast) != length(alternative)) {
    stop('Must specifiy which side to test for each contrast')
  }
  thres.list <- thresholdWrapperLinear(test.val = test.val,
                                       a.vec = a.vec,
                                       cov.mat = cov.mat,
                                       first.threshold.lower = threshold.lower,
                                       first.threshold.upper = threshold.upper,
                                       alpha = alpha,
                                       calc.t2 = 'switch' %in% ci.type)
  if (any(c(pval.type, ci.type, est.type) %in% c('hybrid', 'polyhedral', 'switch', 'moment'))) {
    trunc.para.mat <- apply(contrast, 1, findTruncParametersLinear,
                            cov.mat = cov.mat,
                            a.vec = a.vec,
                            threshold.lower = threshold.lower,
                            threshold.upper = threshold.upper,
                            y = y)
  }
  trunc.para.mat <- cbind(t(trunc.para.mat), 'alternative' = sideToNum(alternative))
  if (any(c(pval.type, ci.type) %in% c('hybrid', 'global.null', 'switch'))) {
    samp.null         <- sampleGNLinear(y = y,
                                        null.mu = rep(0, p),
                                        sigma   = cov.mat,
                                        threshold.lower = threshold.lower,
                                        threshold.upper = threshold.upper,
                                        a.vec = a.vec,
                                        num.samp = samp.para$samp.size)
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
                              zscore.bound          = zscore.bound,
                              est.type              = est.type)
  return(list('Pvalues' = pval.df,
              ## 'Rejection' = reject.df,
              'CI' = ci.df,
              'Point.estimation'  = est.df,
              'Threshold.details' = thres.list))
}



### Finding the parameters of the truncated normal
findTruncParametersLinear <- function(y, cov.mat, a.vec,
                                      threshold.lower, threshold.upper, eta) {
  p         <- length(y)
  I.p       <- diag(p)
  etaSigeta <- drop(t(eta) %*% cov.mat %*% eta)
  c.vec     <- (1 / etaSigeta) * (cov.mat %*% eta)
  ## The effect orthogonal to the linear contrast
  W         <- (I.p - c.vec %*% t(eta)) %*% y
  ## Connstant for truncation bounds
  a.c       <- drop(t(a.vec) %*% (c.vec))
  a.w       <- drop(t(a.vec) %*% W)
  trunc.bound <- c('lower' = (threshold.lower - a.w) / a.c,
                   'upper' = (threshold.upper - a.w) / a.c)
  return(c('x' = t(eta) %*% y, 'var' = etaSigeta, 'bounds' = trunc.bound))
}

### Sampling from global null
sampleGNLinear <- function(y, null.mu, sigma,
                            threshold.lower, threshold.upper, a.vec,
                            num.samp = samp.para$samp.size) {
  lower.sample    <- NULL
  upper.sample    <- NULL
  p               <- length(null.mu)
  omega           <- solve(sigma)
  ### Global null test statistic
  gn.test.mu  <- drop(t(a.vec) %*% null.mu)
  gn.test.var <- drop(t(a.vec) %*% sigma %*% a.vec)
  ### Finding the number of observations per each group
  prob.negative    <- pnorm(threshold.lower, gn.test.mu, sqrt(gn.test.var))
  prob.positive    <- pnorm(threshold.upper, gn.test.mu, sqrt(gn.test.var), lower.tail = FALSE)
  percentage.ratio <- prob.negative / prob.positive
  percentage.upper <- 1 / (1 + percentage.ratio)
  samp.size.upper  <- round(percentage.upper * num.samp)
  samp.size.lower  <- num.samp - samp.size.upper
  ### Conditional variance (conditionaing multivariate normal)
  cov.x.sum        <- sigma %*% a.vec
  conditional.var  <- sigma - gn.test.var^(-1) * cov.x.sum %*% t(cov.x.sum)
  ### Crosses lower threshold
  if (samp.size.lower > 0) { ## returns error if no observation are sampled
    mean.dist <- truncnorm::rtruncnorm(n = samp.size.lower,
                                       a = -Inf, b = threshold.lower,
                                       mean = gn.test.mu, sd = sqrt(gn.test.var))
    conditional.mean <- null.mu + gn.test.var^(-1) * cov.x.sum %*% (mean.dist - gn.test.mu)
    lower.sample     <-  MASS::mvrnorm(samp.size.upper, rep(0, p), conditional.var) + t(conditional.mean)
  }
  ### Crosses upper threshold
  ## returns error if no observation are sampled
  if (samp.size.upper > 0) {
    mean.dist <- truncnorm::rtruncnorm(n = samp.size.upper,
                                       a = threshold.upper, b = Inf,
                                       mean = gn.test.mu, sd = sqrt(gn.test.var))
    conditional.mean <- null.mu + gn.test.var^(-1) * cov.x.sum %*% (mean.dist - gn.test.mu)
    upper.sample     <- MASS::mvrnorm(samp.size.upper, rep(0, p), conditional.var) + t(conditional.mean)
  }
  return(rbind(lower.sample, upper.sample))
}



# Old vesrion (do not delete yet) ----------------------------------------

# ### Finding the parameters of the truncated normal
# findTruncParametersLinear <- function(y, cov.mat, a.vec,
#                                       threshold.lower, threshold.upper, eta) {
#   p         <- length(y)
#   I.p       <- diag(p)
#   etaSigeta <- drop(t(eta) %*% cov.mat %*% eta)
#   c.vec     <- (1 / etaSigeta) * (cov.mat %*% eta)
#   ## The effect orthogonal to the linear contrast
#   W         <- (I.p - c.vec %*% t(eta)) %*% y
#   ## Connstant for truncation bounds
#   a.c       <- drop(t(a.vec) %*% (c.vec))
#   a.w       <- drop(t(a.vec) %*% W)
#   trunc.bound <- c('lower.truncation' = (threshold.lower - a.w) / a.c,
#                    'upper.truncation' = (threshold.upper - a.w) / a.c)
#   return(c('x' = t(eta) %*% y, 'var' = etaSigeta, 'bounds' = trunc.bound))
# }
#
# ### Sampling from global null
# sampleGNLinear <- function(y, null.mu, sigma,
#                            threshold.lower, threshold.upper, a.vec,
#                            num.samp = samp.para$samp.size, burn.in = 100) {
#   lower.sample    <- NULL
#   upper.sample    <- NULL
#   p               <- length(null.mu)
#   omega           <- solve(sigma)
#   ### Global null test statistic
#   gn.test.mu  <- t(a.vec) %*% null.mu
#   gn.test.var <- drop(t(a.vec) %*% cov.mat %*% a.vec)
#   ### Finding the number of observations per each group
#   prob.negative    <- pnorm(threshold.lower, gn.test.mu, sqrt(gn.test.var))
#   prob.positive    <- pnorm(threshold.upper, gn.test.mu, sqrt(gn.test.var), lower.tail = FALSE)
#   percentage.ratio <- prob.negative / prob.positive
#   percentage.upper <- 1 / (1 + percentage.ratio)
#   #percentage.lower <- 1 - percentage.upper
#   samp.size.upper  <- round(percentage.upper * num.samp)
#   samp.size.lower  <- num.samp - samp.size.upper
#
#   ### For some reason constraint.list needs to be a list of list, otherwise rtmg returns error
#   ### Crosses lower threshold
#   constraint.list.lower <- list(list(A = matrix(0, nrow = p, ncol = p),
#                                      B = a.vec,
#                                      C = -threshold.lower))
#   if (samp.size.lower > 0) { ## returns error if no observation are sampled
#     lower.sample <- tmg::rtmg(n = (num.samp - samp.size.upper),
#                               M = omega,
#                               r = drop(omega %*% null.mu),
#                               initial = y,
#                               q = constraint.list.lower,
#                               burn.in = burn.in)
#   }
#   ### Crosses upper threshold
#   ## returns error if no observation are sampled
#   constraint.list.upper <- list(a = list(A = matrix(0, nrow = p, ncol = p),
#                                          B = a.vec,
#                                          C = -threshold.upper))
#   if (samp.size.upper > 0) {
#     upper.sample <- tmg::rtmg(n = samp.size.upper,
#                               M = omega,
#                               r = drop(omega %*% null.mu),
#                               initial = y,
#                               q = constraint.list.upper,
#                               burn.in = burn.in)
#   }
#   return(rbind(lower.sample, upper.sample))
# }
