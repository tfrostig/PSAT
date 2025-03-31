
# Optimization functions  -------------------------------------------------
### Find Uniformly most powerful test

#' #' Inference for Normal Means after Aggregate Testing
#'
#' @description \code{findUMPU} is to find either the UMPU rejection area, or the p-value
#' corresponding to that area.
#'
#' @param test.stat the test statistic. Not required for finding the rejection area.
#' @param alpha type I error to be kept by the rejection area. Not required for finding p-value.
#' @param lower.bound lower bound of the truncated normal distribution.
#' @param upper.bound upper bound of the truncated normal distribution.
#' @param mean the expected value of the truncated normal distribution.
#' @param sd the standard deviation of the truncated normal distribution
#' @param optim.type accepts, either rejection or pvalue, indicates which value to find.
#' @param sd.dist defines the bounds of the search, will not search for rejection area more than sd.dist Z-scores from the mean.
#' @param optim.control allows the user to set some characteristics of the Differential Evolution optimization algorithm implemented in DEoptim.
#' For more details see code{\link[DEoptim]{DEoptim.control}}.
#'
#' @details The function finds either the rejection area or p-value under the conditions of an UMPU test,
#' as described in section 3.1 \href{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12318}{Post-selection estimation and testing following aggregate association tests}.

findUMPU <- function(test.stat = NA,
                      alpha = NA,
                      lower.bound, upper.bound,
                      mean, sd,
                      optim.type = 'rejection',
                      sd.dist = 10,
                      optim.control) {
  target.b <- alpha * calcExpectedTruncNormal(lower.bound,
                                              upper.bound,
                                              mean,
                                              sd)
  if (optim.type == 'rejection') {
    optim.res <- DEoptim::DEoptim(targetFunctionOptimRejection,
                         lower  = -c(mean + sd.dist * sd, mean + sd.dist * sd),
                         upper  = c(mean + sd.dist * sd, mean + sd.dist * sd),
                         #lower = -c(Inf, Inf),
                         #upper = c(Inf, Inf),
                         lower.bound = lower.bound,
                         upper.bound = upper.bound,
                         mean = mean, sd = sd,
                         alpha = alpha,
                         target.b = target.b,
                         control = optim.control)
  }
  if (optim.type == 'pvalue') {
    optim.res <- DEoptim::DEoptim(targetFunctionOptimPval,
                         lower  = c(0 + 10^-6, -(mean + sd.dist * sd)),
                         upper  = c(1 - 10^-6, mean + sd.dist * sd),
                         test.stat   = test.stat,
                         lower.bound = lower.bound,
                         upper.bound = upper.bound,
                         mean = mean, sd = sd,
                         control = optim.control)
  }
  return(optim.res)
}


### Target to minimize (sum of square of two functions)
targetFunctionOptimRejection <- function(c.vec, lower.bound, upper.bound, mean, sd, alpha, target.b) {
  res.eq1 <- calcProbDoubleTrunc(c.vec[1], lower.bound, upper.bound, mean, sd) +
            (1 - calcProbDoubleTrunc(c.vec[2], lower.bound, upper.bound, mean, sd)) -
            alpha
  res.eq2 <- calcExpectedTruncTest(c.vec[1], c.vec[2],
                                   lower.bound, upper.bound,
                                   mean, sd) - target.b
  return(res.eq1^2 + res.eq2[1]^2)
}



### Functions aimed to find the UMVUE test p-value
### Finds p-value according to UMVUE equations (should input c1)
targetFunctionOptimPval <- function(alpha.c, test.stat, lower.bound, upper.bound, mean, sd) {
  ### Initalizing parameters
  c.vec <- c(min(alpha.c[2], test.stat), max(alpha.c[2], test.stat))
  alpha <- alpha.c[1]
  ### Equations used in UMPU
  target.b <- alpha * calcExpectedTruncNormal(lower.bound,
                                              upper.bound,
                                              mean,
                                              sd)
  res.eq1 <- calcProbDoubleTrunc(c.vec[1], lower.bound, upper.bound, mean, sd) +
    (1 - calcProbDoubleTrunc(c.vec[2], lower.bound, upper.bound, mean, sd)) -
    alpha
  res.eq2 <- calcExpectedTruncTest(c.vec[1], c.vec[2],
                                   lower.bound, upper.bound,
                                   mean, sd) - target.b
  return(res.eq1^2 + res.eq2[1]^2)
}


### Finding CIs
### Inverting tests according to the tests
### Finds the minimal or maximal mu that will lead to reject a
### test statistic at probability of type one error prob
findCIoneSide <- function(stat, lower.trunc, upper.trunc, sd, prob, sd.dist) {
  ### Checking the limit of the search defined for the function
  lower.val <- polyCIOneSideToOptim(mu = stat - sd.dist * sd,
                                    stat = stat,
                                    lower.trunc = lower.trunc,
                                    upper.trunc = upper.trunc,
                                    sd = sd,
                                    prob = prob)
  upper.val <- polyCIOneSideToOptim(mu = stat + sd.dist * sd,
                                    stat = stat,
                                    lower.trunc = lower.trunc,
                                    upper.trunc = upper.trunc,
                                    sd = sd,
                                    prob = prob)
  if (sign(lower.val) == sign(upper.val)) {
     return(ifelse(lower.val < upper.val, stat - sd.dist * sd, stat + sd.dist * sd))
  }
  ### Making sure that uniroot can solve the equation
  if (sign(lower.val) != sign(upper.val)) {
    bound <- uniroot(polyCIOneSideToOptim,
                     lower = stat - sd.dist * sd,
                     upper = stat + sd.dist * sd,
                     stat = stat,
                     lower.trunc = lower.trunc,
                     upper.trunc = upper.trunc,
                     sd = sd,
                     prob = prob)
    return(bound$root)
  }
}

#' Returns the difference between the quantile of a truncated normal with a certain probability
#' to the test statistic.
polyCIOneSideToOptim <- function(mu, stat, lower.trunc, upper.trunc, sd, prob) {
  temp.rej <- findTruncQuantile(prob = prob,
                                lower.trunc = lower.trunc,
                                upper.trunc = upper.trunc,
                                mean = mu,
                                sd = sd)
  ans <- (temp.rej - stat)
  return(ans)
}



polyCIUMPUToOptim <- function(prob, stat, lower.trunc, upper.trunc, sd, alpha, sd.dist = 10) {
  ci <- c('lower' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, 1 - alpha +  prob, sd.dist),
          'upper' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, prob, sd.dist))
  ci.length <- ci[2] - ci[1]
  return(ci.length)
}


#' #' Inference for Normal Means after Aggregate Testing
#'
#' @description \code{findUMPU} is to find either the UMPU rejection area, or the p-value
#' corresponding to that area.
#'
#' @param stat the test statistic.
#' @param lower.trunc lower bound of the truncated normal distribution.
#' @param upper.trunc upper bound of the truncated normal distribution.
#' @param sd the standard deviation of the truncated normal distribution.
#' @param alpha \eqn{1 - \alpha} confidence level for the CI.
#' @param sd.dist defines the bounds of the search, will not search for rejection area more than sd.dist Z-scores from the mean.
#' For more details see code{\link[DEoptim]{DEoptim.control}}.
#'
#' @details The functions finds the UMAU confidecne interval.
#' Based on Pratt's theorem \href{https://www.tandfonline.com/doi/abs/10.1080/01621459.1961.10480644}{Pratt Theorem}, instead of directly inverting the UMPU test,
#' finding the shortest CI will yield the same CI.

polyCIUMPU <- function(stat, lower.trunc, upper.trunc, sd, alpha, sd.dist) {
  best.prob <- optimize(polyCIUMPUToOptim,
                        stat = stat,
                        lower.trunc = lower.trunc,
                        upper.trunc = upper.trunc,
                        sd = sd,
                        alpha = alpha,
                        lower = 0, upper = alpha,
                        sd.dist = sd.dist)
  prob   <- best.prob$minimum
  return(c('lower' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, 1 - alpha +  prob, sd.dist),
            'upper' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, prob, sd.dist)))
}
