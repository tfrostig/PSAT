# Truncation functions  ---------------------------------------------------

## Calculates the probability on truncated outwards normal distribution
calcProbDoubleTrunc <- function(x, lower.bound, upper.bound, mean, sd) {
  trunc.part <- pnorm(upper.bound, mean, sd) - pnorm(lower.bound, mean, sd)
  return((pnorm(min(lower.bound, x), mean, sd) +
            pnorm(max(upper.bound, x), mean, sd) -
            pnorm(upper.bound, mean, sd)) / (1 - trunc.part))
}

## Finds the density of the function
calcDensDoubleTrun <- function(x, lower.bound, upper.bound, mean, sd) {
  is.bound   <- !(x < upper.bound & x > lower.bound)
  trunc.part <- pnorm(upper.bound, mean, sd) - pnorm(lower.bound, mean, sd)
  dens       <- (is.bound * dnorm(x, mean, sd)) / (1 - trunc.part)
  return(dens)
}


### Find the expected value of T where is truncated normal
calcExpectedTruncNormal <- function(lower.bound, upper.bound, mean, sd) {
  if (lower.bound == -Inf) {
    return(truncnorm::etruncnorm(a = -Inf, b = upper.bound, mean, sd))
  }
  if (lower.bound == Inf) {
    return(truncnorm::etruncnorm(a = -Inf, b = upper.bound, mean, sd))
  }
  upper.prob <- (1 - pnorm(upper.bound, mean, sd)) ### Right tail prob
  lower.prob <- pnorm(lower.bound, mean, sd) ### Left tail prob
  trunc.part <- pnorm(upper.bound, mean, sd) - pnorm(lower.bound, mean, sd) ### Normal prob removed
  if (trunc.part == 1) {
    stop('Error in calculating probabilites, area under truncation is 1')
  }
  upper.part.expect <- truncnorm::etruncnorm(a = -Inf, b = lower.bound, mean, sd) *
    ((1 - trunc.part - upper.prob) / (1 - trunc.part)) ### Rescale expected value (upper part)
  lower.part.expect <- truncnorm::etruncnorm(a = upper.bound, b = Inf, mean, sd) *
    ((1 - trunc.part - lower.prob) / (1 - trunc.part)) ### Rescale expected value (lower part)
  return(upper.part.expect + lower.part.expect)

}

### Truncated normal quantiles
findTruncQuantile <- function(prob, lower.trunc, upper.trunc, mean, sd) {
  lower.prob    <- pnorm(lower.trunc, mean, sd)
  upper.prob    <- pnorm(upper.trunc, mean, sd, lower.tail = FALSE)
  trun.prob     <- lower.prob + upper.prob
  left.side.ind <- prob < lower.prob / trun.prob
  quant <- qnorm(prob * trun.prob, mean, sd) * left.side.ind +
    (1 - left.side.ind) * qnorm(prob * trun.prob - lower.prob + (1 - upper.prob), mean ,sd)
  return(quant)
}

## Example
# lower <- -2
# upper <- 4
# x <- rnorm(10000000, 2, 5)
# x <- x[x < lower | x > upper]
# quantile(x, c(0.025, 0.975))
# lower.trunc <- lower
# upper.trunc <- upper
# mean <- 2
# sd <- 2
#
# findTruncQuantile(c(0.025, 0.975), -2, 4, 2 , 5)


### Function of (X I(x < c1 | X > c2))
calcTruncTest <- function(x, c1, c2, lower.bound, upper.bound, mean, sd) {
  return(x *
           (x < c1 | x > c2) *
           calcDensDoubleTrun(x, lower.bound, upper.bound, mean, sd))
}


calcExpectedTruncTest <-  function(c1, c2, lower.bound, upper.bound, mean, sd) {
  ### The integrator returns errorouns numbers when not split into two
  ### Left part of truncated gaussian
  part.a <- integrate(calcTruncTest,
                      lower = -Inf,
                      upper = min(c1, lower.bound),
                      c1 = c1,
                      c2 = c2,
                      lower.bound = lower.bound,
                      upper.bound = upper.bound,
                      mean,
                      sd)
  ### Middle parts
  #### If c2 < lower.bound
  part.a.b <- integrate(calcTruncTest,
                        lower = min(c2, lower.bound),
                        upper = lower.bound,
                        c1 = c1,
                        c2 = c2,
                        lower.bound = lower.bound,
                        upper.bound = upper.bound,
                        mean,
                        sd)
  #### If c1 > upper.bound
  part.b.b <- integrate(calcTruncTest,
                        lower = upper.bound,
                        upper = max(c1, upper.bound),
                        c1 = c1,
                        c2 = c2,
                        lower.bound = lower.bound,
                        upper.bound = upper.bound,
                        mean,
                        sd)
  ### Right part of truncated gaussian
  part.b <- integrate(calcTruncTest,
                      lower = max(c2, upper.bound),
                      upper = Inf,
                      c1 = c1,
                      c2 = c2,
                      lower.bound = lower.bound,
                      upper.bound = upper.bound,
                      mean,
                      sd)
  return(c('Estimate' = part.a$value + part.b$value + part.a.b$value + part.b.b$value,
           'Error'    = part.a$abs.error + part.b$abs.error + part.a.b$abs.error + part.b.b$abs.error))
}
