context("Truncations functions")

library(MASS)
library(PSATinference)
library(testthat)

### Simulating data
n <- 10^7
lower.trunc <- -1.1
upper.trunc <- 1.5
mu <- 0.2
sigma <- 1.5

### Emperical sample
x <- rnorm(n, mu, sigma)
x <- x[x < lower.trunc | x > upper.trunc]


test_that('Truncation probabilites return correct output', {
  temp.t <- rnorm(1)
  sim.est  <- mean(x < temp.t)
  psat.est <- PSATinference:::calcProbDoubleTrunc(temp.t, lower.trunc, upper.trunc, mu, sigma)
  testthat::expect_lte(abs(psat.est - sim.est), 10^-2)
})

test_that('Expectation return correct output', {
  sim.est <- mean(x)
  psat.est <- PSATinference:::calcExpectedTruncNormal(lower.trunc, upper.trunc, mu, sigma)
  testthat::expect_lte(abs(psat.est - sim.est), 10^-2)
})

test_that('Expectation of test (adding rejection areas) return correct output', {
  c1 <- -3
  c2 <- 3.1
  y  <- x

  y[y > c1 & y < c2] <- 0
  sim.est <- mean(y)
  psat.est <- PSATinference:::calcExpectedTruncTest(c1, c2, lower.trunc, upper.trunc, mu, sigma)[1]
  testthat::expect_lte(abs(psat.est - sim.est), 10^-2)
})


test_that('Quantiles function return the correct output', {
  temp.quant <- runif(1)
  sim.est    <- quantile(x, temp.quant)
  psat.est   <- PSATinference:::findTruncQuantile(temp.quant, lower.trunc, upper.trunc, mu, sigma)
  testthat::expect_lte(abs(psat.est - sim.est), 10^-2)
})

