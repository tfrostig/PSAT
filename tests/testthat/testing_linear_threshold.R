context("Checking the threshold fucntion for linear testing")

## Simulating data
p       <- 10
cov.mat <- diag(p)

set.seed(959)
y         <- c(15, rnorm(p - 1))
a.vec     <- rep(1, p)

first.threshold.lower <- qnorm(0.025, 0, sqrt(sum(a.vec)))
first.threshold.upper <- qnorm(0.975, 0, sqrt(sum(a.vec)))
alpha <- 0.05



testthat::test_that('PSATQuadratic returns the expected outputs, default parameters', {
  threshold.linear <- PSATinference:::thresholdWrapperLinear(test.val = t(a.vec) %*% y,
                                                             first.threshold.upper = first.threshold.upper,
                                                             first.threshold.lower = first.threshold.lower,
                                                             a.vec = a.vec,
                                                             alpha = 0.05,
                                                             cov.mat = cov.mat,
                                                             calc.t2 = TRUE)

  testthat::expect_equal(threshold.linear$t1.lower, alpha / 2)
  testthat::expect_equal(threshold.linear$t1.upper, alpha / 2)
})
