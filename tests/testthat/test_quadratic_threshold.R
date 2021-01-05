context("Checking the threshold fucntion for linear testing")

## Simulating data
p       <- 10
cov.mat <- diag(p)

set.seed(959)
y         <- c(15, rnorm(p - 1))
a.vec     <- rep(1, p)

first.threshold <- qchisq(0.95, p)
alpha <- 0.05



testthat::test_that('thresholdWrapperQuadratic returns the expected outputs, default parameters', {
  threshold.quadratic <- PSATinference:::thresholdWrapperQuadratic(test.val = t(a.vec) %*% y,
                                                                   first.thres = first.threshold,
                                                                   K = diag(p),
                                                                   alpha = 0.05,
                                                                   cov.mat = cov.mat,
                                                                   calc.t2 = TRUE)
  testthat::expect_equal(threshold.quadratic$t1, alpha)
})
