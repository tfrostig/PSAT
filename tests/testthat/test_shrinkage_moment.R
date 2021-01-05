context("Checking moment estimator shrinkage")

library(MASS)
library(PSATinference)
library(testthat)
### Simulating data
p       <- 100
cov.mat <- diag(p)

set.seed(959)
y         <- c(-4, rnorm(p - 1))
while (sum(y) > qnorm(0.025, 0, sqrt(p)) | sum(y^2) < qchisq(0.95, p)) {
  y         <- c(-4, rnorm(p - 1))
}
cov.mat <- diag(p)
list.psat.quadratic <- PSATQuadratic(y, cov.mat, test.type = 'symmetric')
list.psat.linear    <- PSATLinear(y, cov.mat,
                                  threshold.lower = qnorm(0.025, 0, sqrt(p)),
                                  threshold.upper = qnorm(0.975, 0, sqrt(p)),
                                  test.type = 'symmetric')


test_that('PSATQuadratic returns the expected outputs, default parameters', {
  testthat::expect_true(any(list.psat.linear$Point.estimation$naive < list.psat.linear$Point.estimation$moment))
  testthat::expect_true(any(list.psat.quadratic$Point.estimation$naive < list.psat.quadratic$Point.estimation$moment))

})
