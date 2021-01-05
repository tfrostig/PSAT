context("Returning expected outputs for estimation")

library(MASS)
library(PSATinference)
library(testthat)
### Simulating data
p       <- 100
cov.mat <- diag(p)

set.seed(959)
y         <- c(15, rnorm(p - 1))
list.psat.quadratic <- PSATQuadratic(y, cov.mat,
                                     ci.type = 'polyhedral', pval.type = 'polyhedral', test.type = 'symmetric')
list.psat.linear    <- PSATLinear(y, cov.mat,
                                  threshold.lower = qnorm(0.025, 0, sqrt(p)),
                                  threshold.upper = qnorm(0.975, 0, sqrt(p)),
                                  ci.type = 'polyhedral', pval.type = 'polyhedral', test.type = 'symmetric')


test_that('PSATQuadratic returns the expected outputs for estimation, default parameters', {


  expect_equal(colnames(list.psat.quadratic$Point.estimation), c('naive', 'moment'))

})


test_that('PSATLinear returns the expected outputs for estimation, default parameters', {

  expect_equal(colnames(list.psat.linear$Point.estimation), c('naive', 'moment'))

})



