context("Returning expected outputs")

library(MASS)
library(PSATinference)
library(testthat)
### Simulating data
p       <- 100
cov.mat <- diag(p)

set.seed(959)
y         <- c(15, rnorm(p - 1))
cov.mat <- diag(p)
list.psat.quadratic <- PSATQuadratic(y, cov.mat)
list.psat.linear    <- PSATLinear(y, cov.mat,
                                  threshold.lower = qnorm(0.025, 0, sqrt(p)),
                                  threshold.upper = qnorm(0.975, 0, sqrt(p)))


test_that('PSATQuadratic returns the expected outputs, default parameters', {


  expect_equal(colnames(list.psat.quadratic$Pvalues), c('polyhedral.symmetric',
                                              'polyhedral.umpu',
                                              'naive',
                                              'global',
                                              'hybrid.symmetric',
                                              'hybrid.umpu'))

  expect_equal(colnames(list.psat.quadratic$CI), c('polyhedral.symmetric',
                                         'polyhedral.umpu',
                                         'naive',
                                         'global',
                                         'hybrid.symmetric',
                                         'hybrid.umpu'))
})


test_that('PSATLinear returns the expected outputs, default parameters', {


  expect_equal(colnames(list.psat.linear$Pvalues), c('polyhedral.symmetric',
                                                     'polyhedral.umpu',
                                                     'naive',
                                                     'global',
                                                     'hybrid.symmetric',
                                                     'hybrid.umpu'))

  expect_equal(colnames(list.psat.linear$CI), c('polyhedral.symmetric',
                                                'polyhedral.umpu',
                                                'naive',
                                                'global',
                                                'hybrid.symmetric',
                                                'hybrid.umpu'))
})


test_that('Linear returns warning when not specifying threshold', {
  expect_warning(PSATLinear(y, cov.mat))
})


