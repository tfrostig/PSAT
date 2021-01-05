context("Returning expected outputs for non default parameters")

library(MASS)
library(PSATinference)
library(testthat)
### Simulating data
p       <- 100
cov.mat <- diag(p)

set.seed(959)
y         <- c(10, rnorm(p - 1))
contrast  <- rbind(rep(1, p), c(1, 1, 1, rep(0, p - 3)))
list.psat.quadratic <- PSATQuadratic(y = y, cov.mat = cov.mat,
                                     contrast = contrast, alternative =  c('less', 'greater'))


test_that('PSATQuadratic returns the error when test side are not defined correctly', {
  ### test.side is ill defined
  expect_error(PSATQuadratic(y = y, cov.mat = cov.mat,
                             contrast = contrast, alternative = c('sd', 'upper')))

  ### Not matches the number of contarsts
  expect_error(PSATQuadratic(y = y, cov.mat = cov.mat,
                             contrast = contrast, alternative = c('greater')))

})

