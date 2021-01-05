### Checking confidence intervals
p             <- 10
cov.mat       <- matrix(0, ncol = p, nrow = p)
diag(cov.mat) <- 1
beta.vec      <- c(0.5, rep(0, p - 2), 0.4)
threshold     <- qchisq(1 - 10^-2, p)

iter.num      <- 1000
### WALD test for ease of computation
library(MASS)
res.tzviel.ci   <- matrix(NA, ncol = p, nrow = iter.num)
res.tzviel.pval <- matrix(NA, ncol = p, nrow = iter.num)

res.amit   <- matrix(NA, ncol = p, nrow = iter.num)
res.tzviel.ci.global <- matrix(NA, ncol = p, nrow = iter.num)

for (i in 1:iter.num) {
  ### Creading data
  y          <- mvrnorm(1, beta.vec, cov.mat)
  test.stat  <- t(y) %*% solve(cov.mat) %*% y

  while (test.stat < threshold) {
    y          <- mvrnorm(1, beta.vec, cov.mat)
    test.stat  <- drop(t(y) %*% solve(cov.mat) %*% y)
  }
  ### PSAT
  # amit.psat <- PSAT::mvnQuadratic(y,
  #                                 cov.mat,
  #                                 ci_type = 'polyhedral',
  #                                 pvalue_type = 'polyhedral',
  #                                 threshold = threshold)
  tzvi.psat <- PSATinference::PSATQuadratic(y,
                                            cov.mat,
                                            pval.type = c('polyhedral', 'global.null'),
                                            test.type = 'symmetric',
                                            threshold = threshold)
  ### Coverage rate
  res.tzviel.ci[i, ] <- (tzvi.psat$CI$polyhedral.symmetric[ ,1] < beta.vec) &
    (tzvi.psat$CI$polyhedral.symmetric[ ,2] > beta.vec)
  res.tzviel.ci.global[i, ] <- (tzvi.psat$CI$global[ ,1] < beta.vec) &
    (tzvi.psat$CI$global[ ,2] > beta.vec)
  res.tzviel.pval[i, ] <- tzvi.psat$Pvalues$polyhedral.symmetric < 0.05

  # res.amit[i, ] <- (amit.psat$polyCI[ ,1] < beta.vec) &
  #   (amit.psat$polyCI[ ,2] > beta.vec)
  #res.amit[i, ] <- amit.psat$polyPval
  if (i %% 50 == 0) {
    print(paste('Done with iteration', i, 'out of', iter.num))
  }
}
colMeans(res.tzviel.ci.global)
colMeans(res.amit)
colMeans(res.tzviel.pval)

### Checking linear
set.seed(999)
library(MASS)

iter.num      <- 100
p             <- 10
cov.mat       <- matrix(0, ncol = p, nrow = p)
diag(cov.mat) <- 1
beta.vec      <- c(0.5, rep(0, p - 2), 0.4)
t
### Just a sum test
a         <- rep(1, p)
threshold <- qnorm(1 - 0.05, 0, sqrt(t(a) %*% cov.mat %*% a))


res.tzviel.ci.poly   <- matrix(NA, ncol = p, nrow = iter.num)
res.tzviel.ci.global <- matrix(NA, ncol = p, nrow = iter.num)
res.tzviel.pval      <- matrix(NA, ncol = p, nrow = iter.num)

res.amit   <- matrix(NA, ncol = p, nrow = iter.num)

for (i in 1:iter.num) {
  ### Creading data
  y          <- mvrnorm(1, beta.vec, cov.mat)
  test.stat  <- t(a) %*% y

  while (test.stat < threshold) {
    y          <- mvrnorm(1, beta.vec, cov.mat)
    test.stat  <- drop( t(a) %*% y)
  }
  ### PSAT
  # amit.psat <- PSAT::mvnLinear(y = y,
  #                              sigma = cov.mat,
  #                              testVec = a,
  #                              ci_type = 'polyhedral',
  #                              pvalue_type = 'polyhedral',
  #                              threshold = threshold,
  #                              test_direction = 'upper')
  tzvi.psat <- PSATinference::PSATLinear(y = y,
                                         cov.mat = cov.mat,
                                         a.vec = a,
                                         pval.type = c('polyhedral', 'naive', 'global.null'),
                                         test.type = 'symmetric',
                                         threshold.upper = threshold)
  ### Coverage rate
  # res.tzviel.ci.poly[i, ] <- (tzvi.psat$CI$polyhedral.symmetric[ ,1] < beta.vec) &
  #   (tzvi.psat$CI$polyhedral.symmetric[ ,2] > beta.vec)
  res.tzviel.ci.global[i, ] <- (tzvi.psat$CI$global[ ,1] < beta.vec) &
    (tzvi.psat$CI$global[ ,2] > beta.vec)
  res.tzviel.pval[i, ] <- tzvi.psat$Pvalues$global
  # res.amit[i, ] <- (amit.psat$polyCI[ ,1] < beta.vec) &
  #   (amit.psat$polyCI[ ,2] > beta.vec)
  #res.amit[i, ] <- amit.psat$polyPval
  if (i %% 50 == 0) {
    print(paste('Done with iteration', i, 'out of', iter.num))
  }
}

colMeans(res.tzviel.ci.global, na.rm = TRUE)
#colMeans(res.amit)
colMeans(res.tzviel.pval < 0.05, na.rm = TRUE)
colMeans(res.tzviel.ci.global < 0.05, na.rm = TRUE)
