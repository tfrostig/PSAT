### Generating data
print('Testing PSATLinear testing type I error and coverage rate')

a         <- rep(1, p)
threshold <- qnorm(1 - 0.05, 0, sqrt(t(a) %*% cov.mat %*% a))


pval.list     <- NULL
coverage.list <- NULL
res.tzviel.pval <- matrix(NA, ncol = p, nrow = iter.num)
for (i in 1:iter.num) {
  ### Creading data
  y          <- MASS::mvrnorm(1, beta.vec, cov.mat)
  test.stat  <- drop(t(a) %*% y)
  ### Ensuring that data passes global null test
  while (test.stat < threshold) {
    y          <- mvrnorm(1, beta.vec, cov.mat)
    test.stat  <- drop(t(a) %*% y)
  }
  psat.linear <- PSATinference::PSATLinear(y,
                                           cov.mat = cov.mat,
                                           a.vec = a,
                                           threshold.upper = threshold,
                                           test.type = 'symmetric')
  res.tzviel.pval[i, ] <- psat.linear$Pvalues$polyhedral.symmetric < 0.05
  pval.list[[i]]       <- psat.linear$Pvalues < alpha
  coverage.list[[i]]   <- (sweep(sapply(psat.linear$CI, function(x) x[ ,1]), 1, beta.vec, '-') < 0) +
    (sweep(sapply(psat.linear$CI, function(x) x[ ,2]), 1, beta.vec, '-') > 0) == 2

  if (i %% 50 == 0) {
    print(paste('Done with iteration', i, 'out of', iter.num))
  }
}
temp <- NULL

for (i in 1:iter.num) {
  temp <- c(temp, pval.list[[i]][1,1])
}
mean(temp)
colMeans(res.tzviel.pval)

test.names <- colnames(psat.linear$Pvalues)
ci.names   <- colnames(psat.linear$CI)
### Arranging results pvalue
pval.res <- array(unlist(pval.list),
                  dim = c(nrow(contrast),
                          ncol(psat.linear$Pvalues),
                          length(pval.list)),
                  dimnames = list(NULL, test.names, NULL))
pval.res <- apply(pval.res, c(1, 2), mean)
### Arranging results CI
ci.res <- array(unlist(coverage.list),
                dim = c(nrow(contrast),
                        ncol(psat.linear$Pvalues),
                        length(pval.list)),
                dimnames = list(NULL, ci.names, NULL))
ci.res <- apply(ci.res, c(1, 2), mean)



cat('Results of type I error of tests\n',
    'Non-null contrasts\n')
print(pval.res[null.ind, ])
cat('Null contrasts\n')
print(pval.res[non.null.ind, ])


cat('Corverage rates of \n',
    'Non-null contrasts\n')
print(ci.res[null.ind, ])
cat('Null contrasts\n')
print(ci.res[non.null.ind, ])
