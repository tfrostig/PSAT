### Generating data
print('Testing PSATquadratic testing type I error and coverage rate')

pval.list     <- NULL
coverage.list <- NULL
res.tzviel.pval <- matrix(NA, ncol = p, nrow = iter.num)
for (i in 1:iter.num) {
  ### Creading data
  y          <- MASS::mvrnorm(1, beta.vec, cov.mat)
  test.stat  <- t(y) %*% solve(cov.mat) %*% y
  ### Ensuring that data passes global null test
  count <- 0
  while (test.stat < threshold) {
    y          <- mvrnorm(1, beta.vec, cov.mat)
    test.stat  <- drop(t(y) %*% solve(cov.mat) %*% y)
    count <- count + 1
    print(count)
  }
  psat.quad <- PSATinference::PSATQuadratic(y,
                                            cov.mat,
                                            threshold = threshold,
                                            test.type = 'symmetric')
  res.tzviel.pval[i, ] <- psat.quad$Pvalues$polyhedral.symmetric < 0.05
  pval.list[[i]]       <- psat.quad$Pvalues < alpha
  coverage.list[[i]]   <- (sweep(sapply(psat.quad$CI, function(x) x[ ,1]), 1, beta.vec, '-') < 0) +
    (sweep(sapply(psat.quad$CI, function(x) x[ ,2]), 1, beta.vec, '-') > 0) == 2

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

test.names <- colnames(psat.quad$Pvalues)
ci.names   <- colnames(psat.quad$CI)
### Arranging results pvalue
pval.res <- array(unlist(pval.list),
                  dim = c(nrow(contrast),
                          ncol(psat.quad$Pvalues),
                          length(pval.list)),
                  dimnames = list(NULL, test.names, NULL))
pval.res <- apply(pval.res, c(1, 2), mean)
### Arranging results CI
ci.res <- array(unlist(coverage.list),
                  dim = c(nrow(contrast),
                          ncol(psat.quad$Pvalues),
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




