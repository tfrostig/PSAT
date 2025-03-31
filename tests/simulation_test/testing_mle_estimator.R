### Generating data
print('Testing PSATquadratic MLE estimation')


all.est.my.psat <- list()
all.est.amit.psat <- list()
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
  }
  psat.quad <- PSATinference::PSATQuadratic(y,
                                            cov.mat,
                                            threshold = threshold,
                                            test.type = 'symmetric',
                                            est.type = c('naive', 'mle', 'moment'))


  psat.amit.quad <- PSAT::mvnQuadratic(y,
                                       cov.mat,
                                       threshold = threshold)

  all.est[[i]] <- cbind(psat.quad$Point.estimation, 'amit' = psat.amit.quad$muhat)

  if (i %% 50 == 0) {
    print(paste('Done with iteration', i, 'out of', iter.num))
  }
}

### Variance
temp <- list()
for (i in 1:iter.num) {
  temp[[i]] <- apply(all.est[[i]], 2, var)
}
var.est <- do.call('rbind', temp)

### MSE
temp <- list()
for (i in 1:iter.num) {
  temp[[i]] <- apply(sweep(all.est[[i]], 1, beta.vec), 2, var)

}
mse.est <- do.call('rbind', temp)

cat('Estimators variance \n')
colMeans(var.est)

cat('Estimators MSE \n')
colMeans(mse.est)

