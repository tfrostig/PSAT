library(PSATinference)
library(MASS)

### Simulation parameter
iter.num      <- 1000
alpha         <- 0.05
### Data generation
p             <- 10
cov.mat       <- 0.7^abs(outer(1:10, 1:10, '-'))
diag(cov.mat) <- 1
cov.mat
beta.vec      <- c(0.9, rep(0, p - 2), 0.4)
threshold     <- qchisq(1 - 10^-2, p)
contrast      <- diag(p)

non.null.ind  <- which(!(beta.vec %*% contrast == 0))
null.ind      <- setdiff(1:nrow(contrast), non.null.ind)



