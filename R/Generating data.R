## Generating data for GWAS 
library(MASS)
set.seed(999)

p       <- 150 
rho     <- 0.7
cov.mat <- rho^abs(outer(1:p, 1:p, "-"))


#### Create SNPs transforms normal vector into 0,1,2 while keeping MAF
snpMaker <- function(x.vec, maf) {
  g.vec <- x.vec 
  g.vec[x.vec <= qnorm(1 - maf)] <- 0
  g.vec[x.vec >  qnorm(1 - maf)]  <- 1
  g.vec[x.vec >  qnorm(1 - (1 / 3) * maf)]  <- 2
  return(g.vec)
}

x           <- mvrnorm(2000, rep(0, p), cov.mat)
maf.vec     <- rbeta(p, 1.1 ,3) 
G           <- t(apply(x, 1, snpMaker, maf.vec))
colnames(G) <- paste0(c(rep('Gene_a', 10), rep('Gene_b', 50), rep('Gene_c', 10)), '_', c(1:10, 1:(p - 20), 1:10))

E           <- mvrnorm(2000, rep(10, 2), diag(c(1, 0.5))) 
E[ ,2]      <- round(E[ ,2])
colnames(E) <- c('Avg.Calories', 'Num.Sibilings')


beta.vec <- c(0.15, rep(0, p - 4), rep(0.2, 3)) 
b.vec    <- c(1, 0)
y.pheno  <- 15 + drop(G %*% beta.vec) + drop(E %*% b.vec) + rnorm(2000)


summary(lm(y.pheno ~ G + E))

