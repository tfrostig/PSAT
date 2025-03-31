## code taken from Amit Meir. https://github.com/ammeir2/PSAT/blob/master/R/quadraticSGD.R

quadraticSGD <- function(y, sigma, precision, testmat, threshold,
                         stepRate = 0.75, stepCoef = NULL,
                         delay = 10,
                         sgdSteps = 800, assumeCovergence = 400,
                         mhIters = 40) {

  if(is.null(stepCoef)) {
    stepCoef <- 0.25 / sqrt(diag(sigma))
  } else if(length(stepCoef) == 1) {
    stepCoef <- stepCoef / sqrt(diag(sigma))
  }

  # Pre-processing test matrix ---------
  testEigen <- eigen(testmat)
  sqrtTestMat <- testEigen$vectors %*% diag(sqrt(testEigen$values))
  sqrtTestMat <- sqrtTestMat %*% t(testEigen$vectors)
  invSqrtTestMat <- testEigen$vectors %*% diag(1 / sqrt(testEigen$values))
  invSqrtTestMat <- invSqrtTestMat %*% t(testEigen$vectors)
  # sampSig <- t(sqrtTestMat) %*% sigma %*% sqrtTestMat
  # sampP <- solve(sampSig)

  # Setting up SGD ----------------
  itermu <- y
  mu <- y
  mu <- sqrtTestMat %*% mu
  lastsamp <- as.numeric(sqrtTestMat %*% y)
  solutionPath <- matrix(ncol = length(mu), nrow = sgdSteps)
  sampmat <- matrix(ncol = length(y), nrow = sgdSteps)
  test.list <- list(a = list('A' = testmat, 'B' = rep(0, length(y)), 'C' = -threshold))
  for(i in 1:sgdSteps) {
    # Sampling -----------------
    sampMu <- drop(sqrtTestMat %*% itermu)
    ### replaced with tmg::rtmg (to not carry dependencies)
    sample <- tmg::rtmg (initial = lastsamp,
                         r = sampMu,
                         M = precision,
                         q = test.list,
                         n = 2,
                         burn.in = round(mhIters) / 2)
    lastsamp <- sample[2, ]

    # computing gradient ------------
    grad <- sample %*% invSqrtTestMat
    grad <- colMeans(grad)
    grad <- precision %*% (y - grad) * stepCoef / max(1, i - delay)^stepRate
    grad <- sign(grad) * pmin(abs(grad), sqrt(diag(sigma)) * 0.05)

    # Updating estimate --------------
    itermu <- itermu + grad
    itermu <- pmax(0, itermu * sign(y)) * sign(y)
    itermu <- pmin(abs(y), abs(itermu)) * sign(y)
    solutionPath[i, ] <- itermu
    sampmat[i, ] <- lastsamp
  }

  # Reporting ----------------------
  mle <- colMeans(solutionPath[assumeCovergence:sgdSteps, ])
  return(list(mle = mle, solutionPath = solutionPath, sampMat = sampmat))
}


evaluateQuadraticLikelihood <- function(mu, y, threshold,
                                        precision, lam, deltamat) {
  diff <- y - mu
  dens <- -0.5 * (t(diff) %*%  precision %*% diff)
  delta <- as.numeric(deltamat %*% mu)
  prob <-  CompQuadForm::liu(threshold, lambda = lam, delta = delta^2)
  return(-dens + log(prob))
}

quadraticNM <- function(y, sigma, precision, testMat, threshold) {
  liuParams <- getQudraticLam(testMat, sigma, FALSE)
  lam <- liuParams$lam
  deltamat <- liuParams$deltamat
  result <- optim(par = y, fn = evaluateQuadraticLikelihood,
                  y = y, threshold = threshold, precision = precision,
                  lam = lam, deltamat = deltamat,
                  method = "Nelder-Mead")$par
  return(result)
}
