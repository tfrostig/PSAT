# Point estimation  -------------------------------------------------------

estimationWrapper <- function(truncation.parameters,
                              y,
                              est.type,
                              test.mat,
                              zscore.bound,
                              cov.mat,
                              threshold,
                              optim.method) {
  est.df <- data.frame(matrix(NA, nrow = nrow(truncation.parameters), ncol = 0))
  if ('naive' %in% est.type) {
    est.df['naive'] <- truncation.parameters[ , 'x']
  }
  if ('moment' %in% est.type) {
    est.df['moment'] <- momentEstimation(truncation.parameters = truncation.parameters,
                                         zscore.bound = zscore.bound)
  }
  if( "mle" %in% est.type) {
    if(optim.method == "sgd") {
      quadlam <- getQudraticLam(test.mat, cov.mat = cov.mat)
      if(any(quadlam < 10e-8)) {
        optim.method <- "Nelder-Mead"
        warning("Test matrix nearly singualr: using Nelder-Mead for MLE computation")
      } else {
        sgdfit <- quadraticSGD(y, cov.mat, solve(cov.mat), test.mat, threshold,
                               stepRate = 0.75, stepCoef = NULL,
                               delay = 20, assumeCovergence = 800,
                               mhIters = 40)
        mleMu <- sgdfit$mle
        solutionPath <- sgdfit$solutionPath
      }
    }
    if(optim.method == "Nelder-Mead") {
      mleMu <- quadraticNM(y, cov.mat, solve(cov.mat), test.mat, threshold)
      solutionPath <- NULL
    }
    est.df['mle'] <- mleMu
  }
  return(est.df)
}

momentEstimation <- function(truncation.parameters, zscore.bound, search.jump = 0.25) {
  moment.est <- rep(NA, nrow(truncation.parameters))
  for (i in 1:nrow(truncation.parameters)) {
    temp.x   <- truncation.parameters[i, 'x']
    temp.sig <- sqrt(truncation.parameters[i, 'var'])
    temp.lower.bound <- truncation.parameters[i, 'bounds.lower']
    temp.upper.bound <- truncation.parameters[i, 'bounds.upper']
    if (is.infinite(temp.lower.bound) & is.infinite(temp.upper.bound)) {
      moment.est[i] <- temp.x ### Naive estimator (no adjustment)
    } else {
      ### Ensuring that the search bounds for uniroot are of opposite signs
      search.bounds <- c(0, 0)
      search.vec    <- c(temp.x, temp.x)
      is.search.inf <- is.infinite(search.bounds)
      while (sign(search.bounds[1]) == sign(search.bounds[2]) &
             all(!is.search.inf)) {
        search.vec <- c(search.vec[1] - search.jump, search.vec[2] + search.jump)
        search.bounds <- muFunction(x   = temp.x,
                                    mu  = search.vec,
                                    sig = temp.sig,
                                    lower.bound = temp.lower.bound,
                                    upper.bound = temp.upper.bound)
        is.search.inf <- is.infinite(search.bounds)
      }
      if (any(is.search.inf)) {
        ### muFunction can reach Inf, usually when pnorm/dnorm very small
        ### if so return the closest value to 0 reached.
        moment.est[i] <- search.vec[is.search.inf]
      }
      if (all(!is.search.inf)) {
        moment.est[i] <- uniroot(muFunction,
                                 x = temp.x,
                                 sig = temp.sig,
                                 lower.bound = temp.lower.bound,
                                 upper.bound = temp.upper.bound,
                                 interval =  c(search.vec[1],
                                               search.vec[2]))$root
      }
    }
  }
  return(moment.est)
}


muFunction <- function(x, mu, sig, lower.bound, upper.bound) {
  return(x -
           (mu -
              sig * (dnorm((lower.bound - mu) / sig) - dnorm((upper.bound - mu) / sig)) /
              (pnorm((lower.bound - mu) / sig) + 1 - pnorm((upper.bound - mu) / sig))))
}


