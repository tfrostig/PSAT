# Rejection areas ---------------------------------------------------------
rejectionAreaWrapper <- function(truncation.parameters,
                                 reject.naive,
                                 sample.global.null,
                                 contrast.y, var.contrast.y,
                                 contrast, alpha,
                                 test.type,
                                 rejection.area.type,
                                 optim.control) {
  reject.df <- data.frame(matrix(NA, nrow = ncol(truncation.parameters), ncol = 0))
  if (any(rejection.area.type %in% c('hybrid', 'polyhedral'))) {
    if (any(test.type == 'symmetric')) {
      reject.symm <-  t(apply(truncation.parameters, 2, polyReject,
                              test.type = 'symmetric'))
      colnames(reject.symm) <- c('lower', 'upper')
      reject.df['polyhedral.symmetric'] <- reject.symm
    }
    if (any(test.type == 'UMPU')) {
      reject.umpu <- t(apply(truncation.parameters, 2, polyReject,
                              test.type = 'UMPU',
                              optim.control = optim.control))
      colnames(reject.umpu) <- c('lower', 'upper')
      reject.df['polyhedral.umpu'] <- reject.umpu
    }
  }
  ## Naive
  if (any(rejection.area.type == 'naive')) {
    reject.naive       <- var.contrast.y %*% t(qnorm(c('lower' = alpha / 2,
                                                       'upper' = 1 - alpha / 2)))
    colnames(reject.naive) <- c('lower', 'upper')
    reject.df['naive'] <- reject.naive
  }
  ## Global null test
  if (any(rejection.area.type %in% c('hybrid', 'global.null'))) {
    reject.global <- t(globalReject(sample.mat = sample.global.null, alpha = alpha))
    colnames(reject.global) <- c('lower', 'upper')
    reject.df['global'] <- reject.global
  }
  return(reject.df)
}




### Finding rejection area
polyReject <- function(trunc.parameters, test.type = 'symmetric', optim.control, alpha = 0.05) {
  x     <- trunc.parameters[1]
  sd    <- sqrt(trunc.parameters[2])
  lower <- trunc.parameters[3]
  upper <- trunc.parameters[4]
  ### Condition
  minimal.zscore <- (sign(lower) == sign(upper)) * min(abs(c(lower, upper))) / sd
  ### No selection effect
  if (lower == -Inf | upper == Inf | minimal.zscore > 5) {
    ### Rejection area
    rej.area.symm <- c('lower' = qnorm(alpha / 2, mean = 0, sd = sd),
                       'upper' = qnorm(1 - alpha / 2, mean = 0, sd = sd))
    return(rej.area.symm)
  }
  ### Selection effect, symmetric p-values
  if (test.type == 'symmetric') {
    rej.area.trunc.symm  <- c(findTruncQuantile(c(alpha / 2, 1 - alpha / 2),
                                                lower.trunc = lower.upper,
                                                upper.trunc = upper.trunc,
                                                mean = 0, sd = sd))
    return(rej.area.trunc.symm)
  }
  ### Selection effect, UMPU p-values
  if (test.type == 'UMPU') {
    rej.area.umpu <- findUMPU(test.stat = NA,
                                alpha = alpha,
                                lower.bound = lower,
                                upper.bound = upper,
                                mean = 0,
                                sd = sd,
                                optim.type = 'rejection',
                                optim.control = optim.control)
    rej.area.umpu <- rej.area.umpu$optim$bestmem
    return(rej.area.umpu)
  }
}

## Rejection area for globalNull test
globalReject <- function(sample.mat, alpha) {
  apply(sample.mat, 2, function(x) quantile(x, probs = c('lower' = alpha / 2,
                                                         'upper' = 1 - alpha / 2)))
}

