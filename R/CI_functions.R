# CI Functions  -----------------------------------------------------------
### Applying CI's
CIwrapper <- function(truncation.parameters,
                      reject.naive,
                      sample.global.null,
                      contrast, alpha,
                      test.type,
                      ci.type,
                      thres.list,
                      zscore.bound) {
  ci.df     <- data.frame(matrix(NA, nrow = nrow(truncation.parameters), ncol = 0))
  if (any(ci.type %in% c('hybrid', 'polyhedral'))) {
    if (any(test.type == 'symmetric')) {
      ci.df['polyhedral.symmetric'] <- t(apply(truncation.parameters, 1, polyCI,
                                               test.type = 'symmetric',
                                               alpha = alpha,
                                               zscore.bound = zscore.bound))
      colnames(ci.df$polyhedral.symmetric) <- c('lower', 'upper')
    }
    if (any(test.type == 'UMPU')) {
      ci.df['polyhedral.umpu'] <- t(apply(truncation.parameters, 1, polyCI,
                                          test.type = 'UMPU',
                                          alpha = alpha,
                                          zscore.bound = zscore.bound))
      colnames(ci.df$polyhedral.umpu) <- c('lower', 'upper')
    }
  }
  ## Naive
  if (any(ci.type == 'naive')) {
    conf.interval.naive <- t(apply(truncation.parameters, 1, polyCI,
                                   test.type = 'naive',
                                   alpha = alpha,
                                   zscore.bound))
    ci.df['naive']      <- conf.interval.naive
  }
  if (any(ci.type %in% c('hybrid', 'global.null'))) {
    ci.df['global'] <- GNCI(sample.mat = t(contrast %*% t(sample.global.null)),
                            alpha = alpha,
                            alternative = numToSide(truncation.parameters[ ,'alternative']),
                            truncation.parameters = truncation.parameters)
    colnames(ci.df$global) <- c('lower', 'upper')
  }

  ## Hybrid CI
  if (any(ci.type %in% 'hybrid')) {
    ### Different test.types
    if (any(test.type == 'symmetric')) {
      ci.df['hybrid.symmetric'] <- hybridCI(truncation.parameters = truncation.parameters,
                                            sample.global.null = sample.global.null,
                                            contrast = contrast,
                                            alpha = alpha / 2,
                                            test.type = 'symmetric',
                                            ci.type = ci.type,
                                            zscore.bound = zscore.bound)
    }
    if (any(test.type == 'UMPU')) {
      ci.df['hybrid.umpu']     <- hybridCI(truncation.parameters = truncation.parameters,
                                           sample.global.null = sample.global.null,
                                           contrast = contrast,
                                           alpha = alpha / 2,
                                           test.type = 'UMPU',
                                           ci.type = ci.type,
                                           zscore.bound = zscore.bound)
    }
  }
  ## Switch CI
  if (any(ci.type %in% 'switch')) {
    if (any(test.type == 'symmetric')) {
      ci.df['switch.symmetric'] <- switchCI(truncation.parameters = truncation.parameters,
                                            thres.list            = thres.list,
                                            sample.global.null    = sample.global.null,
                                            contrast              = contrast,
                                            alpha = alpha,
                                            test.type = 'symmetric',
                                            ci.type = ci.type,
                                            zscore.bound = zscore.bound)

    }
    if (any(test.type == 'UMPU')) {
      ci.df['switch.umpu']     <- switchCI(truncation.parameters = truncation.parameters,
                                           thres.list = thres.list,
                                           sample.global.null = sample.global.null,
                                           contrast = contrast,
                                           alpha = alpha,
                                           test.type = 'UMPU',
                                           ci.type = ci.type,
                                           zscore.bound = zscore.bound)
    }
  }
  return(ci.df)
}

### Finding rejection area
polyCI <- function(trunc.parameters, test.type = 'symmetric', alpha,
                   zscore.bound) {
  x     <- trunc.parameters[1]
  sd    <- sqrt(trunc.parameters[2])
  lower <- trunc.parameters[3]
  upper <- trunc.parameters[4]
  alternative <- numToSide(trunc.parameters[5])
  ### Condition
  minimal.zscore <- (sign(lower) == sign(upper)) * min(abs(c(lower, upper))) / sd
  ### No selection effect
  if ((lower == -Inf & upper == Inf) | minimal.zscore > zscore.bound | 'naive' %in% test.type) {
    ci.naive <- naiveCI(stat = x,
                        sd = sd,
                        alpha = alpha,
                        alternative = alternative)
    return(ci.naive)
  }
  ### Selection effect, symmetric p-values
  if ('symmetric' %in% test.type| alternative %in% c('less', 'greater')) {
    ci.trunc.symm <- polyCISymm(stat = x,
                                lower.trunc = lower,
                                upper.trunc = upper,
                                sd = sd,
                                alpha = alpha,
                                alternative = alternative,
                                sd.dist = zscore.bound)
    return(ci.trunc.symm)
  }
  ### Selection effect, UMPU p-values
  if (test.type == 'UMPU' & alternative == 'two.sided') {
    ci.umpu <- polyCIUMPU(stat = x,
                          lower.trunc = lower,
                          upper.trunc = upper,
                          sd = sd,
                          alpha = alpha,
                          sd.dist = zscore.bound)
    return(ci.umpu)
  }
}


polyCISymm <-  function(stat, lower.trunc, upper.trunc, sd, alpha, alternative, sd.dist) {
  if (alternative == 'two.sided') {
    return(c('lower' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, 1 - alpha / 2, sd.dist),
             'upper' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, alpha / 2, sd.dist)))
  }
  if (alternative == 'greater') {
    return(c('lower' = -Inf,
             'upper' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, alpha, sd.dist)))
  }
  if (alternative == 'less') {
    return(c('lower' = findCIoneSide(stat, lower.trunc, upper.trunc, sd, 1 - alpha, sd.dist),
             'upper' = Inf))
  }
}

naiveCI <- function(stat, sd, alpha, alternative) {
  if (alternative == 'two.sided') {
    return(stat + c('lower' = qnorm(alpha / 2, mean = 0, sd = sd),
                    'upper' = qnorm(1 - alpha / 2, mean = 0, sd = sd)))
  }
  if (alternative == 'greater') {
    return(stat + c('lower' = -Inf,
                    'upper' = qnorm(1 - alpha / 2, mean = 0, sd = sd)))
  }
  if (alternative == 'less') {
    return(stat + c('lower' = qnorm(alpha, mean = 0, sd = sd),
                    'upper' = Inf))
  }
}

### Global null confidence interval
GNCI <- function(sample.mat, alpha, alternative, truncation.parameters) {
  ci.bounds.complete <- matrix(NA, nrow = ncol(sample.mat), ncol = 2)

  ci.bounds.two   <- t(apply(sample.mat, 2, function(x)
    quantile(x, probs = c('lower' = alpha / 2,
                          'upper' = 1 - alpha / 2))))

  #### The one sided CI's return the minimal or the maximal value
  ci.bounds.less     <- apply(sample.mat, 2, function(x)
    quantile(x, probs = c('upper' = 1 - alpha)))
  ci.bounds.less    <- cbind('lower' = -Inf,
                             'upper' = ci.bounds.less)

  ci.bounds.greater <- apply(sample.mat, 2, function(x)
    quantile(x, probs = c('lower' = alpha)))
  ci.bounds.greater <- cbind('lower' = ci.bounds.greater,
                             'upper' = Inf)


  ci.bounds.complete[alternative == 'two.sided', ] <- ci.bounds.two[alternative     == 'two.sided', ]
  ci.bounds.complete[alternative == 'greater', ]   <- ci.bounds.greater[alternative == 'greater', ]
  ci.bounds.complete[alternative == 'less', ]      <- ci.bounds.less[alternative    == 'less', ]
  ### Using the bounds
  constrats.ci <- cbind('lower' = truncation.parameters[ ,1] - ci.bounds.complete[ ,2],
                        'upper' = truncation.parameters[ ,1] - ci.bounds.complete[ ,1])
  return(constrats.ci)
}

### Hybrid confidence interval
hybridCI <- function(truncation.parameters,
                     sample.global.null,
                     contrast,
                     alpha,
                     test.type,
                     ci.type,
                     zscore.bound) {
  ### Polyhedral CI (alpha / 2)
  poly.ci <- t(apply(truncation.parameters, 1, polyCI,
                     test.type = test.type,
                     alpha = alpha / 2,
                     zscore.bound = zscore.bound))
  colnames(poly.ci) <- c('lower', 'upper')
  ### Global-null CI (alpha / 2)
  global.ci <- GNCI(sample.mat = t(contrast %*% t(sample.global.null)),
                    alpha = alpha / 2,
                    alternative = numToSide(truncation.parameters[ ,'alternative']),
                    truncation.parameters = truncation.parameters)
  colnames(global.ci) <- c('lower', 'upper')
  ### Combining the CI's
  hybrid.ci <- cbind('lower' = pmax(poly.ci[,'lower'],
                                    global.ci[ ,'lower']),
                     'upper' = pmin(poly.ci[,'upper'],
                                    global.ci[ ,'upper']))
  return(hybrid.ci)
}



### Switch confidence interval
switchCI <- function(truncation.parameters,
                     thres.list,
                     sample.global.null,
                     contrast,
                     alpha,
                     test.type,
                     ci.type,
                     zscore.bound) {
  if (!thres.list$pass.second.threshold) {
    switch.ci <- hybridCI(truncation.parameters = truncation.parameters,
                          sample.global.null = sample.global.null,
                          contrast = contrast,
                          alpha = thres.list$new.alpha,
                          test.type = test.type,
                          ci.type = ci.type,
                          zscore.bound = zscore.bound)
  }
  if (thres.list$pass.second.threshold) {
    switch.ci <- t(apply(truncation.parameters, 1, polyCI,
                                  test.type = 'naive',
                                  alpha = alpha,
                                  zscore.bound))
  }
  colnames(switch.ci) <- c('lower', 'upper')
  return(switch.ci)
}
