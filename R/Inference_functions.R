# Inference functions ------------------------------------------------------

pvalueWrapper <- function(truncation.parameters,
                          sample.global.null,
                          contrast, alpha,
                          pval.type,
                          test.type,
                          optim.control,
                          zscore.bound) {
  pval.df <- data.frame(matrix(NA, nrow = nrow(truncation.parameters), ncol = 0))
  if (any(pval.type %in% c('hybrid', 'polyhedral'))) {
    if (any(test.type == 'symmetric')) {
      pval.df['polyhedral.symmetric'] <- apply(truncation.parameters,
                                               1,
                                               polyPvalues,
                                               test.type = 'symmetric',
                                               zscore.bound = zscore.bound)
    }
    if (any(test.type == 'UMPU')) {
      pval.df['polyhedral.umpu'] <- apply(truncation.parameters,
                                           1,
                                           polyPvalues,
                                           test.type = 'UMPU',
                                           optim.control = optim.control,
                                           zscore.bound = zscore.bound)
    }
  }
  ## Naive
  if (any(pval.type == 'naive')) {
    pval.df['naive'] <- apply(truncation.parameters,
                              1,
                              polyPvalues,
                              test.type = 'naive',
                              zscore.bound = zscore.bound)
  }
  if (any(pval.type %in% c('hybrid', 'global.null'))) {
    contrast.samp <- sample.global.null %*% t(contrast)
    pval.df['global'] <- globalPvalues(contrast.val = truncation.parameters[ ,'x'],
                                       sample.mat   = contrast.samp,
                                       alternative  = numToSide(truncation.parameters[ ,'alternative']))
  }
  ## Hybrid p-value
  if (any(pval.type == 'hybrid')) {
    ### Different test.types
    if (any(test.type == 'symmetric')) {
      pval.df['hybrid.symmetric'] <- pmin(2 * apply(pval.df[c('global',
                                                              'polyhedral.symmetric')],
                                                    1, min), 1)
    }
    if (any(test.type == 'UMPU')) {
      pval.df['hybrid.umpu']     <- pmin(2 * apply(pval.df[c('global',
                                                              'polyhedral.umpu')],
                                                    1, min), 1)
    }
  }
  return(pval.df)
}


polyPvalues <- function(trunc.parameters, test.type = 'symmetric',
                        optim.control, alpha = 0.05,
                        zscore.bound) {
  x         <- trunc.parameters[1]
  sd        <- sqrt(trunc.parameters[2])
  lower     <- trunc.parameters[3]
  upper     <- trunc.parameters[4]
  alternative <- numToSide(trunc.parameters[5])
  ### Condition
  minimal.zscore <- (sign(lower) == sign(upper)) * min(abs(c(lower, upper))) / sd
  ### No selection effect
  if ((lower == -Inf & upper == Inf) | minimal.zscore > zscore.bound | test.type == 'naive') {
    naive.pval    <- naivePvalues(stat = x, sd = sd, alternative = alternative)
    return(naive.pval)
  }
  ### Selection effect, symmetric p-values
  if (test.type == 'symmetric' | alternative %in% c('less', 'greater')) {
      symmetric.pval <- polyPvaluesSymm(stat  = x,
                                        lower = lower,
                                        upper = upper,
                                        sd = sd,
                                        mean = 0,
                                        alternative = alternative)
      return(symmetric.pval)
    }
    ### Selection effect, UMPU p-values
  if (test.type == 'UMPU' & alternative == 'two.sided') {
      umpu.pval   <- findUMPU(test.stat = x,
                                alpha = NA,
                                lower.bound = lower,
                                upper.bound = upper,
                                mean = 0,
                                sd = sd,
                                optim.type = 'pvalue',
                                optim.control = optim.control)
    umpu.pval <- umpu.pval$optim$bestmem
    return(umpu.pval[1])
    }
  }


polyPvaluesSymm <- function(stat, lower, upper, mean, sd, alternative) {
  if (alternative == 'two.sided') {
    symmetric.pval <- calcProbDoubleTrunc(stat, lower, upper, mean = 0, sd = sd)
    symmetric.pval <- pmin(2 * min(symmetric.pval, 1 - symmetric.pval), 1)
  }
  if (alternative == 'greater') {
    symmetric.pval <- calcProbDoubleTrunc(stat, lower, upper, mean = 0, sd = sd)
    symmetric.pval <- 1 - symmetric.pval
  }
  if (alternative == 'less') {
    symmetric.pval <- calcProbDoubleTrunc(stat, lower, upper, mean = 0, sd = sd)
    symmetric.pval <- symmetric.pval
  }
  return(symmetric.pval)
}

### Naive p-values based on normal distribution two sided using same parameters as polyPvalues
naivePvalues <- function(stat, mean = 0, sd, alternative) {
  contrast.pval <- pnorm(stat,
                         mean = 0,
                         sd = sd)
  if (alternative == 'two.sided') {
    contrast.pval <- pmin(2 * pmin(contrast.pval, 1 - contrast.pval), 1)
  }
  if (alternative == 'greater') {
    contrast.pval <- 1 - pnorm(stat, mean = mean, sd = sd)
  }
  if (alternative == 'less') {
    contrast.pval <- pnorm(stat, mean = mean, sd = sd)
  }
  return(contrast.pval)
}

### Global p-values based on sampling from the null distribution
globalPvalues <- function(contrast.val, sample.mat, alternative) {
  compare.mat   <- sweep(sample.mat, 2, contrast.val, '-')
  num.samp      <- nrow(sample.mat)
  global.pvalue <- 2 * apply(cbind(colMeans(compare.mat < 0) + 1 / num.samp,
                                   colMeans(compare.mat > 0) + 1 / num.samp), 1, min)
  global.pvalue[alternative == 'greater'] <- (colMeans(compare.mat > 0) + 1 / num.samp)[alternative ==
                                                                                        'greater']
  global.pvalue[alternative == 'less'] <- (colMeans(compare.mat > 0) + 1 / num.samp)[alternative ==
                                                                                     'less']
  return(global.pvalue)
}
