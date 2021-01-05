context("Checking the switch CI for linear testing")

## Simulating data
p       <- 10
contrast <- diag(p)
a.vec    <- rep(1, p)
cov.mat  <- diag(p)

set.seed(959)
y               <- c(15, rnorm(p - 1))
first.threshold <- threshold <- threshold <- qchisq(0.95, p)
alpha           <- 0.05
test.val        <- t(a.vec) %*% y
test.side       <- rep('two.sided', nrow(contrast))

ci.type         <- 'switch'


zscore.bound    <- 5

test.var  <- drop(t(a.vec) %*% cov.mat %*% a.vec)
testthat::test_that('switch CI for linear testing returns expected values', {
  thres.list <- PSATinference:::thresholdWrapperLinear(test.val = test.val,
                                                       cov.mat  = cov.mat,
                                                       a.vec    = a.vec,
                                                       first.threshold.lower  = -Inf,
                                                       first.threshold.upper  = threshold,
                                                       alpha = alpha,
                                                       calc.t2 = TRUE)

  ### Finding the contrast values
  contrast.y     <- drop(a.vec %*% y)
  var.contrast.y <- (t(a.vec) %*% cov.mat %*% a.vec)
  ### Preparing data according to requested types
  trunc.para.mat <- apply(contrast, 1, PSATinference:::findTruncParametersLinear,
                          cov.mat = cov.mat, a.vec = a.vec,
                          threshold.lower = -Inf,
                          threshold.upper = threshold,
                          y = y)
  trunc.para.mat <- cbind(t(trunc.para.mat), 'alternative' = PSATinference:::sideToNum(test.side))
  samp.null         <- PSATinference:::sampleGNLinear(y = y,
                                                      null.mu = rep(0, p),
                                                      sigma   = cov.mat,
                                                      threshold.lower = -Inf,
                                                      threshold.upper = threshold,
                                                      a.vec = a.vec,
                                                      num.samp = 20000)
  switch.psat  <- PSATinference:::CIwrapper(truncation.parameters = trunc.para.mat,
                                            sample.global.null    = samp.null,
                                            contrast              = contrast,
                                            alpha                 = alpha,
                                            test.type             = 'symmetric',
                                            ci.type               = 'switch',
                                            thres.list            = thres.list,
                                            zscore.bound          = zscore.bound)
  testthat::expect_equal(colnames(switch.psat$switch.symmetric),
                         c('lower', 'upper'))
})
