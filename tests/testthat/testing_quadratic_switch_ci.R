context("Checking the switch CI")

## Simulating data
p       <- 10
contrast <- diag(p)
K        <- diag(p)
cov.mat  <- diag(p)

set.seed(959)
y               <- c(15, rnorm(p - 1))
first.threshold <- threshold <- threshold <- qchisq(0.95, p)
alpha           <- 0.05
test.val        <- t(y) %*% y
test.side       <- rep('two.sided', nrow(contrast))

ci.type         <- 'switch'
zscore.bound    <- 5

testthat::test_that('switch CI for quadratic testing returns expected values', {
  thres.list <- PSATinference:::thresholdWrapperQuadratic(test.val = test.val,
                                                          cov.mat = cov.mat,
                                                          K = K,
                                                          first.threshold = threshold,
                                                          second.threshold = NULL, ## Should add to parameters?
                                                          alpha = alpha,
                                                          calc.t2 = 'switch' %in% ci.type)
  ### Finding the contrast values
  contrast.y     <- drop(contrast %*% y)
  var.contrast.y <- diag(contrast %*% cov.mat %*% t(contrast))
  ### Preparing data according to requested types
  trunc.para.mat <- apply(contrast, 1, PSATinference:::findTruncParametersQuardratic,
                          cov.mat = cov.mat, K = K, threshold = threshold, y = y)
  trunc.para.mat <- cbind(t(trunc.para.mat), 'alternative' = PSATinference:::sideToNum(test.side))
  samp.null         <- PSATinference:::sampleGNQuadratic(y = y,
                                                         null.mu = rep(0, p),
                                                         sigma   = cov.mat,
                                                         threshold = threshold,
                                                         K.mat = K,
                                                         num.samp = 10000,
                                                         burn.in  = 200)
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
