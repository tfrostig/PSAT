context('Testing UMUVE p-value finder')

### Parameters are taken such that error is returned due to too large of a truncation.

test_that('Error for too large values of truncation is known', {
  expect_error(PSATinference:::findUMPU(test.stat   = 30,
                                         lower.bound = -11.15088,
                                         upper.bound = 11.15088,
                                         mean        = 0,
                                         sd          = 1,
                                         optim.type  = 'pvalue',
                                         optim.control = DEoptim::DEoptim.control()))
})
