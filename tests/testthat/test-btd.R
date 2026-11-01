# =================================================================================================
# Test:      Bayesian Test of Deficit (BTD)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2007)
# =================================================================================================

test_that("BTD delivers reproducible Bayesian estimates with fixed seed", {
  set.seed(42)
  res <- singcar::BTD(
    case = 130,
    controls = 100,
    sd = 15,
    sample_size = 10,
    alternative = "two.sided",
    int_level = 0.95,
    iter = 10000,
    na.rm = TRUE
  )
  
  expect_equal(as.numeric(res$estimate[1]), 2.0, tolerance = 1e-4) # Z-CC
  expect_true(res$interval[2] < res$interval[3])                    # Lower CrI < Upper CrI
  
  # Monte Carlo p-value converges toward TD theoretical value ~0.0889
  expect_equal(as.numeric(res$p.value), 0.0889, tolerance = 0.05)
})