# =================================================================================================
# Test:      Bayesian Standardised Difference Test (BSDT)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2007)
# =================================================================================================

test_that("BSDT delivers reproducible results for discrepancy comparisons", {
  set.seed(42)
  res <- singcar::BSDT(
    case_a = 70,
    case_b = 130,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 15,
    r_ab = 0.5,
    alternative = "two.sided",
    int_level = 0.95,
    iter = 10000
  )
  
  expect_true(!is.null(res$p.value))
  expect_true(is.numeric(res$interval))
  expect_equal(length(res$interval), 5) # Level + Z-DCC (Lower, Upper) + Prop (Lower, Upper)
})