# =================================================================================================
# Test:      Bayesian Dissociation Test (B-DISSOCS)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2007)
# =================================================================================================

test_that("B-DISSOCS evaluates Bayesian criteria reproducibly via BSDT", {
  set.seed(42)
  res <- singcar::BSDT(
    case_a = 60,
    case_b = 100,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 20,
    r_ab = 0.6,
    alternative = "less",
    int_level = 0.95,
    iter = 10000
  )
  
  expect_true(!is.null(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})
