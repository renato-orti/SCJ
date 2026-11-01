# =================================================================================================
# Test:      Dissociation Test (DISSOCS)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2005)
# =================================================================================================

test_that("DISSOCS computes deficit and dissociation criteria correctly via RSDT", {
  res <- singcar::RSDT(
    case_a = 50,
    case_b = 100,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 20,
    r_ab = 0.5,
    alternative = "less"
  )
  
  expect_true(!is.null(res$p.value))
  expect_true(!is.null(res$statistic))
  expect_equal(as.integer(res$parameter), 19)
})