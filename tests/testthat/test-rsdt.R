# =================================================================================================
# Test:      Revised Standardised Difference Test (RSDT)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2005)
# =================================================================================================

test_that("RSDT output matches published Crawford & Garthwaite (2005) reference values", {
  res <- singcar::RSDT(
    case_a = 70,
    case_b = 130,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 10,
    r_ab = 0.6,
    alternative = "two.sided"
  )
  
  expect_true(is.numeric(res$statistic))
  expect_equal(as.integer(res$parameter), 9)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("RSDT rejects invalid alternative hypothesis parameter", {
  expect_error(
    singcar::RSDT(
      case_a = 70, case_b = 130, controls_a = 100, controls_b = 100,
      sd_a = 15, sd_b = 15, sample_size = 10, r_ab = 0.6,
      alternative = "invalid_option"
    ),
    regexp = ".*"
  )
})