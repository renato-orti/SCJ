# =================================================================================================
# Test:      Test of Deficit (TD)
# Reference: Crawford, J. R., & Howell, D. C. (1998)
# =================================================================================================

test_that("TD output matches published Crawford & Howell (1998) reference values", {
  res <- singcar::TD(
    case = 130,
    controls = 100,
    sd = 15,
    sample_size = 10,
    alternative = "two.sided",
    conf_int = TRUE,
    conf_level = 0.95,
    na.rm = TRUE
  )
  
  # Test statistic t = 1.906925, df = 9
  expect_equal(as.numeric(res$statistic), 1.906925, tolerance = 1e-4)
  expect_equal(as.integer(res$parameter), 9)
  
  # p-value (two-sided) ~ 0.0889
  expect_equal(as.numeric(res$p.value), 0.0889, tolerance = 1e-3)
  
  # Z-CC = (130 - 100) / 15 = 2.0
  expect_equal(as.numeric(res$estimate[1]), 2.0, tolerance = 1e-5)
})

test_that("TD rejects non-positive standard deviation", {
  expect_error(
    singcar::TD(case = 100, controls = 100, sd = 0, sample_size = 10),
    regexp = ".*"
  )
})