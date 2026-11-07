# =================================================================================================
# Test:      Revised Standardised Difference Test (RSDT)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2005)
# =================================================================================================

test_that("rsdt output matches expected values", {
  res <- SCJ::rsdt(
    case_a = 70,
    case_b = 130,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 10,
    correlation = 0.6,
    alternative = "two.sided"
  )
  
  res_table <- res$resultsTable$asDF
  expect_true(is.numeric(res_table$t_val[1]))
  expect_equal(res_table$df[1], 9)
  expect_true(res_table$p_val[1] >= 0 && res_table$p_val[1] <= 1)
})