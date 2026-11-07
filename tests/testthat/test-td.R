# =================================================================================================
# Test:      Test of Deficit (TD)
# Reference: Crawford, J. R., & Howell, D. C. (1998)
# =================================================================================================

test_that("td output matches published Crawford & Howell (1998) reference values", {
  res <- SCJ::td(
    case = 130,
    controls = 100,
    sd = 15,
    sample_size = 10,
    alternative = "two.sided",
    confidence_interval = TRUE,
    confidence_level = 95
  )
  
  res_table <- res$resultsTable$asDF
  expect_equal(res_table$t_val[1], 1.906925, tolerance = 1e-4)
  expect_equal(res_table$df[1], 9)
  expect_equal(res_table$p_val[1], 0.0889, tolerance = 1e-3)
  expect_equal(res_table$z_cc[1], 2.0, tolerance = 1e-5)
})