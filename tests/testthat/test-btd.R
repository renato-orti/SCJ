# =================================================================================================
# Test:      Bayesian Test of Deficit (BTD)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2007)
# =================================================================================================

test_that("btd delivers reproducible Bayesian estimates with fixed seed", {
  set.seed(42)
  res <- SCJ::btd(
    case = 130,
    controls = 100,
    sd = 15,
    sample_size = 10,
    alternative = "two.sided",
    credible_level = 95,
    credible_interval = TRUE,
    iter = 10000,
    seed = 42
    )
  
  res_table <- res$resultsTable$asDF
  expect_equal(res_table$z_cc[1], 2.0, tolerance = 1e-4)
  expect_true(res_table$ci_lower[1] < res_table$ci_upper[1])
  expect_equal(res_table$p_val[1], 0.0889, tolerance = 0.05)
})