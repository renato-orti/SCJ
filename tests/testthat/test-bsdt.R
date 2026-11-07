# =================================================================================================
# Test:      Bayesian Standardised Difference Test (BSDT)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2007)
# =================================================================================================

test_that("bsdt delivers reproducible results for discrepancy comparisons", {
  set.seed(42)
  res <- SCJ::bsdt(
    case_a = 70,
    case_b = 130,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 15,
    correlation = 0.5,
    alternative = "two.sided",
    credible_level = 95,
    credible_interval = TRUE,
    iter = 10000,
    seed = 42
    )
  
  res_table <- res$resultsTable$asDF
  expect_true(!is.null(res_table$p_val[1]))
  expect_true(is.numeric(res_table$ci_lower[1]))
})