# =================================================================================================
# Test:      Bayesian Dissociation Test (B-DISSOCS)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2007)
# =================================================================================================

test_that("bdissocs evaluates Bayesian criteria reproducibly", {
  set.seed(42)
  res <- SCJ::bdissocs(
    case_a = 60,
    case_b = 100,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 20,
    correlation = 0.6,
    alternative_a = "less",
    alternative_b = "less",
    credible_level = 95,
    credible_interval = TRUE,
    iter = 10000,
    seed = 42
    )
  
  res_table <- res$resultsTable$asDF
  expect_true(!is.null(res_table$p_val_ab[1]))
  expect_true(res_table$p_val_ab[1] >= 0 && res_table$p_val_ab[1] <= 1)
})