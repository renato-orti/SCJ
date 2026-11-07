# =================================================================================================
# Test:      Dissociation Test (DISSOCS)
# Reference: Crawford, J. R., & Garthwaite, P. H. (2005)
# =================================================================================================

test_that("classifyDissociation assigns correct profiles", {
  # No deficit
  expect_equal(SCJ:::classifyDissociation(p_a = 0.10, p_b = 0.20, p_ab = 0.50)$dissoc_type, "None")
  
  # Single deficit but with no significant difference (No dissociation)
  expect_equal(SCJ:::classifyDissociation(p_a = 0.01, p_b = 0.50, p_ab = 0.10)$dissoc_type, "None")
  
  # Single deficit but with significant difference (Classical dissociation)
  expect_equal(SCJ:::classifyDissociation(p_a = 0.01, p_b = 0.50, p_ab = 0.01)$dissoc_type, "Classical Dissociation")
  
  # Double deficits on both tasks with significant differences (Strong dissociation)
  expect_equal(SCJ:::classifyDissociation(p_a = 0.01, p_b = 0.01, p_ab = 0.01)$dissoc_type, "Strong Dissociation")
})

test_that("dissocs computes deficit and dissociation criteria correctly", {
  res <- SCJ::dissocs(
    case_a = 50,
    case_b = 100,
    controls_a = 100,
    controls_b = 100,
    sd_a = 15,
    sd_b = 15,
    sample_size = 20,
    correlation = 0.5,
    alternative_a = "less",
    alternative_b = "less"
  )
  
  res_table <- res$resultsTable$asDF
  expect_true(!is.null(res_table$p_val_ab[1]))
  expect_equal(res_table$df[1], 19)
})