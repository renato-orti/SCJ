# =================================================================================================
# Module: DISSOCS - Single Case Comparison
# File: dissociationBranch.R (Back-end Analysis Script)
# Purpose: Shared decision tree logic for classifying classical vs. strong dissociations
# 
# Authors:
#
#   - Michela Possenti
#     Affiliation: University of Campania "L. Vanvitelli"
#     ORCID: https://orcid.org/0000-0001-8937-0611 | Scopus ID: 57215597427
#
#   - Renato Orti
#     Affiliation: University of Padova
#     ORCID: https://orcid.org/0000-0002-6856-7981 | Scopus ID: 57218222421
#
# References:
#   - Crawford, J. R. & Garthwaite, P.H. (2005). Testing for suspected impairments and dissociations 
#     in single-case studies in neuropsychology: Evaluation of alternatives using Monte Carlo 
#     simulations and revised tests for dissociations”.  Neuropsychology, 19, 318-331.
#     https://doi.org/10.1037/0894-4105.19.3.318
#
#   - Rittmo J.Ö., & McIntosh R.D. (2021). singcar: Comparing single cases 
#     to small samples in R. Journal of Open Source Software, 6(68), 3887. 
#     https://doi.org/10.21105/joss.03887
#
# Note:      Developed to provide a zero-code graphical user interface, 
#            making single-case comparisons easily accessible for clinical 
#            and experimental neuropsychology without compromising scientific rigor.
# =================================================================================================

classifyDissociation <- function(p_a, p_b, p_ab, alpha = 0.05) {
  score_a_bool <- p_a <= alpha
  score_b_bool <- p_b <= alpha
  dissoc_bool  <- p_ab <= alpha
  
  score_a_text <- ifelse(score_a_bool, "Out of normal range", "Normal Range")
  score_b_text <- ifelse(score_b_bool, "Out of normal range", "Normal Range")
  
  if (dissoc_bool) {
    if (score_a_bool && score_b_bool) {
      dissoc_type <- "Strong Dissociation"
    } else if (score_a_bool || score_b_bool) {
      dissoc_type <- "Classical Dissociation"
    } else {
      dissoc_type <- "None"
    }
  } else {
    dissoc_type <- "None"
  }
  
  return(list(
    score_a = score_a_text,
    score_b = score_b_text,
    dissoc_type = dissoc_type
  ))
  }
