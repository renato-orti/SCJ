# =================================================================================================
# Module:    DISSOCS - Single Case Comparison
# File:      dissocs.b.R (Back-end Analysis Script)
# Purpose:   R computational logic for executing single-case deficit tests 
#            using the singcar package (Rittmo & McIntosh, 2021) in jamovi.
#            This script processes UI inputs, runs the statistical functions, 
#            and populates the results table.
# 
# Authors:
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
#     simulations and revised tests for dissociations. Neuropsychology, 19, 318-331.
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

dissocsClass <- if (
  requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "dissocsClass",
    inherit = dissocsBase,
    private = list(
      .run = function() {
        
        # Extract UI input data
        case_a      <- self$options$case_a            # Single Case Score Task A
        case_b      <- self$options$case_b            # Single Case Score Task B
        controls_a  <- self$options$controls_a        # Control Mean Task A
        controls_b  <- self$options$controls_b        # Control Mean Task B
        sd_a        <- self$options$sd_a              # Control SD Task A
        sd_b        <- self$options$sd_b              # Control SD Task B
        sample_size <- self$options$sample_size       # Control Sample Size
        correlation <- self$options$correlation       # Task A-B Correlation in Control Sample
        
        alt_a <- if (!is.null(self$options$alternative_a)) self$options$alternative_a else "less"
        alt_b <- if (!is.null(self$options$alternative_b)) self$options$alternative_b else "less"
        
      # Silent check for empty inputs on initial run
        if (is.null(case_a) || is.null(case_b) ||
            is.null(controls_a) || is.null(controls_b) ||
            is.null(sd_a) || is.null(sd_b) ||
            is.null(sample_size) || is.null(correlation) ||
            is.na(case_a) || is.na(case_b) ||
            is.na(controls_a) || is.na(controls_b) ||
            is.na(sd_a) || is.na(sd_b) ||
            is.na(sample_size) || is.na(correlation)
        ) {
          return()
        }
        
      # Input Range Validation with user-friendly rejection messages
        if (sd_a <= 0 || sd_b <= 0) {
          jmvcore::reject("Standard deviations must be strictly positive (> 0).")
          return()
        }
        if (sample_size <= 1) {
          jmvcore::reject("Sample size must be an integer greater than 1.")
          return()
        }
        if (correlation <= -1 || correlation >= 1) {
          jmvcore::reject("Correlation coefficient (r) must be strictly between -1 and 1.")
          return()
        }
        
      # Execute Revised Standardised Difference Test (RSDT) and Test of Deficit (TD)
        res_rsdt <- singcar::RSDT(
          case_a = case_a, case_b = case_b,
          controls_a = controls_a, controls_b = controls_b, 
          sd_a = sd_a, sd_b = sd_b,
          sample_size = sample_size, r_ab = correlation, 
          alternative = "two.sided", na.rm = TRUE
        )
        
        res_td_a <- singcar::TD(
          case = case_a, controls = controls_a, sd = sd_a, 
          sample_size = sample_size, alternative = alt_a, na.rm = TRUE
        )
        
        res_td_b <- singcar::TD(
          case = case_b, controls = controls_b, sd = sd_b, 
          sample_size = sample_size, alternative = alt_b, na.rm = TRUE
        )
        
      # Populate results tables and generate report
        private$.populateTables(res_rsdt, res_td_a, res_td_b, alt_a, alt_b)
        
      },
      
      .populateTables = function(res_rsdt, res_td_a, res_td_b, alt_a, alt_b) {
        
      # Extract parameters, p-values, and estimates
        df_val      <- as.integer(res_rsdt$parameter)
        pvalue_ab   <- as.numeric(res_rsdt$p.value)
        pvalue_a    <- as.numeric(res_td_a$p.value)
        pvalue_b    <- as.numeric(res_td_b$p.value)
        z_cc_a      <- as.numeric(res_rsdt$estimate[1])
        z_cc_b      <- as.numeric(res_rsdt$estimate[2])
        z_dcc       <- as.numeric(res_rsdt$estimate[3])
        point_estim <- as.numeric(res_rsdt$estimate[4])
        
      # Perform dissociation classification 
        class_res <- classifyDissociation(
          p_a = pvalue_a,
          p_b = pvalue_b,
          p_ab = pvalue_ab,
          alpha = 0.05
        )
        
      # Populate Clinical Classification Table
        table_class <- self$results$classTable
        table_class$setRow(
          rowNo = 1, 
          values = list(
            score_a = class_res$score_a,
            score_b = class_res$score_b,
            dissoc_type = class_res$dissoc_type
          )
        )
        
        dir_note_a <- if (alt_a == "less") "Case < Controls" else "Case > Controls"
        dir_note_b <- if (alt_b == "less") "Case < Controls" else "Case > Controls"
        table_class$setNote(
          "criteria_info",
          paste0("The difference test is two-sided. ",
                 "Deficit criteria evaluated directional (Task A: ", dir_note_a, 
                 "; Task B: ", dir_note_b, "; alpha = .05). Theoretical criteria: Crawford & Garthwaite (2005).")
        )
        
      # Populate Statistical Results Table
        table_stats <- self$results$resultsTable
        table_stats$setRow(
          rowNo = 1, 
          values = list(
            z_cc_a      = z_cc_a,
            p_val_a     = pvalue_a,
            z_cc_b      = z_cc_b,
            p_val_b     = pvalue_b,
            z_dcc       = z_dcc,
            df          = df_val,
            p_val_ab    = pvalue_ab,
            point_estim = point_estim
          )
        )
        
      # Footnote: Hypothesis and Theoretical Criteria Details
        table_stats$setNote(
          "info",
          paste0("Hypothesis: Task A (", dir_note_a, "), Task B (", dir_note_b, "), Difference (Case ≠ Controls).", 
                 "Theoretical criteria: Revised Standardised Difference Test (RSDT; Crawford & Garthwaite, 2005)."))
        },
      
      .plot = function(image, ggtheme, theme, ...) {
        case_a     <- self$options$case_a
        case_b     <- self$options$case_b
        controls_a <- self$options$controls_a
        controls_b <- self$options$controls_b
        sd_a       <- self$options$sd_a
        sd_b       <- self$options$sd_b
        
        if (is.null(case_a) || is.null(case_b) || 
            is.null(controls_a) || is.null(controls_b) || 
            is.null(sd_a) || is.null(sd_b) ||
            is.na(case_a) || is.na(case_b) || 
            is.na(controls_a) || is.na(controls_b) || 
            is.na(sd_a) || is.na(sd_b) || 
            sd_a <= 0 || sd_b <= 0) {
          return()
        }
        
        z_a <- (case_a - controls_a) / sd_a
        z_b <- (case_b - controls_b) / sd_b
        
        df_plot <- data.frame(
          Task = factor(c("Task A", "Task B", "Task A", "Task B"), levels = c("Task A", "Task B")),
          Group = factor(c("Case", "Case", "Control Mean", "Control Mean"), levels = c("Control Mean", "Case")),
          Z_Score = c(z_a, z_b, 0, 0)
        )
        
        profile_plot <- ggplot2::ggplot(
          df_plot, 
          ggplot2::aes(x = Task, y = Z_Score, group = Group, color = Group)) +
          ggplot2::geom_line(linewidth = 1.2) + 
          ggplot2::geom_point(size = 4) +
          ggplot2::scale_color_manual(values = c("Case" = "#d95f02", "Control Mean" = "#8c96c6")) +
          ggplot2::labs(
            x = "Tasks",
            y = "Standardised Score (Z-Score)"
          ) +
          ggtheme +
          ggplot2::theme(
            legend.position = "bottom",
            legend.title = ggplot2::element_blank()
          )
        
        print(profile_plot)
        TRUE
      }
    )
  )