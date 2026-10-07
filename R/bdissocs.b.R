# =================================================================================================
# Module:    DISSOCS - Single Case Comparison
# File:      dissocs.b.R (Back-end Analysis Script)
# Purpose:   R computational logic for executing single-case deficit tests 
#            using the singcar package (Rittmo & McIntosh, 2021) in Jamovi.
#            This script processes UI inputs, runs the statistical functions 
#            and populates the results table.
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

bdissocsClass <- if (
  requireNamespace('jmvcore',
                   quietly=TRUE)) R6::R6Class(
    "bdissocsClass",
    inherit = bdissocsBase,
    private = list(
        .run = function() {

          # Look for UI input data
          case_a <- self$options$case_a
          case_b <- self$options$case_b
          controls_a <- self$options$controls_a
          controls_b <- self$options$controls_b
          sd_a <- self$options$sd_a
          sd_b <- self$options$sd_b
          sample_size <- self$options$sample_size
          correlation <- self$options$correlation
          credible_interval <- self$options$credible_interval
          credible_level <- self$options$credible_level
          iter <- self$options$iter
          unstandardised <- self$options$unstandardised
          
          # Check for not-empty UI input data
          if (is.null(case_a) || is.null(case_b) ||
              is.null(controls_a) || is.null(controls_b) ||
              is.null(sd_a) || is.null(sd_b) ||
              is.null(sample_size) ||
              is.null(correlation) ||
              is.null(credible_level) ||
              is.null(iter)
          )
            return()
          
          # Silent JAMOVI warnings waiting for UI input data
          if (
            length(case_a) > 0 && !is.na(case_a) &&
            length(case_b) > 0 && !is.na(case_b) &&
            length(controls_a) > 0 && !is.na(controls_a) &&
            length(controls_b) > 0 && !is.na(controls_b) &&
            length(sd_a) > 0 && !is.na(sd_a) &&
            length(sd_b) > 0 && !is.na(sd_b) &&
            length(sample_size) > 0 && !is.na(sample_size) && sample_size > 3 &&
            length(correlation) > 0 && !is.na(correlation) && correlation != 0 && 
            correlation >= -1 && correlation <= 1 &&
            length(credible_level) > 0 && !is.na(credible_level) &&
            length(iter) > 0 && !is.na(iter)
          ) {
            
            # Calculate Standardized Z-Scores and Dynamic Direction
            z_a <- (case_a - controls_a) / sd_a
            z_b <- (case_b - controls_b) / sd_b
            
            direction_ab <- ifelse(z_a < z_b, "less", "greater")
            direction_a <- ifelse(z_a < 0, "less", "greater")
            direction_b <- ifelse(z_b < 0, "less", "greater")
            
            # Execute singcar Statistical Tests
            res_bsdt <- singcar::BSDT(
              case_a = case_a,
              case_b = case_b,
              controls_a = controls_a,
              controls_b = controls_b,
              sd_a = sd_a,
              sd_b = sd_b,
              sample_size = sample_size,
              r_ab = correlation,
              int_level = credible_level/100,
              iter = iter,
              alternative = direction_ab,
              unstandardised = unstandardised,
              calibrated = TRUE,
              na.rm = TRUE  
            )
            
            res_btd_a <- singcar::BTD(
              case = case_a, 
              controls = controls_a,
              sd = sd_a,
              sample_size = sample_size,
              alternative = direction_a,
              int_level = credible_level/100,
              iter = iter,
              na.rm = TRUE)
            
            res_btd_b <- singcar::BTD(
              case = case_b, 
              controls = controls_b,
              sd = sd_b,
              sample_size = sample_size,
              alternative = direction_b,
              int_level = credible_level/100,
              iter = iter,
              na.rm = TRUE)
            
            # Extract Values
            df_val <- as.integer(res_bsdt$parameter)
            
            pvalue_ab <- as.numeric(res_bsdt$p.value)
            pvalue_a <- as.numeric(res_btd_a$p.value)
            pvalue_b <- as.numeric(res_btd_b$p.value)
            
            z_cc_a <- as.numeric(res_bsdt$estimate[1])
            z_cc_b <- as.numeric(res_bsdt$estimate[2])
            z_dcc <- as.numeric(res_bsdt$estimate[3])
            point_estim <- as.numeric(res_bsdt$estimate[4])
            
            # Dissociation Logic Tree for Classification
            score_a_bool <- pvalue_a < 0.05
            score_b_bool <- pvalue_b < 0.05
            dissoc_bool <- pvalue_ab <= 0.05
            
            score_a_text <- ifelse(score_a_bool, 
                                     "Out of normal range",
                                     "Normal Range")
            
            score_b_text <- ifelse(score_b_bool,
                                     "Out of normal range",
                                     "Normal Range")
            
            if (dissoc_bool) {
              if (score_a_bool && score_b_bool) {
                dissoc_type_text <- "Strong Dissociation"
              } else if (score_a_bool || score_b_bool) {
                dissoc_type_text <- "Classical Dissociation"
              } else {
                dissoc_type_text <- "None"
              }
            } else {
              dissoc_type_text <- "None"
            }
            
            # Populate Dissociation Classification Table
            table_class <- self$results$classTable
            table_class$setRow(
              rowNo = 1, 
              values = list(
                score_a = score_a_text,
                score_b = score_b_text,
                dissoc_type = dissoc_type_text
              )
            )
            
            # Populate Statistical Results Table
            table_stats <- self$results$resultsTable
            table_stats$setRow(
              rowNo = 1, 
              values = list(
                iter = iter,
                df = df_val,
                z_cc_a = z_cc_a,
                p_val_a = pvalue_a,
                z_cc_b = z_cc_b,
                p_val_b = pvalue_b,
                z_dcc = z_dcc,
                p_val_ab = pvalue_ab,
                point_estim = point_estim
              )
            )
            
            # Text Output Compilation
            out_bsdt <- capture.output(res_bsdt)
            out_bstd_a <- capture.output(res_btd_a)
            out_bstd_b <- capture.output(res_btd_b)
            
            report <- paste(
              "--- TEST DIFFERENCES BETWEEN TASKS (BSDT) ---", 
              paste(out_bsdt, collapse = "\n"),
              "\n--- TEST DIFFERENCE CASE CONTROL TASK A (BTD) ---",
              paste(out_bstd_a, collapse = "\n"),
              "\n--- TEST DIFFERENCE CASE CONTROL TASK B (BTD) ---",
              paste(out_bstd_b, collapse = "\n"),
              sep = "\n"
            )
            self$results$txtOutput$setContent(report)
          }
        },
        
        .plot = function(image, ggtheme, theme, ...) {
          
          case_a <- self$options$case_a
          case_b <- self$options$case_b
          controls_a <- self$options$controls_a
          controls_b <- self$options$controls_b
          sd_a <- self$options$sd_a
          sd_b <- self$options$sd_b
          iter <- self$options$iter
          
          if (is.null(case_a) || is.null(case_b) || 
              is.null(controls_a) || is.null(controls_b) || 
              is.null(sd_a) || is.null(sd_b) ||
              is.na(case_a) || is.na(case_b) || 
              is.na(controls_a) || is.na(controls_b) || 
              is.na(sd_a) || is.na(sd_b) || 
              sd_a <= 0 || sd_b <= 0)
            return()
          
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
              title = "Dissociation Profile: Case vs Controls",
              x = "Tasks",
              y = "Standardized Score (Z-Score)",
              caption = paste("Bayesian Iterations:", iter)
            ) +
            ggtheme +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 16),
              plot.caption = ggplot2::element_text(face = "italic", color = "gray40", size = 11),
              axis.title = ggplot2::element_text(size = 14),
              legend.position = "bottom",
              legend.title = ggplot2::element_blank(),
              legend.text = ggplot2::element_text(size = 12)
            )
          
          print(profile_plot)
          TRUE
          }
        )
    )