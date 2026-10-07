# =================================================================================================
# Module:    RSDT - Single Case Comparison
# File:      rsdt.b.R (Back-end Analysis Script)
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
#     Crawford, J. R., & Garthwaite, P. H. (2005). Testing for Suspected Impairments and Dissociations
#     in Single-Case Studies in Neuropsychology: Evaluation of Alternatives Using Monte Carlo
#     Simulations and Revised Tests for Dissociations. Neuropsychology, 19(3), 318 - 331.
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

rsdtClass <- if (
  requireNamespace('jmvcore', 
                   quietly=TRUE)) R6::R6Class(
    "rsdtClass",
    inherit = rsdtBase,
    private = list(
        .run = function() {
          
          # Look for UI input data
          case_a <- self$options$case_a                     # Single Case Score Task A
          case_b <- self$options$case_b                     # Single Case Score Task B
          controls_a <- self$options$controls_a             # Control Score Task A
          controls_b <- self$options$controls_b             # Control Score Task B
          sd_a <- self$options$sd_a                         # Standard deviation Task A from Control Sample
          sd_b <- self$options$sd_b                         # Standard deviation Task B from Control Sample
          sample_size <- self$options$sample_size           # Sample Size
          correlation <- self$options$correlation           # Correlation in control sample between Task A and B
          alternative <- self$options$alternative           # Alternative Hypothesis
          
          # Check for not-empty UI input data
          if (is.null(case_a) ||
              is.null(case_b) ||
              is.null(controls_a) ||
              is.null(controls_b) ||
              is.null(sd_a) ||
              is.null(sd_b) ||
              is.null(sample_size) ||
              is.null(correlation)
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
            length(sample_size) > 0 && !is.na(sample_size) && sample_size > 1 &&
            length(correlation) > 0 && !is.na(correlation)
            )
            
            {
            
          # Calculate Test of Deficit (TD) with singcar package
            res <- singcar::RSDT(
              case_a = case_a,
              case_b = case_b,
              controls_a = controls_a,
              controls_b = controls_b,
              sd_a = sd_a,
              sd_b = sd_b,
              sample_size = sample_size,
              r_ab = correlation,
              alternative = alternative,
              na.rm = TRUE)
            
            # Extract Z-CC, Z-DCC and Point estimate of the proportion
            z_cc_a <- as.numeric(res$estimate[1])
            z_cc_b <- as.numeric(res$estimate[2])
            z_dcc <- as.numeric(res$estimate[3])
            point_estim <- as.numeric(res$estimate[4])
            
            # Set Variables Type (numeric or integer)
            alt_out <- as.character(
              switch(
                alternative,
                "two.sided" = "Two Sided (Case ≠ Control)",
                "less"      = "Less (Case < Control)",
                "greater"   = "Greater (Case > Control)"))
            t_stat <- as.numeric(res$statistic)
            df_val <- as.integer(res$parameter)
            p_val <- as.numeric(res$p.value)
            
            # Table Output
            table <- self$results$resultsTable
            
            # List RSDT Output values
            table$setRow(
              rowNo = 1, 
              values = list(
                alternative = alt_out,
                t_val = t_stat,
                df = df_val,
                p_val = p_val,
                z_cc_a = z_cc_a,
                z_cc_b = z_cc_b,
                z_dcc = z_dcc,
                point_estim = point_estim)
              )
            
            # Text Output
            txt_output <- paste(capture.output(res), collapse = "\n")
            self$results$txtOutput$setContent(txt_output)
          }
          },
          
          .plot = function(image, ggtheme, theme, ...) {
            
            # Look for UI input data
            case_a <- self$options$case_a
            case_b <- self$options$case_b
            controls_a <- self$options$controls_a
            controls_b <- self$options$controls_b
            sd_a <- self$options$sd_a
            sd_b <- self$options$sd_b
            
            # Check for not-empty UI input data
            if (is.null(case_a) || is.null(case_b) || 
                is.null(controls_a) || is.null(controls_b) || 
                is.null(sd_a) || is.null(sd_b) ||
                is.na(case_a) || is.na(case_b) || 
                is.na(controls_a) || is.na(controls_b) || 
                is.na(sd_a) || is.na(sd_b) || 
                sd_a <= 0 || sd_b <= 0)
              return()
            
            # Trasform raw scores in z-scores
            z_a <- (case_a - controls_a) / sd_a
            z_b <- (case_b - controls_b) / sd_b
            
            # Dataset for plot
            df_plot <- data.frame(
              Task = factor(c("Task A", "Task B", "Task A", "Task B"), levels = c("Task A", "Task B")),
              Group = factor(c("Case", "Case", "Control Mean", "Control Mean"), levels = c("Control Mean", "Case")),
              Z_Score = c(z_a, z_b, 0, 0)
            )
            
            # Design plot with ggplot package
            profile_plot <- ggplot2::ggplot(
              df_plot, 
              ggplot2::aes(x = Task, y = Z_Score, group = Group, color = Group)) +
              
              ggplot2::geom_line(linewidth = 1.2) + 
              ggplot2::geom_point(size = 4) +
              ggplot2::scale_color_manual(values = c("Case" = "#d95f02", "Control Mean" = "#8c96c6")) +
              
              ggplot2::labs(
                title = "Dissociation Profile: Case vs Controls",
                x = "Tasks",
                y = "Standardized Score (Z-Score)"
              ) +
              ggtheme +
              ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 16),
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