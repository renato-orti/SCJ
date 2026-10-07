# =================================================================================================
# Module:    TD - Single Case Comparison
# File:      td.b.R (Back-end Analysis Script)
# Purpose:   R computational logic for executing single-case deficit tests 
#            using the singcar package (Rittmo & McIntosh, 2021) in Jamovi.
#            This script processes UI inputs, runs the statistical functions, 
#            formats the confidence intervals, and populates the results table.
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
#   - Crawford, J. R., & Howell, D. C. (1998). Comparing an Individual's Test 
#     Score Against Norms Derived from Small Samples. The Clinical 
#     Neuropsychologist, 12(4), 482-486.
#     https://doi.org/10.1076/clin.12.4.482.7241
#
#   - Rittmo J.Ö., & McIntosh R.D. (2021). singcar: Comparing single cases 
#     to small samples in R. Journal of Open Source Software, 6(68), 3887. 
#     https://doi.org/10.21105/joss.03887
#
# Note:      Developed to provide a zero-code graphical user interface, 
#            making single-case comparisons easily accessible for clinical 
#            and experimental neuropsychology without compromising scientific rigor.
# =================================================================================================

tdClass <- if (
  requireNamespace('jmvcore', 
                   quietly=TRUE)) 
  R6::R6Class(
    "tdClass",
    inherit = tdBase,
    private = list(
        .run = function() {
          
          # Look for UI input data
          case <- self$options$case                     # Single Case Score Stored in "Case"
          controls <- self$options$controls             # Control Score Stored in "Controls"
          sd <- self$options$sd                         # Standard deviation from Control Sample
          sample_size <- self$options$sample_size       # Sample Size
          alternative <- self$options$alternative       # Alternative Hypothesis
          conf_int <- self$options$confidence_interval  # Check if user needs for CI
          conf_lev <- self$options$confidence_level     # Check for user confidence level (default 95%)
          
          # Check for not-empty UI input data
          if (is.null(case) ||
              is.null(controls) ||
              is.null(sd) ||
              is.null(sample_size)
              )
            return()
          
          # Silent JAMOVI warnings waiting for UI input data
          if (
              length(case) > 0 && !is.na(case) &&
              length(controls) > 0 && !is.na(controls) &&
              length(sd) > 0 && !is.na(sd) &&
              length(sample_size) > 0 && !is.na(sample_size) && sample_size > 1
              )
            {
          
          # Calculate Test of Deficit (TD) with singcar package
            res <- singcar::TD(
              case = case, 
              controls = controls,
              sd = sd,
              sample_size = sample_size,
              alternative = alternative,
              conf_int = conf_int,
              conf_level = conf_lev/100,
              conf_int_spec = 0.01,
              na.rm = TRUE)
          
          # Extract Z-CC and Point estimate of the proportion
          z_cc <- as.numeric(res$estimate[1])
          point_estim <- as.numeric(res$estimate[2])
            
          # Extract CI based on UI (if user wants or not CI)
          val_ci <- ifelse(conf_int, paste0(as.numeric(res$interval[1]), "%"), "")
          val_low <- ifelse(conf_int, as.numeric(res$interval[2]), NA_real_)
          val_up <- ifelse(conf_int, as.numeric(res$interval[3]), NA_real_)
            
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
          ci_int <- val_ci
          lower <- val_low
          upper <- val_up
          
          # Table Output
          table <- self$results$resultsTable
          
          # List TD Output values
          table$setRow(
            rowNo = 1, 
            values = list(
              alternative = alt_out,
              df = df_val,
              t_val = t_stat,
              p_val = p_val,
              z_cc = z_cc,
              point_estim = point_estim,
              ci_interval = ci_int,
              ci_lower = lower,
              ci_upper = upper)
            )
          
          # Text Output
          txt_output <- paste(capture.output(res), collapse = "\n")
          self$results$txtOutput$setContent(txt_output)
          }
          },
        
        .plot = function(image, ggtheme, theme, ...) {
          
       # Look for UI input data
          case <- self$options$case
          controls <- self$options$controls
          sd <- self$options$sd
          
       # Check for not-empty UI input data
          if (is.null(case) || is.null(controls) || is.null(sd) ||
              is.na(case) || is.na(controls) || is.na(sd) || sd <= 0)
            return()
          
       # X axis limits definition 
          min_x <- min(controls - 4 * sd, case - 1 * sd)
          max_x <- max(controls + 4 * sd, case + 1 * sd)
          
       # Control points for curve density
          x_vals <- seq(min_x, max_x, length.out = 1000)
          y_vals <- dnorm(x_vals, mean = controls, sd = sd)
          df_plot <- data.frame(x = x_vals, y = y_vals)
          
        # Design plot with ggplot package
          prob_plot <- ggplot2::ggplot(
            df_plot, 
            ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y), fill = "#bfd3e6", alpha = 0.6) + 
            ggplot2::geom_line(color = "#8c96c6", linewidth = 1) + 
            ggplot2::geom_vline(xintercept = case, color = "#d95f02", linewidth = 1.2, linetype = "dashed") + 
            ggplot2::annotate("text", x = case, y = max(df_plot$y) * 1.05, label = " Case", 
                              color = "#d95f02", hjust = -0.1, fontface = "bold", size = 5) +
            ggplot2::labs(
              title = "Test of Deficit: Case vs Controls Distribution",
              x = "Score",
              y = "Control Distribution"
            ) +
            ggtheme +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 16),
              axis.title = ggplot2::element_text(size = 14)
            )
          
          print(prob_plot)
          TRUE
          }
       )
    )