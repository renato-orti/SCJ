# =================================================================================================
# Module:    TD - Single Case Comparison
# File:      td.b.R (Back-end Analysis Script)
# Purpose:   R computational logic for executing single-case deficit tests 
#            using the singcar package (Rittmo & McIntosh, 2021) in Jamovi.
#            This script processes UI inputs, runs the statistical functions, 
#            formats the confidence intervals, and populates the results table.
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
  requireNamespace('jmvcore', quietly = TRUE)) 
  R6::R6Class(
    "tdClass",
    inherit = tdBase,
    private = list(
      .run = function() {
        
      # Extract UI input data
        case        <- self$options$case                     # Single Case Score
        controls    <- self$options$controls                 # Control Mean
        sd          <- self$options$sd                       # Control SD
        sample_size <- self$options$sample_size              # Control Sample Size
        alternative <- self$options$alternative              # Alternative Hypothesis
        conf_int    <- self$options$confidence_interval      # Confidence Interval toggle
        conf_lev    <- self$options$confidence_level         # Confidence Level (%)
        
      # Silent check for empty inputs on initial run
        if (is.null(case) || is.null(controls) || is.null(sd) || is.null(sample_size) ||
            is.null(conf_lev) ||
            is.na(case) || is.na(controls) || is.na(sd) || is.na(sample_size)
        ) {
          return()
        }
        
      # Input Range Validation with user-friendly rejection messages
        if (sd <= 0) {
          jmvcore::reject("Standard deviation must be strictly positive (> 0).")
          return()
        }
        if (sample_size <= 3) {
          jmvcore::reject("Sample size must be an integer greater than 3.")
          return()
        }
        
      # Calculate Test of Deficit (TD) with singcar package
        res <- singcar::TD(
          case = case, 
          controls = controls,
          sd = sd,
          sample_size = sample_size,
          alternative = alternative,
          conf_int = conf_int,
          conf_level = conf_lev / 100,
          conf_int_spec = 0.01,
          na.rm = TRUE
        )
        
      # Extract Z-CC and Point estimate of the proportion
        z_cc        <- as.numeric(res$estimate[1])
        point_estim <- as.numeric(res$estimate[2])
        
      # Extract CI based on UI settings
        val_low <- if (isTRUE(conf_int) && !is.null(res$interval)) as.numeric(res$interval[2]) else NA_real_
        val_up  <- if (isTRUE(conf_int) && !is.null(res$interval)) as.numeric(res$interval[3]) else NA_real_
        
      # Format alternative hypothesis output string
        alt_out <- switch(
          alternative,
          "two.sided" = "Two Sided (Case ≠ Control)",
          "less"      = "Less (Case < Control)",
          "greater"   = "Greater (Case > Control)",
          alternative
        )
        
        t_stat <- as.numeric(res$statistic)
        df_val <- as.integer(res$parameter)
        p_val  <- as.numeric(res$p.value)
        
      # Results Table Setup
        table_stats <- self$results$resultsTable
        
        ci_title <- paste0(conf_lev, "% Confidence Interval for Z-CC")
        table_stats$getColumn("ci_lower")$setSuperTitle(ci_title)
        table_stats$getColumn("ci_upper")$setSuperTitle(ci_title)
        
      # Populate Results Table
        table_stats$setRow(
          rowNo = 1, 
          values = list(
            t_val       = t_stat,
            df          = df_val,
            p_val       = p_val,
            z_cc        = z_cc,
            point_estim = point_estim,
            ci_lower    = val_low,
            ci_upper    = val_up
          )
        )
        
      # Footnote: Unified Hypothesis and Theoretical Criteria Details
        table_stats$setNote(
          "info",
          paste0("Hypothesis: ", alt_out, ".", "\n",
                 "Theoretical criteria: Modified t-test (TD; Crawford & Howell, 1998).")
        )
        
      # Text Output Compilation
        txt_output <- paste(capture.output(res), collapse = "\n")
        self$results$txtOutput$setContent(txt_output)
      },
      
      .plot = function(image, ggtheme, theme, ...) {
        
        case     <- self$options$case
        controls <- self$options$controls
        sd       <- self$options$sd
        
        if (is.null(case) || is.null(controls) || is.null(sd) ||
            is.na(case) || is.na(controls) || is.na(sd) || sd <= 0) {
          return()
        }
        
        # Definition of X-axis limits
        min_x <- min(controls - 4 * sd, case - 1 * sd)
        max_x <- max(controls + 4 * sd, case + 1 * sd)
        
        # Density curve generation
        x_vals  <- seq(min_x, max_x, length.out = 1000)
        y_vals  <- dnorm(x_vals, mean = controls, sd = sd)
        df_plot <- data.frame(x = x_vals, y = y_vals)
        
        prob_plot <- ggplot2::ggplot(
          df_plot, 
          ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y), fill = "#bfd3e6", alpha = 0.6) + 
          ggplot2::geom_line(color = "#8c96c6", linewidth = 1) + 
          ggplot2::geom_vline(xintercept = case, color = "#d95f02", linewidth = 1.2, linetype = "dashed") + 
          ggplot2::annotate("text", x = case, y = max(df_plot$y) * 1.05, label = " Case", 
                            color = "#d95f02", hjust = -0.1, fontface = "bold") +
          ggplot2::labs(
            x = "Score",
            y = "Control Density"
          ) +
          ggtheme
        
        print(prob_plot)
        TRUE
      }
    )
  )