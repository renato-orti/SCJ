# =================================================================================================
# Module:     BSDT - Single Case Comparison
# File:       bsdt.b.R (Back-end Analysis Script)
# Purpose:    R computational logic for executing single-case deficit tests 
#             using the singcar package (Rittmo & McIntosh, 2021) in jamovi.
#             This script processes UI inputs, runs the statistical functions 
#             and populates the results table.
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
#   - Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case 
#     to a control or normative sample in neuropsychology: Development of a Bayesian approach. 
#     Cognitive Neuropsychology, 24(4), 343-372.
#     https://doi.org/10.1080/02643290701290146
#
#   - Rittmo J.Ö., & McIntosh R.D. (2021). singcar: Comparing single cases 
#     to small samples in R. Journal of Open Source Software, 6(68), 3887. 
#     https://doi.org/10.21105/joss.03887
#
# Note:       Developed to provide a zero-code graphical user interface, 
#             making single-case comparisons easily accessible for clinical 
#             and experimental neuropsychology without compromising scientific rigor.
# =================================================================================================

bsdtClass <- if (
  requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "bsdtClass",
    inherit = bsdtBase,
    private = list(
      .run = function() {
        
        # Set seed for Monte Carlo simulation reproducibility
        seed_val <- if (!is.null(self$options$seed)) self$options$seed else 42
        set.seed(seed_val)
        
        # Extract UI input data
        case_a            <- self$options$case_a            # Single Case Score Task A
        case_b            <- self$options$case_b            # Single Case Score Task B
        controls_a        <- self$options$controls_a        # Control Mean Task A
        controls_b        <- self$options$controls_b        # Control Mean Task B
        sd_a              <- self$options$sd_a              # Control SD Task A
        sd_b              <- self$options$sd_b              # Control SD Task B
        sample_size       <- self$options$sample_size       # Control Sample Size
        correlation       <- self$options$correlation       # Task A-B Correlation in Control Sample
        alternative       <- self$options$alternative       # Alternative Hypothesis
        credible_interval <- self$options$credible_interval # Credible Interval toggle
        credible_level    <- self$options$credible_level    # Credible Interval Level (%)
        iter              <- self$options$iter              # Number of Monte Carlo iterations
        unstandardised    <- self$options$unstandardised    # Score scale option
        
        # Silent check for empty inputs on initial run
        if (is.null(case_a) || is.null(case_b) ||
            is.null(controls_a) || is.null(controls_b) ||
            is.null(sd_a) || is.null(sd_b) ||
            is.null(sample_size) || is.null(correlation) ||
            is.null(credible_level) || is.null(iter) ||
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
        if (sample_size <= 3) {
          jmvcore::reject("The Bayesian Standardised Difference Test requires a control sample of at least 4.")
          return()
        }
        if (correlation <= -1 || correlation >= 1) {
          jmvcore::reject("The correlation must be strictly between −1 and 1.")
          return()
        }
        
        # Execute Bayesian Standardised Difference Test (BSDT)
        res <- singcar::BSDT(
          case_a = case_a,
          case_b = case_b,
          controls_a = controls_a,
          controls_b = controls_b,
          sd_a = sd_a,
          sd_b = sd_b,
          sample_size = sample_size,
          r_ab = correlation,
          int_level = credible_level / 100,
          iter = iter,
          alternative = alternative,
          unstandardised = unstandardised,
          calibrated = TRUE,
          na.rm = TRUE  
        )
        
        # Extract Credible Interval values
        val_low <- if (isTRUE(credible_interval)) as.numeric(res$interval[2]) else NA_real_
        val_up  <- if (isTRUE(credible_interval)) as.numeric(res$interval[3]) else NA_real_
        
        # Extract Z-CC, Z-DCC and Point estimate of the proportion
        z_cc_a      <- as.numeric(res$estimate[1])
        z_cc_b      <- as.numeric(res$estimate[2])
        z_dcc       <- as.numeric(res$estimate[3])
        point_estim <- as.numeric(res$estimate[4])
        
        # Format alternative hypothesis output string
        alt_out <- switch(
          alternative,
          "two.sided" = "Two-sided (A − B in case ≠ controls)",
          "less"      = "Task A lower than Task B, relative to controls",
          "greater"   = "Task A higher than Task B, relative to controls",
          alternative
        )
        
        df_val <- as.integer(res$parameter)
        p_val  <- as.numeric(res$p.value)
        
        # Results Table Setup
        table_stats <- self$results$resultsTable
        
        ci_title <- if (isTRUE(unstandardised)) {
          paste0(credible_level, "% Credible Interval for the Unstandardised Discrepancy")
        } else {
          paste0(credible_level, "% Credible Interval for Z-DCC")
        }
        table_stats$getColumn("ci_lower")$setSuperTitle(ci_title)
        table_stats$getColumn("ci_upper")$setSuperTitle(ci_title)
        
        # Populate Results Table
        table_stats$setRow(
          rowNo = 1, 
          values = list(
            z_cc_a      = z_cc_a,
            z_cc_b      = z_cc_b,
            z_dcc       = z_dcc,
            df          = df_val,
            p_val       = p_val,
            point_estim = point_estim,
            ci_lower    = val_low,
            ci_upper    = val_up
          )
        )
        
        # Footnote: Hypothesis, Theoretical Criteria, and Monte Carlo Simulation Details
        note_text <- paste0(
          "Hypothesis: ", alt_out, ". ", "\n",
          "Theoretical criteria: Bayesian Standardised Difference Test (BSDT; Crawford & Garthwaite, 2007). ", "\n",
          "Monte Carlo simulations based on ", format(iter, big.mark = ","), " iterations (random seed: ", seed_val, ")."
        )
        
        if (isTRUE(unstandardised)) {
          note_text <- paste0(note_text, " Note: The credible interval and p-value use unstandardised scores, while the Z-DCC column does not.")
        }
        
        table_stats$setNote("info", note_text)
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