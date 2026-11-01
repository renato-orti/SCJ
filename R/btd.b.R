# =================================================================================================
# Module:     BTD - Single Case Comparison with Bayes methodology
# File:       btd.b.R (Back-end Analysis Script)
# Purpose:    R computational logic for executing single-case bayes deficit tests 
#             using the singcar package (Rittmo & McIntosh, 2021) in Jamovi.
#             This script processes UI inputs, runs the statistical functions, 
#             formats the confidence/credible intervals, and populates the results table.
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
#   - Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case to a control or normative
#     sample in neuropsychology: Development of a Bayesian approach. Cognitive Neuropsychology,
#     24(4), 343-372.
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

btdClass <- if (
  requireNamespace('jmvcore', quietly = TRUE)) 
  R6::R6Class(
    "btdClass",
    inherit = btdBase,
    private = list(
      .run = function() {
        
      # Set seed for Monte Carlo simulation reproducibility
        seed_val <- if (!is.null(self$options$seed)) self$options$seed else 42
        set.seed(seed_val)
        
      # Extract UI input data
        case        <- self$options$case               # Single Case Score
        controls    <- self$options$controls           # Control Mean
        sd          <- self$options$sd                 # Control SD
        sample_size <- self$options$sample_size        # Control Sample Size
        alternative <- self$options$alternative        # Alternative Hypothesis
        cred_int    <- self$options$credible_interval  # Credible Interval toggle
        cred_lev    <- self$options$credible_level     # Credible Interval Level (%)
        iter        <- self$options$iter               # Number of Monte Carlo iterations
        
      # Silent check for empty inputs on initial run
        if (is.null(case) || is.null(controls) || is.null(sd) || is.null(sample_size) ||
            is.null(cred_lev) || is.null(iter) ||
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
        if (!is.null(iter) && (iter < 1000 || iter > 10000000)) {
          jmvcore::reject("Number of Monte Carlo iterations must be an integer between 1 and 1,000,000.")
          return()
        }
        
      # Execute Bayesian Test of Deficit (BTD) with singcar package
        res <- singcar::BTD(
          case = case, 
          controls = controls,
          sd = sd,
          sample_size = sample_size,
          alternative = alternative,
          int_level = cred_lev / 100,
          iter = iter,
          na.rm = TRUE
        )
        
      # Extract Z-CC and Point estimate of the proportion
        z_cc        <- as.numeric(res$estimate[1])
        point_estim <- as.numeric(res$estimate[2])
        
      # Extract Credible Interval values
        val_low <- if (isTRUE(cred_int)) as.numeric(res$interval[2]) else NA_real_
        val_up  <- if (isTRUE(cred_int)) as.numeric(res$interval[3]) else NA_real_
        
      # Format alternative hypothesis output string
        alt_out <- switch(
          alternative,
          "two.sided" = "Two Sided (Case ≠ Control)",
          "less"      = "Less (Case < Control)",
          "greater"   = "Greater (Case > Control)",
          alternative
        )
        
        df_val <- as.integer(res$parameter)
        p_val  <- as.numeric(res$p.value)
        
      # Results Table setup
        table_stats <- self$results$resultsTable
        
        ci_title <- paste0(cred_lev, "% Credible Interval")
        table_stats$getColumn("ci_lower")$setSuperTitle(ci_title)
        table_stats$getColumn("ci_upper")$setSuperTitle(ci_title)
        
      # Populate Results Table 
        table_stats$setRow(
          rowNo = 1, 
          values = list(
            z_cc        = z_cc,
            df          = df_val,
            p_val       = p_val,
            point_estim = point_estim,
            ci_lower    = val_low,
            ci_upper    = val_up
          )
        )
        
      # Footnote: Hypothesis, Theoretical Criteria, and Monte Carlo Simulation Details
        table_stats$setNote(
          "info",
          paste0("Hypothesis: ", alt_out, ".", "\n",
                 "Theoretical criteria: Bayesian Test of Deficit (BTD; Crawford & Garthwaite, 2007). ",
                 "Monte Carlo simulations based on ", format(iter, big.mark = ","), 
                 " iterations (random seed: ", seed_val, ")."))
        
      # Text Output Compilation
        txt_output <- paste(capture.output(res), collapse = "\n")
        self$results$textOutput$setContent(txt_output)
      },
      
      .plot = function(image, ggtheme, theme, ...) {
        
        case     <- self$options$case
        controls <- self$options$controls
        sd       <- self$options$sd
        iter     <- self$options$iter
        seed_val <- if (!is.null(self$options$seed)) self$options$seed else 42
        
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
        
      # Caption formatting for Monte Carlo simulation details
        caption_txt <- if (!is.null(iter)) {
          paste0("Monte Carlo iterations: ", format(iter, big.mark = ","), " (seed: ", seed_val, ")")
        } else NULL
        
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
            y = "Control Density",
            caption = caption_txt
          ) +
          ggtheme
        
        print(prob_plot)
        TRUE
      }
    )
  )