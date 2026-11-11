# =================================================================================================
# Module:     B-DISSOCS - Single Case Comparison
# File:       bdissocs.b.R (Back-end Analysis Script)
# Purpose:    R computational logic for executing single-case deficit tests 
#             using the singcar package (Rittmo & McIntosh, 2021) in jamovi.
#             This script processes UI inputs, runs the statistical functions, 
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
#   - Crawford, J. R. & Garthwaite, P.H. (2005). Testing for suspected impairments and dissociations 
#     in single-case studies in neuropsychology: Evaluation of alternatives using Monte Carlo 
#     simulations and revised tests for dissociations. Neuropsychology, 19, 318-331.
#     https://doi.org/10.1037/0894-4105.19.3.318
#
#   - Rittmo J.Ö., & McIntosh R.D. (2021). singcar: Comparing single cases 
#     to small samples in R. Journal of Open Source Software, 6(68), 3887. 
#     https://doi.org/10.21105/joss.03887
#
# Note:       Developed to provide a zero-code graphical user interface, 
#             making single-case comparisons easily accessible for clinical 
#             and experimental neuropsychology without compromising scientific rigor.
# =================================================================================================

bdissocsClass <- if (
  requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "bdissocsClass",
    inherit = bdissocsBase,
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
        credible_interval <- self$options$credible_interval # Credible Interval toggle
        credible_level    <- self$options$credible_level    # Credible Interval Level (%)
        iter              <- self$options$iter              # Number of Monte Carlo iterations
        unstandardised    <- self$options$unstandardised    # Score scale option
        
        alt_a <- if (!is.null(self$options$alternative_a)) self$options$alternative_a else "less"
        alt_b <- if (!is.null(self$options$alternative_b)) self$options$alternative_b else "less"
        
        # Silent check for empty inputs on initial run (prevents empty interface errors)
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
          jmvcore::reject("The Bayesian Dissociation Test requires a control sample of at least 4.")
          return()
        }
        if (correlation <= -1 || correlation >= 1) {
          jmvcore::reject("The correlation must be strictly between −1 and 1.")
          return()
        }
        
        # Execute Bayesian Standardised Difference Test (BSDT)
        res_bsdt <- singcar::BSDT(
          case_a = case_a, case_b = case_b,
          controls_a = controls_a, controls_b = controls_b,
          sd_a = sd_a, sd_b = sd_b,
          sample_size = sample_size, r_ab = correlation,
          int_level = credible_level / 100, iter = iter,
          alternative = "two.sided", unstandardised = unstandardised,
          calibrated = TRUE, na.rm = TRUE  
        )
        
        # Execute Bayesian Test of Deficit (BTD) for Task A and Task B
        res_btd_a <- singcar::BTD(
          case = case_a, controls = controls_a, sd = sd_a,
          sample_size = sample_size, alternative = alt_a,
          int_level = credible_level / 100, iter = iter, na.rm = TRUE
        )
        
        res_btd_b <- singcar::BTD(
          case = case_b, controls = controls_b, sd = sd_b,
          sample_size = sample_size, alternative = alt_b,
          int_level = credible_level / 100, iter = iter, na.rm = TRUE
        )
        
        # Populate results tables and generate report
        private$.populateTables(res_bsdt, res_btd_a, res_btd_b, seed_val, iter, alt_a, alt_b, credible_interval, credible_level, unstandardised)
        
      },
      
      .populateTables = function(res_bsdt, res_btd_a, res_btd_b, seed_val, iter, alt_a, alt_b, credible_interval, credible_level, unstandardised) {
        # Extract parameters, p-values, and estimates
        df_val      <- as.integer(res_bsdt$parameter)
        pvalue_ab   <- as.numeric(res_bsdt$p.value)
        pvalue_a    <- as.numeric(res_btd_a$p.value)
        pvalue_b    <- as.numeric(res_btd_b$p.value)
        z_cc_a      <- as.numeric(res_bsdt$estimate[1])
        z_cc_b      <- as.numeric(res_bsdt$estimate[2])
        z_dcc       <- as.numeric(res_bsdt$estimate[3])
        point_estim <- as.numeric(res_bsdt$estimate[4])
        
        # Perform clinical classification of dissociation
        class_res <- classifyDissociation(
          p_a = pvalue_a,
          p_b = pvalue_b,
          p_ab = pvalue_ab,
          alpha = 0.05
        )
        
        # Extract Credible Interval values
        val_low <- if (isTRUE(credible_interval)) as.numeric(res_bsdt$interval[2]) else NA_real_
        val_up  <- if (isTRUE(credible_interval)) as.numeric(res_bsdt$interval[3]) else NA_real_
        
        # Populate Dissociation Classification Table
        table_class <- self$results$classTable
        table_class$setRow(
          rowNo = 1, 
          values = list(
            score_a = class_res$score_a,
            score_b = class_res$score_b,
            dissoc_type = class_res$dissoc_type
          )
        )
        
        dir_note_a <- if (alt_a == "less") "case < control" else "case > control"
        dir_note_b <- if (alt_b == "less") "case < control" else "case > control"
        table_class$setNote(
          "criteria_info",
          paste0("The difference test is two-sided. ",
                 "Deficit tests are one-tailed (Task A: ", dir_note_a, 
                 "; Task B: ", dir_note_b, "; alpha = .05). Theoretical criteria: Crawford & Garthwaite (2005).")
        )
        
        # Results Table Setup
        table_stats <- self$results$resultsTable
        
        ci_title <- if (isTRUE(unstandardised)) {
          paste0(credible_level, "% Credible Interval for the Unstandardised Discrepancy")
        } else {
          paste0(credible_level, "% Credible Interval for Z-DCC")
        }
        table_stats$getColumn("ci_lower")$setSuperTitle(ci_title)
        table_stats$getColumn("ci_upper")$setSuperTitle(ci_title)
        
        # Populate Statistical Results Table
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
            point_estim = point_estim,
            ci_lower    = val_low,
            ci_upper    = val_up
          )
        )
        
        # Footnote: Hypothesis, Theoretical Criteria, and Monte Carlo Simulation Details
        note_text <- paste0(
          "Hypothesis: Task A (", dir_note_a, "), Task B (", dir_note_b,"), Difference (case ≠ controls). ",  "\n",
          "Theoretical criteria: Bayesian Standardised Difference Test (BSDT; Crawford & Garthwaite, 2005, 2007). ", "\n",
          "Monte Carlo simulations based on ", format(iter, big.mark = ","), " iterations (random seed: ", seed_val, ")."
        )
        
        if (isTRUE(unstandardised)) {
          note_text <- paste0(note_text, " The credible interval and p (Task A - B) use unstandardised scores; the Z-DCC column and the Task A/Task B tests do not.")
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