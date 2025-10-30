#' Create a Summary of Marginal Subgroup Treatment Effects
#'
#' This function orchestrates calls to `estimate_subgroup_effects` to generate
#' a complete summary table. It calculates both the overall marginal effect
#' and the effects for specific subgroups.
#'
#' @param brms_fit A fitted `brmsfit` object.
#' @param original_data The original `data.frame` used for model fitting (before
#'   processing by `prepare_formula_model`).
#' @param trt_var The name of the treatment variable (coded 0/1).
#' @param response_type The type of outcome variable. One of "binary", "count",
#'   "continuous", or "survival".
#' @param subgroup_vars A character vector of subgrouping variable names.
#'   Defaults to "auto", which detects subgroups from the model. Set to `NULL`
#'   to get only the overall effect.
#'
#' @return An object of class `subgroup_summary`, which is a list containing:
#' \item{estimates}{A `tibble` combining the overall and subgroup-specific
#'   effect estimates.}
#' \item{response_type}{The specified response type.}
#' \item{ci_level}{The credible interval level (hard-coded to 0.95).}
#' \item{trt_var}{The name of the treatment variable.}
#'
#' @importFrom checkmate assert_class assert_data_frame assert_string assert_subset
#' @importFrom checkmate assert_choice assert check_string check_character test_character
#'
#' @export
#'
#' @examples
#' if (require("brms") && require("survival")) {
#'   # 1. Create Sample Data
#'   set.seed(123)
#'   n <- 100
#'   sim_data <- data.frame(
#'     time = round(runif(n, 1, 100)),
#'     status = sample(0:1, n, replace = TRUE),
#'     trt = sample(0:1, n, replace = TRUE),
#'     age = rnorm(n, 50, 10),
#'     region = sample(c("A", "B"), n, replace = TRUE),
#'     subgroup = sample(c("S1", "S2", "S3"), n, replace = TRUE)
#'   )
#'   sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
#'   sim_data$region <- as.factor(sim_data$region)
#'   sim_data$subgroup <- as.factor(sim_data$subgroup)
#'
#'   # 2. Run the full analysis
#'   \donttest{
#'   full_fit <- run_brms_analysis(
#'     data = sim_data,
#'     response_formula_str = "Surv(time, status) ~ trt",
#'     response_type = "survival",
#'     shrunk_predictive_formula_str = "~ trt:subgroup",
#'     unshrunk_prognostic_formula_str = "~ age",
#'     shrunk_prognostic_formula_str = "~ region",
#'     chains = 1, iter = 50, warmup = 10, refresh = 0
#'   )
#'
#'   # 3. Get the summary
#'   # This function calls estimate_subgroup_effects internally
#'   eff_summary <- summary_subgroup_effects(
#'     brms_fit = full_fit,
#'     original_data = sim_data,
#'     trt_var = "trt",
#'     response_type = "survival",
#'     subgroup_vars = c("subgroup", "region")
#'   )
#'
#'   print(eff_summary)
#'   }
#' }
summary_subgroup_effects <- function(brms_fit,
                                     original_data,
                                     trt_var,
                                     response_type = c("binary", "count", "continuous", "survival"),
                                     subgroup_vars = "auto") {

  # --- 1. Argument Validation ---
  checkmate::assert_class(brms_fit, "brmsfit")
  checkmate::assert_data_frame(original_data, min.rows = 1)
  checkmate::assert_string(trt_var, min.chars = 1)
  checkmate::assert_subset(trt_var, names(original_data))

  # This replaces match.arg()
  response_type <- checkmate::assert_choice(response_type,
                                            choices = c("binary", "count", "continuous", "survival")
  )

  # Complex check: subgroup_vars must be NULL, the string "auto", or a character vector
  checkmate::assert(
    checkmate::check_string(subgroup_vars, pattern = "^auto$"),
    checkmate::check_character(subgroup_vars, null.ok = TRUE, min.len = 1, unique = TRUE)
  )

  # If it is a character vector (and not "auto"), check that all vars exist
  if (checkmate::test_character(subgroup_vars, null.ok = FALSE, min.len = 1) && !identical(subgroup_vars, "auto")) {
    checkmate::assert_subset(subgroup_vars, names(original_data))
  }

  # --- 2. Always calculate the overall effect ---
  message("--- Calculating overall marginal effect... ---")

  overall_effect_list <- estimate_subgroup_effects(
    brms_fit = brms_fit,
    original_data = original_data,
    trt_var = trt_var,
    subgroup_vars = NULL, # Gets the overall effect
    response_type = response_type
  )

  # Calculate metrics for the overall effect draws
  #overall_metrics <- .calculate_draw_metrics(overall_effect_list$draws)

  # Join metrics. Assumes the estimates table has a column named "Subgroup"
  # that contains the group name (e.g., "Overall").
  overall_effect_tbl <- overall_effect_list$estimates

  subgroup_effects_tbl <- NULL

  # --- 3. Calculate subgroup effects if requested ---
  if (!is.null(subgroup_vars)) {
    message("\n--- Calculating specific subgroup effects... ---")
    subgroup_effects_list <- estimate_subgroup_effects(
      brms_fit = brms_fit,
      original_data = original_data,
      trt_var = trt_var,
      subgroup_vars = subgroup_vars,
      response_type = response_type
    )

    # Calculate metrics for the subgroup effect draws
    #subgroup_metrics <- .calculate_draw_metrics(subgroup_effects_list$draws)

    # Join metrics. Assumes the subgroup estimates table also has a "Subgroup" column.
    subgroup_effects_tbl <- subgroup_effects_list$estimates
  }

  # --- 4. Combine the augmented overall and subgroup estimates ---
  estimates_tbl <- dplyr::bind_rows(overall_effect_tbl, subgroup_effects_tbl)

  # --- 5. Finalize the output object ---
  summary_output <- list(
    estimates = estimates_tbl,
    response_type = response_type,
    ci_level = 0.95,
    trt_var = trt_var
  )

  class(summary_output) <- "subgroup_summary"

  return(summary_output)
}


#' Plot Marginal Subgroup Treatment Effects
#'
#' Creates a forest plot from a `subgroup_summary` object. This version uses
#' simplified labels of the format "variable: level".
#'
#' @param x An object of class `subgroup_summary`.
#' @param x_lab A character string for the x-axis label. If `NULL`, a default
#'   is chosen based on `response_type`.
#' @param title A character string for the plot title.
#' @param ... Additional arguments (not currently used).
#'
#' @return A `ggplot` object representing the forest plot.
#'
#' @importFrom dplyr %>% mutate if_else arrange
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbarh geom_point geom_text
#'   annotate labs theme_classic theme element_text element_blank margin
#'   coord_cartesian
#' @importFrom checkmate assert_class assert_data_frame assert_string
#'
#' @exportS3Method graphics::plot
#'
#' @examples
#' if (require("brms") && require("survival")) {
#'   # 1. Create Sample Data (as in previous examples)
#'   set.seed(123)
#'   n <- 100
#'   sim_data <- data.frame(
#'     time = round(runif(n, 1, 100)),
#'     status = sample(0:1, n, replace = TRUE),
#'     trt = sample(0:1, n, replace = TRUE),
#'     age = rnorm(n, 50, 10),
#'     region = sample(c("A", "B"), n, replace = TRUE),
#'     subgroup = sample(c("S1", "S2", "S3"), n, replace = TRUE)
#'   )
#'   sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
#'   sim_data$region <- as.factor(sim_data$region)
#'   sim_data$subgroup <- as.factor(sim_data$subgroup)
#'
#'   # 2. Run the full analysis (fast example)
#'   \donttest{
#'   full_fit <- run_brms_analysis(
#'     data = sim_data,
#'     response_formula_str = "Surv(time, status) ~ trt",
#'     response_type = "survival",
#'     shrunk_predictive_formula_str = "~ trt:subgroup",
#'     chains = 1, iter = 50, warmup = 10, refresh = 0
#'   )
#'
#'   # 3. Get the summary object
#'   eff_summary <- summary_subgroup_effects(
#'     brms_fit = full_fit,
#'     original_data = sim_data,
#'     trt_var = "trt",
#'     response_type = "survival",
#'     subgroup_vars = c("subgroup", "region")
#'   )
#'
#'   # 4. Plot the object
#'   plot(eff_summary, title = "My Subgroup Analysis")
#'   }
#' }
plot.subgroup_summary <- function(x, x_lab = NULL, title = NULL, ...) {

  # --- 1. Input Validation (with checkmate) ---
  checkmate::assert_class(x, "subgroup_summary")
  checkmate::assert_string(x_lab, null.ok = TRUE)
  checkmate::assert_string(title, null.ok = TRUE)

  # Check that the estimates component is a data.frame with at least one row
  checkmate::assert_data_frame(x$estimates, min.rows = 1)

  # --- 2. Data Preparation (SIMPLIFIED LOGIC) ---
  message("Preparing data for plotting...")

  # Determine null effect and default x-axis label
  null_effect_line <- switch(x$response_type, continuous = 0, binary = 1, count = 1, survival = 1, 0)
  if (is.null(x_lab)) {
    x_lab <- switch(x$response_type, continuous = "Difference in Mean Outcome", binary = "Odds Ratio", count = "Rate Ratio", survival = "Average Hazard Ratio (AHR)", "Treatment Effect")
  }

  # Prepare plot data: format labels and create a clean ordering for the y-axis
  plot_data <- x$estimates %>%
    mutate(
      estimate_label = sprintf("%.2f (%.2f to %.2f)", .data$Median, .data$CI_Lower, .data$CI_Upper),
      sort_group = if_else(.data$Subgroup == "Overall", " Overall", str_extract(.data$Subgroup, "^[^:]+"))
    ) %>%
    arrange(.data$sort_group, .data$Subgroup) %>%
    # Create the y-axis factor. The order is now correct. We reverse it so "Overall" is at the top of the plot.
    mutate(y_axis_label = factor(.data$Subgroup, levels = rev(unique(.$Subgroup))))

  # --- 3. Determine Plot Limits & Table Positions ---
  min_ci <- min(plot_data$CI_Lower, na.rm = TRUE)
  max_ci <- max(plot_data$CI_Upper, na.rm = TRUE)
  x_range <- max_ci - min_ci
  pos_estimate <- max_ci + (x_range * 0.15) # Position for the text table

  # --- 4. Build the ggplot ---
  message("Generating plot...")
  p <- ggplot(plot_data, aes(y = .data$y_axis_label, x = .data$Median)) +

    # --- Additions & Fixes ---

    # 1. (ADDED) Add the vertical null effect line
    geom_vline(xintercept = null_effect_line, linetype = "dashed", color = "grey50") +

    # 2. (EXISTING) The error bars
    geom_errorbar(aes(xmin = .data$CI_Lower, xmax = .data$CI_Upper), height = 0.2, color = "black", orientation = "y") +

    # 3. (ADDED) The point estimate (the missing dot)
    geom_point(shape = 22, size = 3, fill = "black", color = "black") +

    # --- Text Annotations ---

    # Text for the estimate values
    geom_text(aes(label = .data$estimate_label), x = pos_estimate, hjust = 0, size = 3.5)  +

    # Text for the column header
    annotate("text", x = pos_estimate, y = nrow(plot_data) + 0.8, label = "Estimate (95% CI)", hjust = 0, fontface = "bold", size = 3.5) +

    # --- Labels & Theme ---
    labs(title = title, x = x_lab, y = "") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(hjust = 0, color = "black"),
      legend.position = "none",
      # Ensure right margin is large enough for the text table
      plot.margin = margin(1, 8, 1, 1, "lines")
    ) +

    # 4. (ADDED) This is the key fix to prevent text from being cut off
    coord_cartesian(clip = "off")

  message("Done.")
  return(p)
}

