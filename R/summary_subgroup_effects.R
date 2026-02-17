#' Create a Summary of Marginal Subgroup Treatment Effects
#'
#' This function orchestrates calls to `estimate_subgroup_effects` to generate
#' a summary table for specific subgroups.
#'
#' @param brms_fit A `brmsfit` object. Fitted model object from `fit_brms_model()` or
#'   `run_brms_analysis()`. Must contain the necessary attributes for extracting treatment
#'   variable and response type information.
#' @param trt_var A character string or `NULL`. Treatment variable name (coded 0/1). If `NULL`,
#'   automatically extracted from model attributes (set by `fit_brms_model()`). This should
#'   match the treatment variable used in model fitting.
#' @param response_type A character string or `NULL`. Outcome type, one of `"binary"`, `"count"`,
#'   `"continuous"`, or `"survival"`. If `NULL`, automatically extracted from model attributes
#'   (set by `fit_brms_model()`). This determines the appropriate scale for summarizing effects.
#' @param subgroup_vars A character vector or `"auto"`. Subgrouping variable names for which
#'   to generate summaries. Defaults to `"auto"`, which automatically detects treatment
#'   interactions from all formula components (`unshrunktermeffect`, `shprogeffect`,
#'   `shpredeffect`). Cannot be `NULL`.
#'
#' @return `subgroup_summary`. S3 object containing:
#'   \describe{
#'     \item{`estimates`}{`tibble` where each row represents a subgroup, with columns
#'       for `Subgroup`, `Median`, `CI_Lower`, and `CI_Upper` from posterior distribution}
#'     \item{`response_type`}{`character(1)` indicating outcome type}
#'     \item{`ci_level`}{`numeric(1)` credible interval level (default 0.95)}
#'     \item{`trt_var`}{`character(1)` treatment variable name}
#'   }
#'
#' @importFrom checkmate assert_class assert_data_frame assert_string assert_subset
#' @importFrom checkmate assert_choice assert check_string check_character test_character
#' @importFrom dplyr bind_rows
#' @export
summary_subgroup_effects <- function(brms_fit,
                                     trt_var = NULL,
                                     response_type = NULL,
                                     subgroup_vars = "auto") {

  # --- 1. Argument Validation ---
  # Validate inputs to ensure compatibility with brms model structure
  checkmate::assert_class(brms_fit, "brmsfit")

  # Extract trt_var from model attributes if not provided
  if (is.null(trt_var)) {
    trt_var <- attr(brms_fit, "trt_var")
    if (is.null(trt_var)) {
      stop("trt_var must be specified or stored in the model attributes (via fit_brms_model()).")
    }
    message("Using trt_var from model attributes: ", trt_var)
  }
  checkmate::assert_string(trt_var, min.chars = 1)
  checkmate::assert_subset(trt_var, names(brms_fit$data))

  # Extract response_type from model attributes if not provided
  if (is.null(response_type)) {
    response_type <- attr(brms_fit, "response_type")
    if (is.null(response_type)) {
      stop("response_type must be specified or stored in the model attributes (via fit_brms_model()).")
    }
    message("Using response_type from model attributes: ", response_type)
  }
  response_type <- checkmate::assert_choice(response_type,
                                            choices = c("binary", "count", "continuous", "survival")
  )

  # Validate subgroup specification: must be either "auto" or explicit character vector
  checkmate::assert(
    checkmate::check_string(subgroup_vars, pattern = "^auto$"),
    checkmate::check_character(subgroup_vars, null.ok = FALSE, min.len = 1, unique = TRUE)
  )

  # Check that variables exist in data if not "auto"
  if (checkmate::test_character(subgroup_vars, min.len = 1) && !identical(subgroup_vars, "auto")) {
    checkmate::assert_subset(subgroup_vars, names(brms_fit$data))
  }

  # --- 2. Calculate subgroup effects ---
  # estimate_subgroup_effects will detect treatment interactions from all formula components:
  # - unshrunktermeffect: all unshrunk terms (if "auto" is used)
  # - shprogeffect: shrunk prognostic effects (if "auto" is used)  
  # - shpredeffect: shrunk predictive effects (if "auto" is used)
  # Or use the explicitly provided subgroup_vars

  message("--- Calculating specific subgroup effects... ---")
  subgroup_effects_list <- estimate_subgroup_effects(
    brms_fit = brms_fit,
    trt_var = trt_var,
    subgroup_vars = subgroup_vars,
    response_type = response_type
  )

  estimates_tbl <- subgroup_effects_list$estimates

  # --- 3. Finalize the output object ---
  summary_output <- list(
    estimates = estimates_tbl,
    response_type = response_type,
    ci_level = 0.95,
    trt_var = trt_var
  )

  class(summary_output) <- "subgroup_summary"

  return(summary_output)
}


#' Combine Multiple Subgroup Summaries for Comparison Plotting
#'
#' Combines estimates from multiple subgroup summary objects into a single
#' summary object that can be plotted to compare models.
#'
#' @param summary_list Named list of `subgroup_summary` objects.
#'
#' @return `subgroup_summary`. Combined object with a "Model" column in estimates.
#'
#' @examples
#' \dontrun{
#' summary1 <- summary_subgroup_effects(brms_fit = model1)
#' summary2 <- summary_subgroup_effects(brms_fit = model2)
#' 
#' combined <- combine_summaries(list(
#'   "One-way model" = summary1,
#'   "Global" = summary2
#' ))
#' plot(combined)
#' }
#'
#' @importFrom dplyr bind_rows mutate
#' @export
combine_summaries <- function(summary_list) {
  checkmate::assert_list(summary_list, min.len = 1)
  
  # Validate all elements are subgroup_summary objects
  for (i in seq_along(summary_list)) {
    checkmate::assert_class(summary_list[[i]], "subgroup_summary", 
                           .var.name = paste0("summary_list[[", i, "]]"))
  }
  
  # Ensure list is named
  if (is.null(names(summary_list))) {
    names(summary_list) <- paste0("Model", seq_along(summary_list))
  }
  
  # Combine all estimates with model labels
  combined_estimates <- bind_rows(lapply(names(summary_list), function(model_name) {
    summary_list[[model_name]]$estimates %>%
      mutate(Model = model_name)
  }))
  
  # Create combined summary object (use first element as template)
  combined_summary <- list(
    estimates = combined_estimates,
    response_type = summary_list[[1]]$response_type,
    ci_level = summary_list[[1]]$ci_level,
    trt_var = summary_list[[1]]$trt_var,
    is_comparison = TRUE  # Flag to indicate this is a comparison plot
  )
  
  class(combined_summary) <- "subgroup_summary"
  return(combined_summary)
}


#' Plot Marginal Subgroup Treatment Effects
#'
#' Creates a forest plot from a `subgroup_summary` object or list of objects.
#'
#' @param x `subgroup_summary` or `list`. Object created by `summary_subgroup_effects()`,
#'   or a named list of such objects for comparison plots.
#' @param x_lab `character(1)` or `NULL`. Custom label for x-axis.
#' @param title `character(1)` or `NULL`. Custom title for plot.
#' @param ... Additional arguments (currently unused).
#'
#' @return `ggplot`. Forest plot visualization of subgroup treatment effects.
#'
#' @details
#' This function creates forest plots for subgroup treatment effects. It supports two modes:
#' 
#' **Single Model Plot**: When `x` is a single `subgroup_summary` object, creates a 
#' traditional forest plot with estimates and confidence intervals displayed as text.
#' 
#' **Multiple Model Comparison**: When `x` is a named list of `subgroup_summary` objects,
#' creates a comparative plot where different models are distinguished by colors and shapes.
#' This is useful for comparing one-way models vs global models, different prior specifications, or
#' sensitivity analyses.
#'
#' @examples
#' \dontrun{
#' # Single model plot
#' summary1 <- summary_subgroup_effects(brms_fit = model1)
#' plot(summary1, title = "Subgroup Effects")
#' 
#' # Multiple model comparison
#' summary2 <- summary_subgroup_effects(brms_fit = model2)
#' comparison <- list(
#'   "Model 1" = summary1,
#'   "Model 2" = summary2
#' )
#' plot(comparison, title = "Model Comparison")
#' }
#'
#' @importFrom dplyr %>% mutate arrange bind_rows
#' @importFrom stringr str_extract
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbar geom_point geom_text
#'   annotate labs theme_classic theme element_text element_blank margin
#'   coord_cartesian scale_shape_manual scale_color_manual position_dodge
#' @importFrom checkmate assert_class assert_data_frame assert_string assert_list
#' @exportS3Method graphics::plot
plot.subgroup_summary <- function(x, x_lab = NULL, title = NULL, ...) {

  # --- 1. Input Validation & Detection of Multiple Models ---
  checkmate::assert_class(x, "subgroup_summary")
  checkmate::assert_string(x_lab, null.ok = TRUE)
  checkmate::assert_string(title, null.ok = TRUE)
  checkmate::assert_data_frame(x$estimates, min.rows = 1)
  
  # Check if this is a comparison plot (has Model column)
  is_comparison <- "Model" %in% names(x$estimates)
  response_type <- x$response_type

  # --- 2. Data Preparation ---
  message("Preparing data for plotting...")

  null_effect_line <- switch(response_type, continuous = 0, binary = 1, count = 1, survival = 1, 0)
  if (is.null(x_lab)) {
    x_lab <- switch(response_type, continuous = "Difference in Mean Outcome", binary = "Odds Ratio", count = "Rate Ratio", survival = "Average Hazard Ratio (AHR)", "Treatment Effect")
  }

  # Prepare plot data
  plot_data <- x$estimates %>%
    mutate(
      estimate_label = sprintf("%.2f (%.2f to %.2f)", .data$Median, .data$CI_Lower, .data$CI_Upper),
      sort_group = stringr::str_extract(.data$Subgroup, "^[^:]+")
    )
  
  # Add Model as factor if this is a comparison plot
  if (is_comparison) {
    plot_data <- plot_data %>%
      arrange(.data$sort_group, .data$Subgroup, .data$Model) %>%  # Sort by subgroup, then model
      mutate(
        Model = factor(.data$Model, levels = unique(.data$Model)),
        y_axis_label = factor(.data$Subgroup, levels = rev(unique(.data$Subgroup)))
      )
  } else {
    plot_data <- plot_data %>%
      arrange(.data$sort_group, .data$Subgroup) %>%
      mutate(y_axis_label = factor(.data$Subgroup, levels = rev(unique(.data$Subgroup))))
  }

  # --- 3. Determine Plot Limits & Table Positions ---
  min_ci <- min(plot_data$CI_Lower, na.rm = TRUE)
  max_ci <- max(plot_data$CI_Upper, na.rm = TRUE)
  x_range <- max_ci - min_ci
  
  # Adjust position for estimate column based on whether we have multiple models
  if (is_comparison) {
    pos_estimate <- max_ci + (x_range * 0.25)  # More space for legend
  } else {
    pos_estimate <- max_ci + (x_range * 0.15)
  }

  # --- 4. Build the ggplot ---
  message("Generating plot...")
  
  if (is_comparison) {
    # Multiple models: use colors and shapes to distinguish
    n_models <- length(unique(plot_data$Model))
    shapes <- c(22, 21, 24, 23, 25)[1:n_models]
    colors <- c("#F8766D", "#00BA38", "#619CFF", "#C77CFF", "#FF9999")[1:n_models]
    
    dodge_width <- 0.5
    
    p <- ggplot(plot_data, aes(y = .data$y_axis_label, x = .data$Median, 
                                color = .data$Model, shape = .data$Model)) +
      geom_vline(xintercept = null_effect_line, linetype = "dashed", color = "grey50") +
      geom_errorbar(aes(xmin = .data$CI_Lower, xmax = .data$CI_Upper), 
                    width = 0.2, 
                    position = position_dodge(width = dodge_width)) +
      geom_point(size = 3, 
                 position = position_dodge(width = dodge_width)) +
      scale_shape_manual(values = shapes) +
      scale_color_manual(values = colors) +
      labs(title = title, x = x_lab, y = "", color = "Model", shape = "Model") +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(hjust = 0, color = "black"),
        legend.position = "bottom",
        plot.margin = margin(1, 2, 1, 1, "lines")
      ) +
      coord_cartesian(clip = "off")
    
  } else {
    # Single model: original plot style
    p <- ggplot(plot_data, aes(y = .data$y_axis_label, x = .data$Median)) +
      geom_vline(xintercept = null_effect_line, linetype = "dashed", color = "grey50") +
      geom_errorbar(aes(xmin = .data$CI_Lower, xmax = .data$CI_Upper), 
                    width = 0.2, color = "black", orientation = "y") +
      geom_point(shape = 22, size = 3, fill = "black", color = "black") +
      geom_text(aes(label = .data$estimate_label), x = pos_estimate, hjust = 0, size = 3.5) +
      annotate("text", x = pos_estimate, y = nrow(plot_data) + 0.8, 
               label = "Estimate (95% CI)", hjust = 0, fontface = "bold", size = 3.5) +
      labs(title = title, x = x_lab, y = "") +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(hjust = 0, color = "black"),
        legend.position = "none",
        plot.margin = margin(1, 8, 1, 1, "lines")
      ) +
      coord_cartesian(clip = "off")
  }

  message("Done.")
  return(p)
}
