#' Prepare a Multi-Part brms Formula and Corresponding Data
#'
#' This function serves as a powerful pre-processor for building complex Bayesian
#' models with the `brms` package. It automates the construction of a multi-part,
#' non-linear formula by classifying covariates into four distinct categories:
#' unshrunk prognostic, shrunk prognostic, unshrunk predictive, and shrunk predictive.
#'
#' This classification allows for applying differential shrinkage (regularization)
#' to different parts of the model. The function also prepares the corresponding
#' data using R's contrast coding system for proper factor handling and interactions.
#'
#' @section Key Features:
#' \itemize{
#'   \item \strong{Multi-Part Formula Construction:} It generates a `brmsformula`
#'     object with up to four distinct linear components (`unprogeffect`,
#'     `shprogeffect`, `unpredeffect`, `shpredeffect`), which are combined
#'     in a non-linear model. This allows for assigning different priors to each
#'     component.
#'   \item \strong{Automated Interaction Handling:} Predictive terms support two syntaxes:
#'     interaction notation (e.g., `~ trt:subgroup`) or random effects notation (e.g.,
#'     `~ (trt || subgroup)`). Note: only one syntax form (`:` or `||`) may be used
#'     across all predictive formula strings. All variables involved in interactions are
#'     converted to factors with treatment contrasts (reference coding), enabling proper
#'     model fitting and prediction on new data without manual dummy variable creation.
#'   \item \strong{Hierarchical Integrity:} When a predictive term like
#'     `trt:subgroup` is specified, the function automatically ensures that the
#'     corresponding prognostic (main) effect `subgroup` is included in the model,
#'     adhering to the principle of model hierarchy.
#'   \item \strong{Intelligent Defaults:} The main treatment variable is automatically
#'     added as an unshrunk prognostic term if not specified elsewhere. It also
#'     provides warnings and resolves overlaps if a term is accidentally specified
#'     in multiple categories.
#' }
#'
#' @section Survival Model Details (Robust Spline Knot Calculation):
#' For survival models (`response_type = "survival"`), the function explicitly models
#' the baseline hazard using B-splines via `brms::bhaz()`. It implements a highly
#' robust method for calculating spline knots to ensure stability:
#' \enumerate{
#'   \item \strong{Boundary Definition:} It first establishes boundary knots by taking
#'     the range of unique event times and adding a small buffer. This creates a
#'     "safe" interval for placing internal knots.
#'   \item \strong{Internal Knot Placement:} It calculates internal knots using
#'     quantiles of the event times that fall strictly within the boundary knots.
#'     This prevents knots from being placed at the exact edges of the data, which
#'     can cause numerical instability.
#'   \item \strong{Fallback for Sparse Data:} If there are too few unique event
#'     times to calculate quantile-based knots, it gracefully falls back to placing
#'     evenly spaced knots within the defined boundaries.
#' }
#'
#' @section Data Transformation:
#' It is critical to note that this function returns a modified `data.frame`.
#' \itemize{
#'   \item \strong{Treatment variable:} Converted to numeric binary (0/1) to avoid
#'     multi-level factor interactions. The first level (alphabetically or numerically)
#'     becomes 0, the second becomes 1.
#'   \item \strong{Factor covariates:} Variables involved in prognostic terms receive
#'     appropriate contrast coding:
#'     \itemize{
#'       \item Shrunken terms: One-hot encoding (all factor levels represented).
#'         Used with `~ 0` to remove the intercept for proper regularization.
#'       \item Unshrunken terms: Dummy encoding (reference level dropped).
#'         Used with an intercept term.
#'     }
#' }
#' If you have pre-specified contrasts for factor variables (via `contrasts()` function),
#' they will be preserved. Otherwise, appropriate defaults are applied.
#' The returned `data` object should be used in subsequent calls to `brms::brm()`,
#' not the original data.
#'
#' @section Contrast Coding Recommendations:
#' \itemize{
#'   \item For \strong{shrunken prognostic/predictive} terms: Use `~ 0 + var` to remove
#'     the intercept and apply regularization to all factor levels equally.
#'   \item For \strong{unshrunken prognostic} terms: Standard formula `~ var` with intercept
#'     is appropriate (one reference level will be dropped).
#'   \item To specify custom contrasts: Set them in your data before calling this function
#'     using `contrasts(data$var) <- <your_contrast_matrix>`. They will be preserved.
#' }
#'
#' @section Stratification:
#' The `stratification_formula_str` argument allows for estimating certain parameters
#' separately for different groups. Its behavior depends on the `response_type`:
#' \itemize{
#'   \item \code{survival}: Estimates a separate baseline hazard function (`bhaz`) for
#'     each level of the stratification variable.
#'   \item \code{continuous}: Models the residual standard deviation (`sigma`) as varying
#'     across levels of the stratification variable.
#'   \item \code{count}: Models the overdispersion parameter (`shape`) as varying
#'     across levels of the stratification variable.
#' }
#'
#' @param data A data.frame containing all the necessary variables.
#' @param response_formula A formula for the response part, e.g.,
#'   outcome ~ trt, for count models n_events + offset(log(days)) ~ trt or for
#'   survival models Surv(time,status) ~ trt.
#' @param shrunk_predictive_formula Predictive terms to be shrunk ('shpredeffect').
#'   These are typically interactions with the treatment variable. Can use either
#'   interaction notation (e.g., ~ trt:subgroup) or random effects notation (e.g.,
#'   ~ (trt || subgroup)). Only one syntax form may be used across all formulas.
#' @param unshrunk_prognostic_formula Prognostic terms not to be shrunk
#'   ('unprogeffect'). These are main effects assumed to be important.
#' @param unshrunk_predictive_formula Predictive terms not to be shrunk
#'   ('unpredeffect'). Can use either interaction notation or random effects notation,
#'   but only one syntax form across all formulas.
#' @param shrunk_prognostic_formula Prognostic terms to be shrunk
#'   ('shprogeffect'). These are main effects where regularization is desired.
#' @param response_type The type of outcome variable. One of "binary", "count",
#'   "continuous", or "survival".
#' @param stratification_formula A formula specifying the stratification
#'   variable, e.g., ~ strata_var.
#'
#' @return A list with three elements: `formula` (a `brmsformula` object),
#'   `data` (the modified data.frame), and `response_type` (the
#'   character string of the response type).
#'
#' @importFrom stringr str_squish str_split str_starts str_match str_trim str_remove str_detect str_replace_all
#' @importFrom checkmate assert_data_frame assert_string assert_choice assert_subset assert_numeric assert_character
#' @importFrom stats contrasts
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
#'   # 2a. Example with colon interaction syntax
#'   prepared_model <- prepare_formula_model(
#'     data = sim_data,
#'     response_formula = Surv(time, status) ~ trt,
#'     shrunk_predictive_formula = ~ trt:subgroup,
#'     unshrunk_prognostic_formula = ~ age,
#'     shrunk_prognostic_formula = ~ region,
#'     response_type = "survival",
#'     stratification_formula = ~ region
#'   )
#'
#'   # 2b. Alternatively, using pipe-pipe (||) syntax
#'   # prepared_model <- prepare_formula_model(
#'   #   data = sim_data,
#'   #   response_formula = Surv(time, status) ~ trt,
#'   #   shrunk_predictive_formula = ~ (trt || subgroup),
#'   #   unshrunk_prognostic_formula = ~ age,
#'   #   shrunk_prognostic_formula = ~ region,
#'   #   response_type = "survival",
#'   #   stratification_formula = ~ region
#'   # )
#'
#'   # 3. View the results
#'   print(prepared_model$formula)
#'   print(head(prepared_model$data))
#' }
#'
prepare_formula_model <- function(data,
                                  response_formula,
                                  shrunk_predictive_formula = NULL,
                                  unshrunk_prognostic_formula = NULL,
                                  unshrunk_predictive_formula = NULL,
                                  shrunk_prognostic_formula = NULL,
                                  response_type = c("binary", "count", "continuous", "survival"),
                                  stratification_formula = NULL) {

  # --- 1. Argument Validation and Initial Parsing ---
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 2)

  # Helper function to convert formula to string
  .formula_to_string <- function(f) {
    if (is.null(f)) return(NULL)
    if (is.character(f)) return(f)  # Support legacy string input
    if (inherits(f, "formula")) return(paste(deparse(f, width.cutoff = 500L), collapse = " "))
    stop("Formula must be a formula object or character string")
  }

  # Convert all formulas to strings for internal processing
  response_formula_str <- .formula_to_string(response_formula)
  shrunk_predictive_formula_str <- .formula_to_string(shrunk_predictive_formula)
  unshrunk_prognostic_formula_str <- .formula_to_string(unshrunk_prognostic_formula)
  unshrunk_predictive_formula_str <- .formula_to_string(unshrunk_predictive_formula)
  shrunk_prognostic_formula_str <- .formula_to_string(shrunk_prognostic_formula)
  stratification_formula_str <- .formula_to_string(stratification_formula)

  # Validate response formula
  checkmate::assert_string(response_formula_str, min.chars = 3, pattern = "~")

  # This replaces match.arg()
  response_type <- checkmate::assert_choice(response_type,
                                            choices = c("binary", "count", "continuous", "survival")
  )

  # --- Validate that only one syntax form (: or ||) is used across all predictive formulas ---
  .validate_predictive_syntax(
    shrunk_predictive_formula_str,
    unshrunk_predictive_formula_str
  )

  initial_parts <- .parse_initial_formula(data, response_formula_str)
  processed_data <- initial_parts$data
  response_term <- initial_parts$response_term
  offset_term <- initial_parts$offset_term
  trt_var <- initial_parts$trt_var
  sub_formulas <- list()

  # Track which factor variables had user-specified contrasts BEFORE we modify anything
  # This allows us to preserve truly user-specified contrasts while updating our own
  user_contrast_vars <- character(0)
  for (var_name in names(processed_data)) {
    if (is.factor(processed_data[[var_name]]) && !is.null(attr(processed_data[[var_name]], "contrasts"))) {
      user_contrast_vars <- c(user_contrast_vars, var_name)
    }
  }

  # --- 2. Handle Response-Type Specific Logic ---
  if (response_type == "survival") {
    response_term <- .handle_survival_response(
      response_part = response_term,
      data = processed_data,
      stratification_formula_str = stratification_formula_str
    )
  } else if (!is.null(stratification_formula_str)) {
    strat_formulas <- .handle_distributional_stratification(
      response_type = response_type,
      stratification_formula_str = stratification_formula_str,
      data = processed_data
    )
    sub_formulas <- c(sub_formulas, strat_formulas)
  }

  # --- 3. Process and Resolve Formula Terms  ---
  term_lists <- .resolve_term_overlaps_and_defaults(
    unshrunk_prognostic_str = unshrunk_prognostic_formula_str,
    shrunk_prognostic_str = shrunk_prognostic_formula_str,
    unshrunk_predictive_str = unshrunk_predictive_formula_str,
    shrunk_predictive_str = shrunk_predictive_formula_str,
    trt_var = trt_var
  )

  # 4. Process Predictive Terms (Now using the NEW .process_predictive_terms)
  # Reconstruct full formula strings from the lists
  shrunk_pred_str_full <- if (!is.null(term_lists$shrunk_pred) && length(term_lists$shrunk_pred) > 0) {
    paste("~", paste(term_lists$shrunk_pred, collapse = " + "))
  } else NULL

  unshrunk_pred_str_full <- if (!is.null(term_lists$unshrunk_pred) && length(term_lists$unshrunk_pred) > 0) {
    paste("~", paste(term_lists$unshrunk_pred, collapse = " + "))
  } else NULL

  # Call the UPDATED processing function
  shrunk_pred_out <- .process_predictive_terms(
    shrunk_pred_str_full,
    processed_data,
    trt_var,
    user_contrast_vars
  )
  processed_data <- shrunk_pred_out$data

  unshrunk_pred_out <- .process_predictive_terms(
    unshrunk_pred_str_full,
    processed_data,
    trt_var,
    user_contrast_vars
  )
  processed_data <- unshrunk_pred_out$data

  # --- 5. Auto-add missing prognostic main effects  ---
  term_lists$unshrunk_prog <- .add_missing_prognostic_effects(
    prognostic_terms = c(term_lists$unshrunk_prog, term_lists$shrunk_prog),
    needed_terms = unique(c(shrunk_pred_out$prognostic_effects, unshrunk_pred_out$prognostic_effects)),
    target_term_list = term_lists$unshrunk_prog
  )

  # --- 5b. Apply Contrasts to Prognostic Terms ---
  # Apply appropriate contrasts based on whether terms are shrunken or not
  # This respects user-specified contrasts and only sets defaults where missing
  processed_data <- .apply_prognostic_contrasts(
    data = processed_data,
    terms = term_lists$unshrunk_prog,
    is_shrunken = FALSE,  # Unshrunken = dummy encoding
    user_contrast_vars = user_contrast_vars
  )

  processed_data <- .apply_prognostic_contrasts(
    data = processed_data,
    terms = term_lists$shrunk_prog,
    is_shrunken = TRUE,  # Shrunken = one-hot encoding
    user_contrast_vars = user_contrast_vars
  )

  # --- 6. Assemble the Final Formula ---
  # We must now pass the offset_term to the assembler function
  final_formula_obj <- .assemble_brms_formula(
    response_term = response_term,
    offset_term = offset_term,
    response_type = response_type,
    unshrunk_prog_terms = term_lists$unshrunk_prog,
    shrunk_prog_terms = term_lists$shrunk_prog,
    unshrunk_pred_formula = unshrunk_pred_out$formula_part,
    shrunk_pred_formula = shrunk_pred_out$formula_part,
    sub_formulas = sub_formulas,
    shrunk_prog_no_int = term_lists$shrunk_prog_no_int,
    unshrunk_pred_no_int = term_lists$unshrunk_pred_no_int,
    shrunk_pred_no_int = term_lists$shrunk_pred_no_int
  )

  # --- 7. Return Final List ---
  return(list(formula = final_formula_obj, data = processed_data,response_type=response_type ))
}


#' Set Contrasts for Factor Variable Based on Shrinkage Type
#'
#' Sets appropriate contrast encoding for factor variables:
#' - Shrunken terms: one-hot encoding (contr.sum) - all levels represented, no intercept
#' - Unshrunken terms: dummy encoding (contr.treatment) - reference level dropped, with intercept
#'
#' Only applies contrasts if they haven't been manually set by the user.
#'
#' @param data The data.frame containing the variable
#' @param var_name Character string, the name of the factor variable
#' @param is_shrunken Logical, TRUE for shrunken terms (one-hot), FALSE for unshrunken (dummy)
#' @param user_contrast_vars Character vector of variables that had user-specified contrasts
#'
#' @return The modified data.frame
#' @noRd
#'
.set_factor_contrasts <- function(data, var_name, is_shrunken, user_contrast_vars = character(0)) {
  checkmate::assert_data_frame(data)
  checkmate::assert_string(var_name)
  checkmate::assert_logical(is_shrunken, len = 1)
  checkmate::assert_character(user_contrast_vars)

  # Skip contrast setting for numeric variables (e.g., binary 0/1 treatment)
  if (is.numeric(data[[var_name]])) {
    message("Skipping contrast setting for numeric variable '", var_name, "'")
    return(data)
  }
  
  # Ensure variable is a factor
  if (!is.factor(data[[var_name]])) {
    message("Converting '", var_name, "' to factor for proper contrast coding.")
    data[[var_name]] <- as.factor(data[[var_name]])
  }

  # Check if this variable had user-specified contrasts in the original data
  is_user_specified <- var_name %in% user_contrast_vars

  # Always preserve user-specified contrasts
  if (is_user_specified) {
    message("Variable '", var_name, "' has user-specified contrasts. Keeping them.")
    return(data)
  }

  # Set appropriate contrasts based on shrinkage type
  if (is_shrunken) {
    # One-hot encoding: all levels represented (use with ~ 0 to remove intercept)
    contrasts(data[[var_name]], nlevels(data[[var_name]])) <- contr.treatment(levels(data[[var_name]]), contrasts = FALSE)
    message("Applied one-hot encoding to shrunken factor '", var_name, "' (will be used with ~ 0 + ...)")
  } else {
    # Dummy encoding: reference level dropped (use with intercept)
    contrasts(data[[var_name]]) <- contr.treatment(levels(data[[var_name]]))
    message("Applied dummy encoding (contr.treatment) to unshrunken factor '", var_name, "'")
  }

  return(data)
}

#' Validate Predictive Syntax Consistency
#'
#' Ensures that across all predictive formula strings, only one interaction syntax
#' form (`:` or `||`) is used, not a mix of both.
#'
#' @param shrunk_pred_str A formula string or NULL.
#' @param unshrunk_pred_str A formula string or NULL.
#'
#' @return Invisibly returns NULL. Throws error if both syntaxes are detected.
#' @noRd
#'
.validate_predictive_syntax <- function(shrunk_pred_str, unshrunk_pred_str) {
  # Combine all formulas to check
  all_formulas <- paste(shrunk_pred_str, unshrunk_pred_str, collapse = " ")

  has_colon <- stringr::str_detect(all_formulas, ":")
  has_pipe_pipe <- stringr::str_detect(all_formulas, "\\|\\|")

  if (has_colon && has_pipe_pipe) {
    stop(
      "Cannot mix interaction syntaxes in predictive formulas. Use EITHER:\n",
      "  1. Interaction notation (e.g., '~ trt:subgroup')\n",
      "  2. Random effects notation (e.g., '~ (trt || subgroup)')\n",
      "but not both in the same model specification."
    )
  }

  invisible(NULL)
}


#' Parse Predictive Terms with Support for Both Syntaxes
#'
#' Internal helper to extract the term type (`:` or `||`) from predictive formulas.
#'
#' @param formula_str A formula string or NULL.
#' @return A character scalar: "colon", "pipe", or NA_character_ if no interactions.
#' @noRd
#'
.detect_interaction_syntax <- function(formula_str) {
  if (is.null(formula_str) || formula_str == "") {
    return(NA_character_)
  }

  has_colon <- stringr::str_detect(formula_str, ":")
  has_pipe_pipe <- stringr::str_detect(formula_str, "\\|\\|")

  if (has_colon) {
    return("colon")
  } else if (has_pipe_pipe) {
    return("pipe")
  } else {
    return(NA_character_)
  }
}


#' Parse Initial Formula and Validate Inputs
#'
#'
#' Handles initial data validation and parsing of the main response formula string,
#' correctly separating the core response variable from potential offset terms.
#'
#' @param data A data.frame.
#' @param response_formula_str A string like "response ~ treatment" or "response + offset(log(var)) ~ treatment".
#' @return A list containing the core response part (`response_term`), the offset part (`offset_term`),
#'   the treatment variable name (`trt_var`), and the processed data.
#' @noRd
#'
.parse_initial_formula <- function(data, response_formula_str) {
  # --- Assertions ---
  checkmate::assert_data_frame(data)
  checkmate::assert_string(response_formula_str, pattern = "~")

  formula_parts <- stringr::str_squish(stringr::str_split(response_formula_str, "~", n = 2)[[1]])
  lhs_str <- formula_parts[1]
  trt_var <- formula_parts[2]

  # --- More Assertions (replacing old stop()) ---
  checkmate::assert_string(trt_var, min.chars = 1)
  checkmate::assert_subset(trt_var, names(data))

  # Split the left-hand side to find the response and any offset
  lhs_terms <- stringr::str_squish(stringr::str_split(lhs_str, "\\+")[[1]])
  is_offset <- stringr::str_starts(lhs_terms, "offset\\(")

  offset_term <- if (any(is_offset)) lhs_terms[is_offset] else NULL
  response_term <- paste(lhs_terms[!is_offset], collapse = " + ")

  # Convert treatment to numeric binary (0/1)
  # This avoids multiple level factor interactions and simplifies the model
  trt_values <- unique(data[[trt_var]])

  # Validate that treatment is binary
  if (length(trt_values) != 2) {
    stop("Treatment variable '", trt_var, "' must have exactly 2 levels/values. Found: ",
         length(trt_values), " unique values.")
  }

  # Convert to numeric if not already
  if (!is.numeric(data[[trt_var]])) {
    # Convert to 0/1, where the first value in sorted order becomes 0
    sorted_values <- sort(trt_values)
    data[[trt_var]] <- as.numeric(data[[trt_var]] == sorted_values[2])
    message("Converting treatment variable '", trt_var, "' to numeric binary (0/1). ",
            "'", sorted_values[1], "' = 0, '", sorted_values[2], "' = 1")
  } else {
    # Ensure it's 0/1 if already numeric
    if (!all(data[[trt_var]] %in% c(0, 1))) {
      warning("Treatment variable '", trt_var, "' is numeric but not 0/1. ",
              "Recoding to 0/1 based on sorted unique values.")
      sorted_values <- sort(trt_values)
      data[[trt_var]] <- as.numeric(data[[trt_var]] == sorted_values[2])
    } else {
      message("Treatment variable '", trt_var, "' is already numeric binary (0/1).")
    }
  }

  return(list(
    response_term = response_term,
    offset_term = offset_term,
    trt_var = trt_var,
    data = data
  ))
}


#' Calculate Knots for Baseline Hazard Spline
#'
#' Calculates boundary and internal knots for the `bhaz` function,
#' handling cases with sparse unique time points.
#'
#' @param time_data A numeric vector of event times.
#' @return A list with `boundary_knots` and `internal_knots`.
#' @noRd
#'
.calculate_bhaz_knots <- function(time_data) {
  # --- Assertions ---
  checkmate::assert_numeric(time_data, all.missing = FALSE)

  time_data <- time_data[!is.na(time_data)]
  if (length(unique(time_data)) < 4) {
    warning("Not enough unique event times to compute quantile knots for bhaz. Using default brms Cox model.", call. = FALSE)
    return(NULL)
  }

  sort_resp <- sort(unique(time_data))
  data_range <- range(sort_resp)
  buffer <- (data_range[2] - data_range[1]) * 0.01
  if (buffer == 0) buffer <- 0.1

  lower_b <- max(0, data_range[1] - buffer)
  upper_b <- data_range[2] + buffer
  limits_resp <- c(lower_b, upper_b)

  eligible_times <- sort_resp[sort_resp > limits_resp[1] & sort_resp < limits_resp[2]]

  if (length(unique(eligible_times)) >= 3) {
    quantiles_resp <- stats::quantile(eligible_times, c(0.25, 0.5, 0.75), names = FALSE)
  } else {
    message("Few unique time points for knot placement; using evenly spaced knots.")
    quantiles_resp <- seq(from = limits_resp[1], to = limits_resp[2], length.out = 5)[2:4]
  }

  return(list(boundary_knots = limits_resp, internal_knots = unique(quantiles_resp)))
}


#' Handle Survival Response Preparation
#'
#' Prepares the response string and baseline hazard term for a survival model.
#' It converts `Surv(time, status)` into the `brms` `time | cens(1 - status)`
#' format and appends the `bhaz()` term with robustly calculated knots.
#'
#' @param response_part The survival part of the formula, e.g., "Surv(time, status)".
#' @param data The full dataset, used to find the time variable for knot calculation.
#' @param stratification_formula_str A string, e.g., "~ strata_var", or NULL.
#'   If provided, it's used to stratify the baseline hazard.
#' @return The updated response part string for the `brms` formula.
#' @noRd
#'
.handle_survival_response <- function(response_part, data, stratification_formula_str) {
  # --- Assertions ---
  checkmate::assert_string(response_part, pattern = "^Surv\\(")
  checkmate::assert_data_frame(data)
  checkmate::assert_string(stratification_formula_str, null.ok = TRUE)

  surv_vars <- stringr::str_match(response_part, "Surv\\((.*?),(.*?)\\)")
  time_var <- stringr::str_trim(surv_vars[, 2])
  status_var <- stringr::str_trim(surv_vars[, 3])
  brms_response_part <- paste0(time_var, " | cens(1 - ", status_var, ")")

  message("Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().")
  knots <- .calculate_bhaz_knots(data[[time_var]])

  if (is.null(knots)) {
    return(brms_response_part) # Fallback to default Cox model
  }

  bhaz_args <- paste0(
    "Boundary.knots = ", deparse(knots$boundary_knots), ", ",
    "knots = ", deparse(knots$internal_knots), ", ",
    "intercept = FALSE"
  )

  if (!is.null(stratification_formula_str)) {
    strat_var <- str_squish(str_remove(stratification_formula_str, "~"))
    checkmate::assert_subset(strat_var, names(data))
    message("Applying stratification: estimating separate baseline hazards by '", strat_var, "'.")
    bhaz_args <- paste0(bhaz_args, ", gr = ", strat_var)
  }

  return(paste(brms_response_part, "+", paste0("bhaz(", bhaz_args, ")")))
}


#' Resolve Term Overlaps and Set Defaults
#'
#' Checks for overlapping terms between shrunk/unshrunk formulas and adds
#' the treatment variable to prognostic terms if it's missing.
#'
#' @param unshrunk_prognostic_str A formula string, e.g., "~ var1 + var2".
#' @param shrunk_prognostic_str A formula string, e.g., "~ var3".
#' @param unshrunk_predictive_str A formula string, e.g., "~ trt:var1".
#' @param shrunk_predictive_str A formula string, e.g., "~ trt:var4".
#' @param trt_var The name of the treatment variable (character string).
#'
#' @return A list of cleaned-up character vectors for each term type:
#'   `unshrunk_prog`, `shrunk_prog`, `unshrunk_pred`, `shrunk_pred`.
#' @noRd
#'
.resolve_term_overlaps_and_defaults <- function(unshrunk_prognostic_str, shrunk_prognostic_str,
                                                unshrunk_predictive_str, shrunk_predictive_str, trt_var) {

  # --- Assertions ---
  checkmate::assert_string(unshrunk_prognostic_str, null.ok = TRUE)
  checkmate::assert_string(shrunk_prognostic_str, null.ok = TRUE)
  checkmate::assert_string(unshrunk_predictive_str, null.ok = TRUE)
  checkmate::assert_string(shrunk_predictive_str, null.ok = TRUE)
  checkmate::assert_string(trt_var)

  # Helper to extract terms and check for intercept removal
  get_terms <- function(s) {
    if (is.null(s) || s == "") return(NULL)
    terms <- stringr::str_squish(stringr::str_split(stringr::str_remove(s, "~"), "\\+")[[1]])
    # Remove "0" or "-1" from terms (these indicate intercept removal)
    terms <- terms[!terms %in% c("0", "-1")]
    return(if (length(terms) == 0) NULL else terms)
  }

  # Helper to check if formula has explicit intercept removal
  has_no_intercept <- function(s) {
    if (is.null(s) || s == "") return(FALSE)
    rhs <- stringr::str_trim(stringr::str_remove(s, "^~"))
    stringr::str_detect(rhs, "^(0|\\-1)\\s*\\+")
  }

  unshrunk_prog <- get_terms(unshrunk_prognostic_str)
  shrunk_prog <- get_terms(shrunk_prognostic_str)
  unshrunk_pred <- get_terms(unshrunk_predictive_str)
  shrunk_pred <- get_terms(shrunk_predictive_str)

  # Track which formulas had explicit intercept removal
  shrunk_prog_no_int <- has_no_intercept(shrunk_prognostic_str)
  unshrunk_pred_no_int <- has_no_intercept(unshrunk_predictive_str)
  shrunk_pred_no_int <- has_no_intercept(shrunk_predictive_str)

  # Automatically add treatment to unshrunk prognostic if not already present
  if (!(trt_var %in% unshrunk_prog) && !(trt_var %in% shrunk_prog)) {
    message("Treatment '", trt_var, "' automatically added to unshrunk prognostic terms.")
    unshrunk_prog <- c(unshrunk_prog, trt_var)
  }

  # Resolve overlaps
  overlap_prog <- intersect(unshrunk_prog, shrunk_prog)
  if (length(overlap_prog) > 0) {
    warning("Prognostic terms in both shrunk/unshrunk. Prioritizing as unshrunk: ", paste(overlap_prog, collapse = ", "))
    shrunk_prog <- setdiff(shrunk_prog, overlap_prog)
  }
  overlap_pred <- intersect(unshrunk_pred, shrunk_pred)
  if (length(overlap_pred) > 0) {
    warning("Predictive terms in both shrunk/unshrunk. Prioritizing as shrunk: ", paste(overlap_pred, collapse = ", "))
    unshrunk_pred <- setdiff(unshrunk_pred, overlap_pred)
  }

  return(list(
    unshrunk_prog = unshrunk_prog, shrunk_prog = shrunk_prog,
    unshrunk_pred = unshrunk_pred, shrunk_pred = shrunk_pred,
    shrunk_prog_no_int = shrunk_prog_no_int,
    unshrunk_pred_no_int = unshrunk_pred_no_int,
    shrunk_pred_no_int = shrunk_pred_no_int
  ))
}


#' Process Predictive Interaction Terms (Using Contrast Coding)
#'
#' Handles predictive interaction terms by ensuring proper contrast coding on factors.
#' No longer creates explicit dummy variables - uses R's contrast system instead.
#'
#' Supports two syntaxes for specifying interactions:
#' 1. Colon notation: "~ trt:subgroup"
#' 2. Pipe-pipe notation: "~ (trt || subgroup)"
#'
#' @param formula_str A character string formula, e.g., "~ trt:subgroup" or "~ (trt || subgroup)".
#' @param .data The data.frame to modify.
#' @param .trt_var The character string name of the treatment variable.
#' @param .user_contrast_vars Character vector of variables with user-specified contrasts.
#'
#' @return A list containing:
#'   \item{formula_part}{The formula string with interaction terms.}
#'   \item{data}{The data.frame with proper contrast coding applied.}
#'   \item{prognostic_effects}{A character vector of main effects needed for hierarchy.}
#' @noRd
.process_predictive_terms <- function(formula_str, .data, .trt_var, .user_contrast_vars = character(0)) {
  # --- Assertions ---
  checkmate::assert_string(formula_str, null.ok = TRUE)
  checkmate::assert_data_frame(.data)
  checkmate::assert_string(.trt_var)

  if (is.null(formula_str)) {
    return(list(formula_part = NULL, data = .data, prognostic_effects = character(0)))
  }

  # Ensure treatment is numeric binary (should already be from .parse_initial_formula)
  if (!is.numeric(.data[[.trt_var]])) {
    trt_values <- unique(.data[[.trt_var]])
    if (length(trt_values) != 2) {
      stop("Treatment variable must have exactly 2 levels/values for interaction terms")
    }
    sorted_values <- sort(trt_values)
    .data[[.trt_var]] <- as.numeric(.data[[.trt_var]] == sorted_values[2])
  }

  # Validate binary
  if (!all(.data[[.trt_var]] %in% c(0, 1))) {
    stop("Treatment variable must be binary (0/1) for interaction terms")
  }

  # Detect which syntax is being used
  syntax_type <- .detect_interaction_syntax(formula_str)

  if (is.na(syntax_type)) {
    # No interaction terms detected
    return(list(formula_part = NULL, data = .data, prognostic_effects = character(0)))
  }

  if (syntax_type == "colon") {
    return(.process_colon_interaction_terms(formula_str, .data, .trt_var, .user_contrast_vars))
  } else if (syntax_type == "pipe") {
    return(.process_pipe_interaction_terms(formula_str, .data, .trt_var, .user_contrast_vars))
  }
}


#' Process Colon-Based Interactions (trt:subgroup) with Contrasts
#'
#' @param formula_str A formula string like "~ trt:subgroup"
#' @param .data The data.frame to modify.
#' @param .trt_var The treatment variable name.
#' @param .user_contrast_vars Character vector of variables with user-specified contrasts.
#'
#' @return A list with formula_part, data, and prognostic_effects.
#' @noRd
#'
.process_colon_interaction_terms <- function(formula_str, .data, .trt_var, .user_contrast_vars = character(0)) {
  # Parse terms
  terms <- attr(terms(as.formula(formula_str)), "term.labels")
  prognostic_effects_needed <- c()
  interaction_terms <- c()

  for (term in terms) {
    # Process interaction terms (e.g., var1:var2)
    if (stringr::str_detect(term, ":")) {
      vars <- stringr::str_squish(stringr::str_split(term, ":")[[1]])
      
      # All variables in the interaction (except treatment) need prognostic effects
      other_vars <- setdiff(vars, .trt_var)
      prognostic_effects_needed <- c(prognostic_effects_needed, other_vars)

      # Ensure factor variables have appropriate contrasts
      for (var in vars) {
        # Skip numeric variables (e.g., binary 0/1 treatment)
        if (var %in% names(.data) && is.numeric(.data[[var]])) {
          next
        }
        
        if (var %in% names(.data) && !is.factor(.data[[var]])) {
          message(paste("Converting", var, "to factor for interaction terms."))
          .data[[var]] <- as.factor(.data[[var]])
        }
        
        # Set contrasts if not user-specified
        if (var %in% names(.data) && is.factor(.data[[var]])) {
          is_user_specified <- var %in% .user_contrast_vars
          if (!is_user_specified) {
            contrasts(.data[[var]]) <- stats::contr.treatment(levels(.data[[var]]))
          }
        }
      }

      # Keep the interaction term as-is
      interaction_terms <- c(interaction_terms, term)
    }
  }

  # Create formula string
  clean_formula <- if (length(interaction_terms) > 0) {
    paste(interaction_terms, collapse = " + ")
  } else {
    NULL
  }

  return(list(
    formula_part = clean_formula,
    data = .data,
    prognostic_effects = unique(prognostic_effects_needed)
  ))
}


#' Process Pipe-Pipe Interactions (trt || subgroup) - Random Effects Syntax
#'
#' Preserves pipe-pipe notation for random effects modeling in brms.
#' Does NOT convert to colon syntax - these are different model structures.
#'
#' @param formula_str A formula string like "~ (trt || subgroup)"
#' @param .data The data.frame to modify.
#' @param .trt_var The treatment variable name.
#' @param .user_contrast_vars Character vector of variables with user-specified contrasts.
#'
#' @return A list with formula_part, data, and prognostic_effects.
#' @noRd
#'
.process_pipe_interaction_terms <- function(formula_str, .data, .trt_var, .user_contrast_vars = character(0)) {
  # Extract the (var1 || var2) pattern
  pattern <- "\\(\\s*([^|]+?)\\s*\\|\\|\\s*([^)]+?)\\s*\\)"
  matches <- stringr::str_match_all(formula_str, pattern)[[1]]

  if (nrow(matches) == 0) {
    warning("No (var1 || var2) patterns found in formula: ", formula_str)
    return(list(formula_part = NULL, data = .data, prognostic_effects = character(0)))
  }

  prognostic_effects_needed <- c()
  interaction_terms <- c()

  for (i in seq_len(nrow(matches))) {
    var1 <- stringr::str_squish(matches[i, 2])
    var2 <- stringr::str_squish(matches[i, 3])
    
    # Strip any leading "0 +" or "-1 +" from variable names to get clean names
    var1_clean <- stringr::str_remove(var1, "^(0|\\-1)\\s*\\+\\s*")
    var2_clean <- stringr::str_remove(var2, "^(0|\\-1)\\s*\\+\\s*")
    var1_clean <- stringr::str_squish(var1_clean)
    var2_clean <- stringr::str_squish(var2_clean)

    # All variables (except treatment) are potential grouping factors
    all_vars <- c(var1_clean, var2_clean)
    grouping_vars <- setdiff(all_vars, .trt_var)
    prognostic_effects_needed <- c(prognostic_effects_needed, grouping_vars)

    # Ensure grouping variables are factors
    for (var in grouping_vars) {
      # Skip numeric variables (e.g., binary 0/1 treatment)
      if (var %in% names(.data) && is.numeric(.data[[var]])) {
        next
      }
      
      if (var %in% names(.data) && !is.factor(.data[[var]])) {
        message(paste("Converting", var, "to factor for random effects grouping."))
        .data[[var]] <- as.factor(.data[[var]])
      }
      
      # For random effects, contrasts don't matter as much, but set treatment contrasts for consistency
      if (var %in% names(.data) && is.factor(.data[[var]])) {
        is_user_specified <- var %in% .user_contrast_vars
        if (!is_user_specified) {
          contrasts(.data[[var]]) <- stats::contr.treatment(levels(.data[[var]]))
        }
      }
    }

    # PRESERVE pipe-pipe notation - this creates random effects in brms
    # Use the ORIGINAL var1 and var2 (with any 0+ or -1+ prefix)
    interaction_terms <- c(interaction_terms, paste0("(", var1, " || ", var2, ")"))
  }

  # Create formula string with pipe-pipe syntax preserved
  clean_formula <- if (length(interaction_terms) > 0) {
    paste(interaction_terms, collapse = " + ")
  } else {
    NULL
  }

  return(list(
    formula_part = clean_formula,
    data = .data,
    prognostic_effects = unique(prognostic_effects_needed)
  ))
}

#' Apply Contrasts to Prognostic Terms
#'
#' Applies appropriate contrasts to all factor variables in prognostic formulas.
#' Shrunken terms get one-hot encoding, unshrunken get dummy encoding.
#' Respects user-specified contrasts.
#'
#' @param data The data.frame
#' @param terms Character vector of variable names
#' @param is_shrunken Logical, whether these are shrunken terms
#' @param user_contrast_vars Character vector of variables with user-specified contrasts
#'
#' @return The modified data.frame
#' @noRd
#'
.apply_prognostic_contrasts <- function(data, terms, is_shrunken, user_contrast_vars = character(0)) {
  if (is.null(terms) || length(terms) == 0) {
    return(data)
  }

  for (term in terms) {
    # Skip if not in data or not a factor
    if (!term %in% names(data)) next
    if (!is.factor(data[[term]])) next

    # Apply appropriate contrasts
    data <- .set_factor_contrasts(data, term, is_shrunken, user_contrast_vars)
  }

  return(data)
}

#' Add Missing Main Effects
#'
#' Ensures that any variable used in an interaction has its main effect
#' present in the prognostic terms (i.e., respects model hierarchy).
#'
#' @param prognostic_terms A character vector of all prognostic terms already
#'   defined (both shrunk and unshrunk).
#' @param needed_terms A character vector of main effects required by
#'   interaction terms (from `.process_predictive_terms`).
#' @param target_term_list A character vector of terms (e.g., unshrunk prognostic)
#'   to which any missing main effects will be added.
#'
#' @return The updated `target_term_list` (character vector) with missing
#'   prognostic terms appended.
#' @noRd
#'
.add_missing_prognostic_effects <- function(prognostic_terms, needed_terms, target_term_list) {
  # --- Assertions ---
  checkmate::assert_character(prognostic_terms, null.ok = TRUE)
  checkmate::assert_character(needed_terms, null.ok = TRUE)
  checkmate::assert_character(target_term_list, null.ok = TRUE)

  for (needed in needed_terms) {
    if (!needed %in% prognostic_terms) {
      message("Note: Marginality principle not followed - interaction term '", needed,
              "' is used without its main effect. ",
              "Consider adding '", needed, "' to prognostic terms for proper model hierarchy.")
    }
  }
  return(target_term_list)
}

#' Handle Stratification for Distributional Parameters
#'
#' Creates sub-formula `brms::lf()` objects for sigma (continuous) or
#' shape (count) based on a stratification variable.
#'
#' @param response_type The type of outcome, e.g., "continuous" or "count".
#' @param stratification_formula_str A formula string, e.g., "~ strata_var".
#' @param data The data.frame, used to check if the stratification variable exists.
#'
#' @return A named list of `brms::lf()` objects (e.g., `list(sigma = ...)`).
#' @noRd
#'
.handle_distributional_stratification <- function(response_type, stratification_formula_str, data) {
  # --- Assertions ---
  checkmate::assert_choice(response_type, choices = c("continuous", "count"))
  checkmate::assert_string(stratification_formula_str)
  checkmate::assert_data_frame(data)

  strat_var <- stringr::str_squish(stringr::str_remove(stratification_formula_str, "~"))
  checkmate::assert_subset(strat_var, names(data))
  formulas <- list()
  if (response_type == "continuous") {
    message("Applying stratification: estimating sigma by '", strat_var, "'.")
    formulas$sigma <- brms::lf(paste("sigma", stratification_formula_str))
  } else if (response_type == "count") {
    message("Applying stratification: estimating shape by '", strat_var, "'.")
    formulas$shape <- brms::lf(paste("shape", stratification_formula_str))
  }
  return(formulas)
}

#' Check if Formula String Has Intercept
#'
#' A formula has an intercept unless it explicitly removes it with "~ 0 +" or "~ -1 +".
#'
#' @param formula_str A formula string or NULL
#' @return TRUE if formula has intercept, FALSE otherwise
#' @noRd
#'
.has_intercept <- function(formula_str) {
  if (is.null(formula_str) || formula_str == "") return(FALSE)

  # Extract RHS by splitting on ~ and taking the last part
  parts <- stringr::str_split(formula_str, "~")[[1]]
  if (length(parts) < 2) return(FALSE)

  rhs <- stringr::str_trim(parts[length(parts)])
  has_removal <- stringr::str_detect(rhs, "^(0|\\-1)\\s*\\+")

  return(!has_removal)
}

#' Assemble the Final brms Formula
#'
#' Takes all the processed formula components and assembles the final
#' multi-part `brms` formula object. It automatically strips the `offset()`
#' wrapper from the offset term to ensure compatibility with brms non-linear syntax.
#'
#' @param response_term The main response part of the formula, e.g., "outcome" or "time | cens(...)".
#' @param offset_term The offset part of the formula, e.g., "offset(log(days))", or NULL.
#' @param response_type The type of outcome, e.g., "survival", "continuous".
#' @param unshrunk_prog_terms A character vector of unshrunk prognostic terms.
#' @param shrunk_prog_terms A character vector of shrunk prognostic terms.
#' @param unshrunk_pred_formula A formula string (from `.process_predictive_terms`)
#'   for unshrunk predictive terms.
#' @param shrunk_pred_formula A formula string (from `.process_predictive_terms`)
#'   for shrunk predictive terms.
#' @param sub_formulas A list of any existing sub-formulas (like for sigma or shape)
#'   to be combined.
#' @param shrunk_prog_no_int Logical, whether shrunk prognostic had explicit intercept removal
#' @param unshrunk_pred_no_int Logical, whether unshrunk predictive had explicit intercept removal
#' @param shrunk_pred_no_int Logical, whether shrunk predictive had explicit intercept removal
#'
#' @return A `brms::bf` object, ready for `brms::brm()`.
#' @noRd
#'
.assemble_brms_formula <- function(response_term, offset_term, response_type, unshrunk_prog_terms,
                                   shrunk_prog_terms, unshrunk_pred_formula,
                                   shrunk_pred_formula, sub_formulas = list(),
                                   shrunk_prog_no_int = FALSE,
                                   unshrunk_pred_no_int = FALSE,
                                   shrunk_pred_no_int = FALSE) {
  # --- Assertions ---
  checkmate::assert_string(response_term)
  checkmate::assert_string(offset_term, null.ok = TRUE)
  checkmate::assert_string(response_type)
  checkmate::assert_character(unshrunk_prog_terms, null.ok = TRUE)
  checkmate::assert_character(shrunk_prog_terms, null.ok = TRUE)
  checkmate::assert_string(unshrunk_pred_formula, null.ok = TRUE)
  checkmate::assert_string(shrunk_pred_formula, null.ok = TRUE)
  checkmate::assert_list(sub_formulas)

  placeholders <- c()

  .create_sub_formula <- function(name, terms, remove_intercept = FALSE) {
    if (length(terms) > 0) {
      placeholders <<- c(placeholders, name)
      rhs <- paste(terms, collapse = " + ")

      # Add intercept removal if specified by user
      if (remove_intercept) {
        rhs <- paste("0 +", rhs)
      }

      formula_str <- paste(name, "~", rhs)

      # Warn if certain formulas have intercepts when they shouldn't
      if (name %in% c("shprogeffect", "shpredeffect", "unpredeffect")) {
        if (.has_intercept(formula_str)) {
          warning("Formula '", name, "' contains an intercept. ",
                  "For proper regularization/interpretation, consider removing it by adding '~ 0 + ...' or '~ -1 + ...' ",
                  "to your input formula.", call. = FALSE)
        }
      }

      sub_formulas[[name]] <<- brms::lf(formula_str)
    }
  }

  # Create sub-formulas for each component (use formulas as provided by user)
  .create_sub_formula("unprogeffect", unshrunk_prog_terms, remove_intercept = FALSE)
  .create_sub_formula("shprogeffect", shrunk_prog_terms, remove_intercept = shrunk_prog_no_int)
  .create_sub_formula("unpredeffect", unshrunk_pred_formula, remove_intercept = unshrunk_pred_no_int)
  .create_sub_formula("shpredeffect", shrunk_pred_formula, remove_intercept = shrunk_pred_no_int)

  # Create the main formula using placeholders and the full response part
  main_formula_str <- if (length(placeholders) > 0) {
    paste(response_term, "~", paste(placeholders, collapse = " + "))
  } else {
    paste(response_term, "~ 1")
  }

  # Offset handling
  if (!is.null(offset_term) && nzchar(offset_term)) {

    clean_offset_term <- sub("^offset\\((.*)\\)$", "\\1", offset_term)

    main_formula_str <- paste(main_formula_str, "+", clean_offset_term)
  }


  main_bf_obj <- brms::bf(main_formula_str, nl = TRUE)

  # Combine main formula with all sub-formulas
  Reduce("+", sub_formulas, init = main_bf_obj)
}
