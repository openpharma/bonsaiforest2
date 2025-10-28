#' Prepare a Multi-Part brms Formula and Corresponding Data
#'
#' This function serves as a powerful pre-processor for building complex Bayesian
#' models with the `brms` package. It automates the construction of a multi-part,
#' non-linear formula by classifying covariates into four distinct categories:
#' unshrunk prognostic, shrunk prognostic, unshrunk predictive, and shrunk predictive.
#'
#' This classification allows for applying differential shrinkage (regularization)
#' to different parts of the model. The function also prepares the corresponding
#' data, including the automatic creation of dummy variables for predictive effects (interactions).
#'
#' @section Key Features:
#' \itemize{
#'   \item \strong{Multi-Part Formula Construction:} It generates a `brmsformula`
#'     object with up to four distinct linear components (`unprogeffect`,
#'     `shprogeffect`, `unpredeffect`, `shpredeffect`), which are combined
#'     in a non-linear model. This allows for assigning different priors to each
#'     component.
#'   \item \strong{Automated Interaction Handling:} Predictive terms (e.g.,
#'     `~ trt:subgroup`) are automatically expanded. For each level of the
#'     `subgroup` variable, a new dummy column is created in the dataset,
#'     representing the interaction, simplifying the modeling of treatment effect
#'     heterogeneity.
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
#' The treatment variable is converted to an integer (0/1), and new columns for
#' interaction terms are created. The returned `data` object should be used in the
#' subsequent call to `brms::brm()`, not the original data.
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
#' @param response_formula_str A character string for the response part, e.g.,
#'   "outcome ~ trt", for count models "n_events + offset(log(days)) ~ trt" or for
#'   survival models "Surv(time,status) ~ trt".
#' @param shrunk_predictive_formula_str Predictive terms to be shrunk ('shpredeffect').
#'   These are typically interactions with the treatment variable.
#' @param unshrunk_prognostic_formula_str Prognostic terms not to be shrunk
#'   ('unprogeffect'). These are main effects assumed to be important.
#' @param unshrunk_predictive_formula_str Predictive terms not to be shrunk
#'   ('unpredeffect').
#' @param shrunk_prognostic_formula_str Prognostic terms to be shrunk
#'   ('shprogeffect'). These are main effects where regularization is desired.
#' @param response_type The type of outcome variable. One of "binary", "count",
#'   "continuous", or "survival".
#' @param stratification_formula_str A formula string specifying the stratification
#'   variable, e.g., "~ strata_var".
#'
#' @return A list with `formula` (a `brmsformula` object) and `data` (the
#'   modified data.frame).
#'
#' @importFrom stringr str_squish str_split str_starts str_match str_trim str_remove str_detect str_replace_all
#' @importFrom checkmate assert_data_frame assert_string assert_choice assert_subset assert_numeric assert_character
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
#'   # 2. Run the function
#'   prepared_model <- prepare_formula_model(
#'     data = sim_data,
#'     response_formula_str = "Surv(time, status) ~ trt",
#'     shrunk_predictive_formula_str = "~ trt:subgroup",
#'     unshrunk_prognostic_formula_str = "~ age",
#'     shrunk_prognostic_formula_str = "~ region",
#'     response_type = "survival",
#'     stratification_formula_str = "~ region"
#'   )
#'
#'   # 3. View the results
#'   print(prepared_model$formula)
#'   print(head(prepared_model$data))
#' }
#'
prepare_formula_model <- function(data,
                                  response_formula_str,
                                  shrunk_predictive_formula_str = NULL,
                                  unshrunk_prognostic_formula_str = NULL,
                                  unshrunk_predictive_formula_str = NULL,
                                  shrunk_prognostic_formula_str = NULL,
                                  response_type = c("binary", "count", "continuous", "survival"),
                                  stratification_formula_str = NULL) {

  # --- 1. Argument Validation and Initial Parsing ---
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 2)
  checkmate::assert_string(response_formula_str, min.chars = 3, pattern = "~")

  # Check optional formula strings. They must be NULL or a string starting with "~"
  checkmate::assert_string(shrunk_predictive_formula_str, null.ok = TRUE, pattern = "^~")
  checkmate::assert_string(unshrunk_prognostic_formula_str, null.ok = TRUE, pattern = "^~")
  checkmate::assert_string(unshrunk_predictive_formula_str, null.ok = TRUE, pattern = "^~")
  checkmate::assert_string(shrunk_prognostic_formula_str, null.ok = TRUE, pattern = "^~")
  checkmate::assert_string(stratification_formula_str, null.ok = TRUE, pattern = "^~")

  # This replaces match.arg()
  response_type <- checkmate::assert_choice(response_type,
                                            choices = c("binary", "count", "continuous", "survival")
  )

  initial_parts <- .parse_initial_formula(data, response_formula_str)
  processed_data <- initial_parts$data
  response_term <- initial_parts$response_term
  offset_term <- initial_parts$offset_term
  trt_var <- initial_parts$trt_var
  sub_formulas <- list()

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

  # --- 4. Process Predictive Terms  ---
  shrunk_pred_str_resolved <- if (!is.null(term_lists$shrunk_pred) && length(term_lists$shrunk_pred) > 0) {
    # Reconstruct the formula string, starting with ~
    paste("~", paste(term_lists$shrunk_pred, collapse = " + "))
  } else {
    NULL # Pass NULL if the resolved list is empty
  }

  unshrunk_pred_str_resolved <- if (!is.null(term_lists$unshrunk_pred) && length(term_lists$unshrunk_pred) > 0) {
    # Reconstruct the formula string, starting with ~
    paste("~", paste(term_lists$unshrunk_pred, collapse = " + "))
  } else {
    NULL # Pass NULL if the resolved list is empty
  }

  # Now call the processing functions with the resolved formula strings
  shrunk_pred_out <- .process_predictive_terms(
    shrunk_pred_str_resolved, # Use the resolved string
    processed_data,
    trt_var
  )
  # Update processed_data sequentially, as .process_predictive_terms modifies it
  processed_data <- shrunk_pred_out$data

  unshrunk_pred_out <- .process_predictive_terms(
    unshrunk_pred_str_resolved, # Use the resolved string
    processed_data,             # Use the data updated by the previous step
    trt_var
  )
  # Final update to processed_data
  processed_data <- unshrunk_pred_out$data

  # --- 5. Auto-add missing prognostic main effects  ---
  term_lists$unshrunk_prog <- .add_missing_prognostic_effects(
    prognostic_terms = c(term_lists$unshrunk_prog, term_lists$shrunk_prog),
    needed_terms = unique(c(shrunk_pred_out$prognostic_effects, unshrunk_pred_out$prognostic_effects)),
    target_term_list = term_lists$unshrunk_prog
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
    sub_formulas = sub_formulas
  )

  # --- 7. Return Final List ---
  return(list(formula = final_formula_obj, data = processed_data))
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

  formula_parts <- str_squish(str_split(response_formula_str, "~", n = 2)[[1]])
  lhs_str <- formula_parts[1]
  trt_var <- formula_parts[2]

  # --- More Assertions (replacing old stop()) ---
  checkmate::assert_string(trt_var, min.chars = 1)
  checkmate::assert_subset(trt_var, names(data))

  # Split the left-hand side to find the response and any offset
  lhs_terms <- str_squish(str_split(lhs_str, "\\+")[[1]])
  is_offset <- str_starts(lhs_terms, "offset\\(")

  offset_term <- if (any(is_offset)) lhs_terms[is_offset] else NULL
  response_term <- paste(lhs_terms[!is_offset], collapse = " + ")

  if (is.factor(data[[trt_var]])) {
    data[[trt_var]] <- as.integer(data[[trt_var]]) - 1
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

  surv_vars <- str_match(response_part, "Surv\\((.*?),(.*?)\\)")
  time_var <- str_trim(surv_vars[, 2])
  status_var <- str_trim(surv_vars[, 3])
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

  get_terms <- function(s) if (is.null(s) || s == "") NULL else str_squish(str_split(str_remove(s, "~"), "\\+")[[1]])

  unshrunk_prog <- get_terms(unshrunk_prognostic_str)
  shrunk_prog <- get_terms(shrunk_prognostic_str)
  unshrunk_pred <- get_terms(unshrunk_predictive_str)
  shrunk_pred <- get_terms(shrunk_predictive_str)

  # Add treatment to prognostic if missing
  if (!(trt_var %in% unshrunk_prog) && !(trt_var %in% shrunk_prog)) {
    message("Treatment '", trt_var, "' added to unshrunk prognostic terms by default.")
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
    unshrunk_pred = unshrunk_pred, shrunk_pred = shrunk_pred
  ))
}


#' Process Predictive Interaction Terms
#'
#' Creates new dummy variables in the data for interaction terms (e.g.,
#' `trt:subgroup`) and returns the corresponding formula part.
#'
#' This function parses terms like `trt:subgroup` and creates new columns
#' in the data (e.g., `subgroup_S1_x_trt`, `subgroup_S2_x_trt`) for each
#' level of the factor, representing the interaction as a set of dummy variables.
#'
#' @param formula_str A character string formula, e.g., "~ trt:subgroup1", or NULL.
#' @param .data The data.frame to modify.
#' @param .trt_var The character string name of the treatment variable (e.g., "trt").
#'
#' @return A list containing:
#' \item{formula_part}{A string for the new formula component (e.g., "subgroup_S1_x_trt + ...").}
#' \item{data}{The modified data.frame with new dummy interaction columns.}
#' \item{prognostic_effects}{A character vector of main effects (e.g., "subgroup")
#'   that are needed to maintain model hierarchy.}
#' @noRd
#'
.process_predictive_terms <- function(formula_str, .data, .trt_var) {
  # --- Assertions ---
  checkmate::assert_string(formula_str, null.ok = TRUE)
  checkmate::assert_data_frame(.data)
  checkmate::assert_string(.trt_var)

  if (is.null(formula_str)) {
    return(list(formula_part = NULL, data = .data, prognostic_effects = character(0)))
  }
  terms <- attr(terms(as.formula(formula_str)), "term.labels")
  new_formula_terms <- c()
  prognostic_effects_needed <- c()

  for (term in terms) {
    if (str_detect(term, ":")) {
      vars <- str_squish(str_split(term, ":")[[1]])
      subgroup_var <- setdiff(vars, .trt_var)
      if (length(subgroup_var) != 1) {
        warning(paste("Skipping term:", term)); next
      }
      prognostic_effects_needed <- c(prognostic_effects_needed, subgroup_var)
      if (!is.factor(.data[[subgroup_var]])) .data[[subgroup_var]] <- as.factor(.data[[subgroup_var]])

      for (level in levels(.data[[subgroup_var]])) {
        clean_level <- level %>% str_replace_all(c(" " = "_", "<" = "lt", ">" = "gt", "=" = "eq")) %>%
          str_replace_all("[^a-zA-Z0-9_]", "")
        new_col_name <- paste0(subgroup_var, "_", clean_level, "_x_", .trt_var)
        .data[[new_col_name]] <- as.integer(.data[[subgroup_var]] == level & .data[[.trt_var]] == 1)
        new_formula_terms <- c(new_formula_terms, new_col_name)
      }
    }
  }
  final_formula_part <- if (length(new_formula_terms) > 0) paste(new_formula_terms, collapse = " + ") else NULL
  return(list(formula_part = final_formula_part, data = .data, prognostic_effects = unique(prognostic_effects_needed)))
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
      message("Auto-adding missing prognostic effect for interaction: ", needed)
      target_term_list <- c(target_term_list, needed)
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

  strat_var <- str_squish(str_remove(stratification_formula_str, "~"))
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
#'
#' @return A `brms::bf` object, ready for `brms::brm()`.
#' @noRd
#'
.assemble_brms_formula <- function(response_term, offset_term, response_type, unshrunk_prog_terms,
                                   shrunk_prog_terms, unshrunk_pred_formula,
                                   shrunk_pred_formula, sub_formulas = list()) {
  # --- Assertions ---
  checkmate::assert_string(response_term)
  checkmate::assert_string(offset_term, null.ok = TRUE)
  checkmate::assert_string(response_type)
  checkmate::assert_character(unshrunk_prog_terms, null.ok = TRUE)
  checkmate::assert_character(shrunk_prog_terms, null.ok = TRUE)
  checkmate::assert_string(unshrunk_pred_formula, null.ok = TRUE)
  checkmate::assert_string(shrunk_pred_formula, null.ok = TRUE)
  checkmate::assert_list(sub_formulas)

  # Define intercept rules based on response type
  has_intercept <- response_type != "survival"

  placeholders <- c()
  .create_sub_formula <- function(name, terms, intercept = FALSE) {
    if (length(terms) > 0 || intercept) {
      placeholders <<- c(placeholders, name)
      rhs <- paste(terms, collapse = " + ")
      formula_str <- if (nchar(rhs) > 0) paste(name, "~", rhs) else paste(name, "~ 1")
      if (!intercept) formula_str <- paste0(formula_str, " + 0")
      sub_formulas[[name]] <<- brms::lf(formula_str)
    }
  }

  # Create sub-formulas for each component
  .create_sub_formula("unprogeffect", setdiff(unshrunk_prog_terms, "1"), has_intercept)
  .create_sub_formula("shprogeffect", setdiff(shrunk_prog_terms, "1"), FALSE)
  .create_sub_formula("unpredeffect", unshrunk_pred_formula, FALSE)
  .create_sub_formula("shpredeffect", shrunk_pred_formula, FALSE)

  # Create the main formula using placeholders and the full response part
  main_formula_str <- if (length(placeholders) > 0) {
    paste(response_term, "~", paste(placeholders, collapse = " + "))
  } else {
    paste(response_term, "~", if (has_intercept) "1" else "0")
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
