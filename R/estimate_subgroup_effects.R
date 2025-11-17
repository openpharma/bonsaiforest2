#' Estimate Marginal Subgroup Treatment Effects
#'
#' This post-processing function estimates marginal treatment effects for
#' specified subgroups from a fitted `brmsfit` object.
#'
#' @description
#' The function uses a counterfactual, marginal approach based on the posterior
#' predictive distribution. It averages over all other covariates to provide
#' robust estimates of subgroup-specific effects.
#'
#' The estimation follows these steps:
#' \enumerate{
#'   \item Creates two counterfactual datasets based on the model data: one where
#'     all subjects receive "treatment" and one where all receive "control".
#'   \item Generates posterior predictions (e.g., `brms::posterior_epred`) for
#'     both scenarios.
#'   \item Calculates the treatment effect (e.g., difference or ratio) for each
#'     subject and posterior draw.
#'   \item Averages these individual-level effects within each subgroup to
#'     produce the final marginal subgroup estimates.
#' }
#'
#' @param brms_fit A fitted `brmsfit` object from `fit_brms_model()` or
#'   `run_brms_analysis()`.
#' @param original_data A `data.frame` containing the original data before
#'   it was processed by `prepare_formula_model()`. This is essential for
#'   correctly identifying the original subgroup factor levels.
#' @param trt_var A character string specifying the name of the treatment variable.
#' @param subgroup_vars A character vector of subgroup variable names found in
#'   `original_data`. If set to `"auto"` (the default), the function
#'   attempts to automatically identify subgroup variables from the model's
#'   interaction terms.
#' @param response_type The type of outcome variable. One of "binary", "count",
#'   "continuous", or "survival". This determines the scale of the effect
#'   (e.g., risk difference, rate ratio).
#' @param ndraws An integer specifying the number of posterior draws to use.
#'   If `NULL` (default), all available draws are used.
#'
#' @return
#' A `data.frame` (tibble) where each row corresponds to a subgroup
#' (or the "Overall" effect), providing the estimated marginal effect and
#' posterior summaries (e.g., `mean`, `sd`, `q2.5`, `q97.5`).
#'
#' @importFrom checkmate assert_class assert_data_frame assert_string assert_choice
#' @importFrom checkmate assert_character assert_count assert_logical assert_subset
#' @importFrom checkmate assert_list assert_matrix assert_integerish
#' @importFrom checkmate assert check_string
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
#'   # 3. Estimate subgroup effects
#'   # Note: We pass the original `sim_data`, not the processed model data.
#'   effects <- estimate_subgroup_effects(
#'     brms_fit = full_fit,
#'     original_data = sim_data,
#'     trt_var = "trt",
#'     subgroup_vars = c("subgroup", "region"),
#'     response_type = "survival"
#'   )
#'
#'   print(effects)
#'   }
#' }
estimate_subgroup_effects <- function(brms_fit,
                                      original_data,
                                      trt_var,
                                      subgroup_vars = "auto",
                                      response_type = c("continuous", "binary", "count", "survival"),
                                      ndraws = NULL) {

  # --- 1. Validate inputs and determine which subgroups to analyze ---
  checkmate::assert_class(brms_fit, "brmsfit")
  checkmate::assert_data_frame(original_data, min.rows = 1)
  checkmate::assert_string(trt_var, min.chars = 1)

  # Assert trt_var is in both data frames
  checkmate::assert_subset(trt_var, names(brms_fit$data))
  checkmate::assert_subset(trt_var, names(original_data))

  # This is a complex check: subgroup_vars must be NULL, the string "auto", or a character vector.
  checkmate::assert(
    checkmate::check_string(subgroup_vars, pattern = "^auto$"),
    checkmate::check_character(subgroup_vars, null.ok = TRUE, min.len = 1, unique = TRUE))

  # If it is a character vector, check that all variables are in the original_data
  if (checkmate::test_character(subgroup_vars, null.ok = FALSE, min.len = 1) && !identical(subgroup_vars, "auto")) {
    checkmate::assert_subset(subgroup_vars, names(original_data))
  }

  # This replaces match.arg()
  response_type <- checkmate::assert_choice(response_type,
                                            choices = c("continuous", "binary", "count", "survival")
  )

  checkmate::assert_count(ndraws, null.ok = TRUE, positive = TRUE)

  prep <- .prepare_subgroup_vars(
    brms_fit = brms_fit,
    original_data = original_data,
    trt_var = trt_var,
    subgroup_vars = subgroup_vars
  )



  # --- 2. Create counterfactual datasets for "all treated" vs "all control" ---
  message("Step 1: Creating counterfactual datasets...")
  counterfactual_data <- .create_counterfactual_datasets(
    model_data = brms_fit$data,
    original_data = prep$data,
    trt_var = trt_var
  )

  # --- 3. Generate posterior predictions for each counterfactual scenario
  message("Step 2: Generating posterior predictions...")
  posterior_preds <- .get_posterior_predictions(
    brms_fit = brms_fit,
    data_control = counterfactual_data$control,
    data_treatment = counterfactual_data$treatment,
    response_type = response_type,
    original_data = prep$data,
    ndraws = ndraws
  )

  # --- 4. Calculate marginal effects for each subgroup based on posterior predictions
  message("Step 3: Calculating marginal effects...")
  results <- .calculate_and_summarize_effects(
    posterior_preds = posterior_preds,
    original_data = prep$data,
    subgroup_vars = prep$subgroup_vars,
    is_overall = prep$is_overall,
    response_type = response_type
  )

  message("Done.")
  return(results)
}


#' Prepare and Validate Subgroup Variables
#'
#' This helper function validates the subgroup variable inputs and handles the
#' "auto" detection logic. It parses interaction terms from the `brms_fit`
#' data (e.g., `subgroup_S1_x_trt`) to identify the base subgroup
#' variables (e.g., `subgroup`).
#'
#' If no subgroups are provided (`NULL`) or detected, it prepares the data
#' for an "Overall" marginal effect calculation by adding an "Overall" column.
#'
#' @param brms_fit The fitted `brmsfit` object.
#' @param original_data The original user-supplied data.frame.
#' @param trt_var The character string name of the treatment variable.
#' @param subgroup_vars The user's input: `NULL`, `"auto"`, or a character
#'   vector of subgroup variable names.
#'
#' @return A list containing:
#' \item{subgroup_vars}{A character vector of the validated or detected subgroup
#'   variable names (or "Overall").}
#' \item{data}{The `original_data`, potentially modified to include an
#'   "Overall" column.}
#' \item{is_overall}{A logical flag, `TRUE` if calculating only the overall effect.}
#'
#' @noRd
#'
.prepare_subgroup_vars <- function(brms_fit, original_data, trt_var, subgroup_vars) {
  # --- Assertions (as a safety net) ---
  checkmate::assert_class(brms_fit, "brmsfit")
  checkmate::assert_data_frame(original_data)
  checkmate::assert_string(trt_var)

  is_overall <- FALSE
  if (is.null(subgroup_vars)) {
    is_overall <- TRUE
  } else if (identical(subgroup_vars, "auto")) {
    message("`subgroup_vars` set to 'auto'. Detecting from model interaction terms...")
    interaction_cols <- names(brms_fit$data)[stringr::str_ends(names(brms_fit$data), paste0("\\_x_", trt_var))]

    if (length(interaction_cols) == 0) {
      message("...no interaction terms found. Calculating overall marginal effect.")
      is_overall <- TRUE
    } else {
      detected_vars <- interaction_cols %>%
        stringr::str_remove(paste0("\\_x_", trt_var, "$")) %>%
        stringr::str_replace("_[^_]*$", "")

      subgroup_vars <- unique(detected_vars)
      message(paste("...detected subgroup variable(s):", paste(subgroup_vars, collapse = ", ")))
    }
  }

  if (is_overall) {
    subgroup_vars <- "Overall"
    original_data$Overall <- "Overall"
  }

  return(list(
    subgroup_vars = subgroup_vars,
    data = original_data,
    is_overall = is_overall
  ))
}


#' Create Counterfactual Datasets
#'
#' Creates two "counterfactual" datasets based on the model data: one
#' representing an "all subjects received control" scenario and one representing
#' an "all subjects received treatment" scenario.
#'
#' @description
#' This function is crucial for marginal effect estimation.
#' 1.  Control Dataset: Sets the treatment variable to 0 and all related
#'     interaction dummy columns to 0.
#' 2.  Treatment Dataset: Sets the treatment variable to 1. Critically, it
#'     then reconstructs the dummy interaction columns (e.g.,
#'     `subgroup_S1_x_trt`). It uses the `original_data` (which contains the
#'     original factors) to determine which dummy should be 1, based on the
#'     subject's original subgroup.
#'
#' @param model_data The data.frame from the `brms_fit` object (`brms_fit$data`).
#'   This data has the dummy interaction columns.
#' @param original_data The original data.frame, which contains the
#'   categorical subgroup variables (e.g., `subgroup` as a factor).
#' @param trt_var The character string name of the treatment variable.
#'
#' @return A named list containing two data.frames: `control` and `treatment`.
#' @noRd
#'
.create_counterfactual_datasets <- function(model_data, original_data, trt_var) {

  # --- Assertions ---
  checkmate::assert_data_frame(model_data)
  checkmate::assert_data_frame(original_data)
  checkmate::assert_string(trt_var)

  # Check that data frames have the same number of rows
  if (nrow(model_data) != nrow(original_data)) {
    stop("model_data and original_data must have the same number of rows.")
  }


  # --- 1. Identify the dummy interaction columns ---
  interaction_cols <- names(model_data)[stringr::str_ends(names(model_data), paste0("_x_", trt_var))]

  # --- 2. Create the "All Control" Dataset ---
  data_control <- model_data
  data_control[[trt_var]] <- 0
  if (length(interaction_cols) > 0) {
    data_control[, interaction_cols] <- 0
  }

  # --- 3. Create the "All Treatment" Dataset ---
  data_treatment <- model_data
  data_treatment[[trt_var]] <- 1

  if (length(interaction_cols) > 0) {
    message("...setting interaction dummy variables for the 'all treatment' scenario.")

    # Pre-parse all interaction columns at once
    pattern <- paste0("^(.+)_([^_]+)_x_", trt_var, "$")
    parsed_interactions <- stringr::str_match(interaction_cols, pattern)

    # Create lookup tables for vectorized operations
    subgroup_vars <- parsed_interactions[, 2]
    levels <- parsed_interactions[, 3]
    names(subgroup_vars) <- interaction_cols
    names(levels) <- interaction_cols

    # Vectorized assignment instead of loop with repeated string operations
    for (col in interaction_cols) {
      subgroup_var <- subgroup_vars[col]
      level <- levels[col]

      if (!subgroup_var %in% names(original_data)) {
        warning(paste("Original subgroup variable '", subgroup_var, "' not found in data for column '", col, "'"))
        next
      }

      data_treatment[[col]] <- as.integer(original_data[[subgroup_var]] == level)
    }
  }

  return(list(control = data_control, treatment = data_treatment))
}

#' Get Posterior Predictions for Counterfactual Data
#'
#' This is the core prediction engine for marginal effects. It handles two
#' distinct model types:
#' 1.  Non-survival models: Calls `brms::posterior_epred` on the
#'     combined control and treatment datasets.
#' 2.  Survival models: This is more complex. It calls
#'     `brms::posterior_linpred` to get the linear predictor (`lp`),
#'     and manually reconstructs the cumulative baseline hazard (`H0`)
#'     from the model's spline parameters (`sbhaz`).
#'
#' @param brms_fit The fitted `brmsfit` object.
#' @param data_control The counterfactual "all control" dataset.
#' @param data_treatment The counterfactual "all treatment" dataset.
#' @param response_type The outcome type ("survival", "continuous", etc.).
#' @param original_data The original data, used for survival models to get
#'   event times and stratification levels.
#' @param ndraws An integer for the number of posterior draws, or `NULL`
#'   to use all draws.
#'
#' @return A list. The contents depend on `response_type`:
#' \itemize{
#'   \item For non-survival: `list(pred_control, pred_treatment)` where
#'     each is a [draws x observations] matrix from `posterior_epred`.
#'   \item For survival: `list(H0_posterior, linpred_control,
#'     linpred_treatment, strat_var, original_data)`. `H0_posterior`
#'     is a list of [draws x times] matrices (one per stratum), and the
#'     `linpred` objects are [draws x observations] matrices.
#' }
#' @noRd
#'
.get_posterior_predictions <- function(brms_fit, data_control, data_treatment, response_type, original_data, ndraws = NULL) {

  # --- Assertions ---
  checkmate::assert_class(brms_fit, "brmsfit")
  checkmate::assert_data_frame(data_control)
  checkmate::assert_data_frame(data_treatment)
  checkmate::assert_string(response_type)
  checkmate::assert_data_frame(original_data)
  checkmate::assert_count(ndraws, null.ok = TRUE)

  if (response_type == "survival") {
    message("... (reconstructing baseline hazard and getting linear predictors)...")

    # --- Parse formula to get details ---
    lhs_formula_str_vec <- deparse(brms_fit$formula$formula[[2]])
    bhaz_term <- paste(lhs_formula_str_vec, collapse = " ")
    resp_match <- stringr::str_match(bhaz_term, "(\\w+)\\s*\\|\\s*cens\\(1\\s*-\\s*(\\w+)\\)")
    if (is.na(resp_match[1,1])) {
      stop("Could not parse 'time | cens(1 - status)' structure from model formula.")
    }
    time_var <- resp_match[1, 2]
    status_var <- resp_match[1, 3]

    # Detect stratification
    strat_match <- stringr::str_match(bhaz_term, "gr\\s*=\\s*(\\w+)")
    strat_var <- if (!is.na(strat_match[1, 2])) strat_match[1, 2] else NULL

    # Extract spline components
    bknots_str <- stringr::str_extract(bhaz_term, "Boundary\\.knots = c\\(.*?\\)")
    knot_str <- stringr::str_extract(bhaz_term, "(?<!Boundary\\.)knots = c\\(.*?\\)")
    if (is.na(knot_str) || is.na(bknots_str)) {
      stop("Could not extract spline knot information from the model formula's bhaz() term.")
    }
    knots <- eval(parse(text = gsub("knots = ", "", knot_str)))
    bknots <- eval(parse(text = gsub("Boundary.knots = ", "", bknots_str)))

    # --- Reconstruct Baseline Hazard (H0) ---
    times_to_predict <- sort(unique(original_data[[time_var]][original_data[[status_var]] == 1]))

    # OPTIMIZATION: Use ndraws parameter if provided
    sbhaz_draws_df <- if (is.null(ndraws)) {
      posterior::as_draws_df(brms_fit, variable = "^sbhaz", regex = TRUE)
    } else {
      posterior::as_draws_df(brms_fit, variable = "^sbhaz", regex = TRUE, ndraws = ndraws)
    }

    all_sbhaz_names <- names(sbhaz_draws_df)[grep("^sbhaz", names(sbhaz_draws_df))]
    i_spline_basis <- splines2::iSpline(times_to_predict, knots = knots, Boundary.knots = bknots, intercept = FALSE)

    H0_posterior_list <- list()

    if (is.null(strat_var)) {
      # --- Non-stratified case ---
      sbhaz_matrix <- as.matrix(sbhaz_draws_df[, all_sbhaz_names])
      H0_posterior_list[["_default_"]] <- as.matrix(i_spline_basis %*% t(sbhaz_matrix))
    } else {
      # --- Stratified case ---
      message(paste("... model is stratified by '", strat_var, "'. Reconstructing hazard for each stratum.", sep = ""))
      strat_levels <- levels(factor(original_data[[strat_var]]))

      for (level in strat_levels) {
        param_pattern <- paste0("^sbhaz\\[", level, ",")
        level_sbhaz_cols <- grep(param_pattern, all_sbhaz_names, value = TRUE)

        if (length(level_sbhaz_cols) == 0) {
          warning(paste("Could not find sbhaz parameters for stratum '", level, "'. Skipping.", sep = ""))
          next
        }

        # Ensure correct numerical order of parameters
        numeric_indices <- as.numeric(stringr::str_extract(level_sbhaz_cols, "(?<=,)\\d+(?=\\$$)"))
        sorted_cols <- level_sbhaz_cols[order(numeric_indices)]

        sbhaz_matrix_level <- as.matrix(sbhaz_draws_df[, sorted_cols])

        if (ncol(i_spline_basis) != ncol(sbhaz_matrix_level)) {
          stop(paste0("Dimension mismatch for stratum '", level, "'."))
        }
        H0_posterior_list[[level]] <- as.matrix(i_spline_basis %*% t(sbhaz_matrix_level))
      }
    }

    # OPTIMIZATION: Combined dataset approach for linear predictors
    data_control$group <- "control"
    data_treatment$group <- "treatment"
    data_combined <- rbind(data_control, data_treatment)

    # Single call with ndraws parameter
    linpred_combined <- if (is.null(ndraws)) {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = NA)
    } else {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = NA, ndraws = ndraws)
    }

    # Split results
    n_control <- nrow(data_control)
    linpred_control <- linpred_combined[, 1:n_control]
    linpred_treatment <- linpred_combined[, (n_control + 1):ncol(linpred_combined)]

    return(list(
      H0_posterior = H0_posterior_list,
      linpred_control = linpred_control,
      linpred_treatment = linpred_treatment,
      strat_var = strat_var,
      original_data = original_data
    ))

  } else {
    message("... (predicting expected outcomes)...")
    # OPTIMIZATION: Combined dataset approach
    data_control$group <- "control"
    data_treatment$group <- "treatment"
    data_combined <- rbind(data_control, data_treatment)

    # Single prediction call with ndraws parameter
    pred_combined <- if (is.null(ndraws)) {
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = NA)
    } else {
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = NA, ndraws = ndraws)
    }

    # Split results
    n_control <- nrow(data_control)
    pred_control <- pred_combined[, 1:n_control]
    pred_treatment <- pred_combined[, (n_control + 1):ncol(pred_combined)]

    return(list(
      pred_control = pred_control,
      pred_treatment = pred_treatment
    ))
  }
}


#' Calculate and Summarize Marginal Effects
#'
#' This helper is the final aggregation step. It iterates through each
#' subgroup variable (e.g., "region") and each of its levels (e.g., "North"),
#' calculates the marginal treatment effect for that specific group, and
#' summarizes its posterior distribution.
#'
#' @description
#' The function operates as follows:
#' 1.  It loops through each variable in `subgroup_vars`.
#' 2.  For each variable, it loops through all its `levels`.
#' 3.  It finds all subjects (`subgroup_indices`) belonging to that level.
#' 4.  For non-survival models, it takes the `rowMeans` of predictions for
#'     those subjects to get a marginal posterior outcome for control and
#'     treatment, then computes the effect (e.g., difference, log-odds ratio, rate ratio).
#' 5.  For survival models, it passes the indices to
#'     `.calculate_survival_ahr_draws` to get the marginal effect draws.
#'
#' @param posterior_preds The list of posterior predictions from
#'   `.get_posterior_predictions`. Its structure differs for survival vs.
#'   non-survival models.
#' @param original_data The original data.frame, used to get subgroup factor levels.
#' @param subgroup_vars A character vector of subgroup variables to analyze.
#' @param is_overall A logical flag. If TRUE, only one "Overall" group is computed.
#' @param response_type The outcome type ("continuous", "binary", etc.). This
#'   controls the effect measure (e.g., difference vs. ratio).
#'
#' @return A named list with two elements:
#' \item{estimates}{A `tibble` with one row per subgroup level, containing
#'   the summary statistics (`Subgroup`, `Median`, `CI_Lower`, `CI_Upper`).}
#' \item{draws}{A wide `data.frame` (tibble) where each column represents the
#'   full posterior draws for a single subgroup effect.}
#'
#' @noRd
#'
.calculate_and_summarize_effects <- function(posterior_preds, original_data, subgroup_vars, is_overall, response_type) {
  # --- Assertions ---
  checkmate::assert_list(posterior_preds)
  checkmate::assert_data_frame(original_data)
  checkmate::assert_character(subgroup_vars, min.len = 1)
  checkmate::assert_logical(is_overall, len = 1)
  checkmate::assert_string(response_type)

  all_results_list <- list()
  all_draws_list <- list()

  # OPTIMIZATION: Pre-compute all factor information
  subgroup_factors <- list()
  if (!is_overall) {
    for (var in subgroup_vars) {
      factor_var <- as.factor(original_data[[var]])
      subgroup_factors[[var]] <- list(
        factor = factor_var,
        levels = levels(factor_var)
      )
    }
  }

  for (current_subgroup_var in subgroup_vars) {
    if (!is_overall) message(paste("... processing", current_subgroup_var))

    # Use pre-computed factor information
    if (is_overall) {
      current_data_subgroups <- as.factor(rep("Overall", nrow(original_data)))
      subgroup_levels <- "Overall"
    } else {
      current_data_subgroups <- subgroup_factors[[current_subgroup_var]]$factor
      subgroup_levels <- subgroup_factors[[current_subgroup_var]]$levels
    }

    level_results_list <- list()


    # Sequential processing for non-survival or single subgroup
    for (level in subgroup_levels) {
      subgroup_indices <- which(current_data_subgroups == level)

      # Calculate raw effect draws based on response type
      if (response_type == "survival") {
        effect_draws <- .calculate_survival_ahr_draws(
          linpred_control = posterior_preds$linpred_control,
          linpred_treatment = posterior_preds$linpred_treatment,
          H0_posterior_list = posterior_preds$H0_posterior,
          indices = subgroup_indices,
          strat_var = posterior_preds$strat_var,
          original_data = posterior_preds$original_data
        )
      } else {
        # OPTIMIZATION: Vectorized calculations for non-survival outcomes
        marginal_outcome_control <- if (length(subgroup_indices) == 1) {
          posterior_preds$pred_control[, subgroup_indices]
        } else {
          rowMeans(posterior_preds$pred_control[, subgroup_indices, drop = FALSE])
        }

        marginal_outcome_treatment <- if (length(subgroup_indices) == 1) {
          posterior_preds$pred_treatment[, subgroup_indices]
        } else {
          rowMeans(posterior_preds$pred_treatment[, subgroup_indices, drop = FALSE])
        }

        effect_draws <- switch(
          response_type,
          continuous = marginal_outcome_treatment - marginal_outcome_control,
          binary = qlogis(marginal_outcome_treatment) - qlogis(marginal_outcome_control),
          count = marginal_outcome_treatment / marginal_outcome_control
        )
      }

      # Summarize and store draws
      subgroup_name <- if (is_overall) "Overall" else paste0(current_subgroup_var, ": ", level)
      all_draws_list[[subgroup_name]] <- effect_draws

      point_estimate <- median(effect_draws, na.rm = TRUE)
      ci <- quantile(effect_draws, probs = c(0.025, 0.975), na.rm = TRUE)

      level_results_list[[level]] <- tibble::tibble(
        Subgroup = subgroup_name,
        Median = point_estimate,
        CI_Lower = ci[1],
        CI_Upper = ci[2]
      )

    }

    all_results_list[[current_subgroup_var]] <- dplyr::bind_rows(level_results_list)
  }

  final_results <- dplyr::bind_rows(all_results_list)
  draws_df <- dplyr::bind_cols(all_draws_list)

  return(list(estimates = final_results, draws = draws_df))
}


#' Calculate Average Hazard Ratio (AHR) Draws for Survival Models
#'
#' This function orchestrates the calculation of the Average Hazard Ratio (AHR)
#' for a specific subgroup. It is designed to be memory-efficient by processing
#' the posterior draws in chunks.
#'
#' @description
#' For each chunk of posterior draws, this function:
#' 1.  Calls `.get_marginal_survival_vectorized` to get the marginal survival
#'     curve for the control and treatment scenarios.
#' 2.  Calls `.calculate_ahr_vectorized` to compute the AHR from those
#'     two survival curves.
#'
#' @param linpred_control A [draws x observations] matrix of linear predictors
#'   for the control scenario.
#' @param linpred_treatment A [draws x observations] matrix of linear predictors
#'   for the treatment scenario.
#' @param H0_posterior_list A list (named by stratum) of [draws x times]
#'   matrices for the cumulative baseline hazard.
#' @param indices A numeric vector of row indices specifying which subjects
#'   belong to the current subgroup.
#' @param strat_var The character string name of the stratification variable,
#'   or `NULL` if not stratified.
#' @param original_data The original data.frame, passed down to the helper
#'   to match subjects in `indices` to their stratum.
#'
#' @return A numeric vector of length `n_draws`, where each element is the
#'   calculated AHR for that posterior draw.
#' @noRd
#'
.calculate_survival_ahr_draws <- function(linpred_control, linpred_treatment, H0_posterior_list, indices, strat_var, original_data) {
  # --- Assertions ---
  checkmate::assert_matrix(linpred_control)
  checkmate::assert_matrix(linpred_treatment)
  checkmate::assert_list(H0_posterior_list, names = "named")
  checkmate::assert_integerish(indices, min.len = 1, unique = TRUE)
  checkmate::assert_string(strat_var, null.ok = TRUE)
  checkmate::assert_data_frame(original_data)

  # Check dimensions
  if (nrow(linpred_control) != nrow(linpred_treatment)) {
    stop("linpred_control and linpred_treatment must have the same number of rows (draws).")
  }

  # OPTIMIZATION: Process in chunks to reduce memory usage for large datasets
  n_draws <- nrow(linpred_control)
  chunk_size <- min(1000, n_draws)
  n_chunks <- ceiling(n_draws / chunk_size)

  ahr_draws <- numeric(n_draws)

  for (chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_draws)
    chunk_indices <- start_idx:end_idx

    # Get survival curves for this chunk
    S_control_chunk <- .get_marginal_survival_vectorized(
      linpred_control[chunk_indices, , drop = FALSE],
      H0_posterior_list, indices, strat_var, original_data
    )
    S_treatment_chunk <- .get_marginal_survival_vectorized(
      linpred_treatment[chunk_indices, , drop = FALSE],
      H0_posterior_list, indices, strat_var, original_data
    )

    # Calculate AHR for this chunk using vectorized operations
    ahr_draws[chunk_indices] <- .calculate_ahr_vectorized(S_control_chunk, S_treatment_chunk)
  }

  return(ahr_draws)
}

#' Calculate Vectorized Marginal Survival Curve
#'
#' Calculates the marginal survival curve S(t) for a specific subgroup
#' for each posterior draw, using vectorized operations.
#'
#' @description
#' This function implements the survival calculation at the core of the AHR estimation.
#' It uses two different methods depending on stratification:
#'
#' 1.  Non-stratified: It computes S(t) using the "average hazard" method:
#'     $S(t|\text{draw}) = \exp(-H_0(t|\text{draw}) \times \text{mean}(\exp(\eta_{\text{subgroup}}|\text{draw})))$
#'
#' 2.  Stratified: It calculates the individual survival curve for each
#'     subject in the subgroup, using their specific stratum's baseline hazard:
#'     $S_i(t) = \exp(-H_{0, \text{stratum}_i}(t|\text{draw}) \times \exp(\eta_i|\text{draw}))$
#'     It then computes the marginal survival by averaging these individual curves:
#'     $S(t|\text{draw}) = \text{mean}(S_i(t|\text{draw}))$
#'
#' @param linpred_posterior A [draws x total_observations] matrix of linear predictors.
#' @param H0_post_list A list (named by stratum, or "_default_") of [times x draws]
#'   matrices for the cumulative baseline hazard.
#' @param sub_indices A numeric vector of row indices for the current subgroup.
#' @param strat_variable The character string name of the stratification variable,
#'   or `NULL`.
#' @param full_data The original data.frame, used to map subjects in
#'   `sub_indices` to their stratum.
#'
#' @return A [draws x times] matrix, where each row is the marginal survival
#'   curve for the subgroup for that posterior draw.
#' @noRd
#'
.get_marginal_survival_vectorized <- function(linpred_posterior, H0_post_list, sub_indices, strat_variable, full_data) {
  # --- Assertions ---
  checkmate::assert_matrix(linpred_posterior)
  checkmate::assert_list(H0_post_list, names = "named")
  checkmate::assert_integerish(sub_indices, min.len = 1, unique = TRUE)
  checkmate::assert_string(strat_variable, null.ok = TRUE)
  checkmate::assert_data_frame(full_data)

  subgroup_linpred <- linpred_posterior[, sub_indices, drop = FALSE]
  n_draws <- nrow(subgroup_linpred)

  # Non-stratified case
  if (is.null(strat_variable)) {
    H0_post <- H0_post_list[["_default_"]]
    n_times <- nrow(H0_post)

    # OPTIMIZATION: Vectorized matrix operations instead of loops
    exp_eta <- exp(subgroup_linpred)  # n_draws x n_subjects

    # Calculate marginal survival using matrix operations
    S_marginal <- matrix(NA, nrow = n_draws, ncol = n_times)
    for (i in 1:n_draws) {
      h0_i <- H0_post[, i]
      eta_exp_i <- exp_eta[i, ]
      # Vectorized calculation: exp(-h0 * mean(exp(eta)))
      S_marginal[i, ] <- exp(-h0_i * mean(eta_exp_i))
    }

    return(S_marginal)
  }

  # Stratified case
  subgroup_strata <- full_data[[strat_variable]][sub_indices]
  unique_strata_in_subgroup <- unique(subgroup_strata)
  n_times <- nrow(H0_post_list[[unique_strata_in_subgroup[1]]])

  S_marginal <- matrix(NA, nrow = n_draws, ncol = n_times)

  # OPTIMIZATION: Pre-compute stratum indices
  stratum_indices <- lapply(unique_strata_in_subgroup, function(s) which(subgroup_strata == s))
  names(stratum_indices) <- unique_strata_in_subgroup

  for (i in 1:n_draws) {
    S_individual_by_time <- matrix(NA, nrow = n_times, ncol = length(sub_indices))

    for (stratum in unique_strata_in_subgroup) {
      indices_this_stratum <- stratum_indices[[stratum]]
      H0_post_stratum <- H0_post_list[[stratum]]

      if (is.null(H0_post_stratum)) {
        warning(paste("No baseline hazard found for stratum:", stratum))
        next
      }

      eta_i_stratum <- subgroup_linpred[i, indices_this_stratum, drop = FALSE]
      h0_i_stratum <- H0_post_stratum[, i]

      # Vectorized calculation for this stratum
      S_individual_by_time[, indices_this_stratum] <- exp(-outer(h0_i_stratum, exp(eta_i_stratum)))
    }

    # Calculate marginal survival curve
    S_marginal[i, ] <- rowMeans(S_individual_by_time, na.rm = TRUE)
  }

  return(S_marginal)
}

#' Calculate Average Hazard Ratio from Survival Curves
#'
#' Calculates the Average Hazard Ratio (AHR) for each posterior draw, given
#' paired control and treatment survival curves.
#'
#' @description
#' This function iterates row-by-row (i.e., per posterior draw) and applies
#' the formula for the Average Hazard Ratio (AHR). This involves:
#' 1.  Approximating the discrete hazard (event probability) for each time
#'     interval, $h(t)$, as the drop in survival: $h(t) \approx S(t-1) - S(t)$.
#' 2.  Calculating the AHR as the ratio of expected hazards:
#'     $\frac{\sum S_{\text{control}}(t) \times h_{\text{treatment}}(t)}{\sum S_{\text{treatment}}(t) \times h_{\text{control}}(t)}$
#'
#' @param S_control A [draws x times] matrix of marginal survival
#'   probabilities for the control group. Each row is one draw, each column is a
#'   time point.
#' @param S_treatment A [draws x times] matrix of marginal survival
#'   probabilities for the treatment group.
#'
#' @return A numeric vector of length `n_draws`, where each element is the
#'   calculated AHR for that posterior draw.
#' @noRd
#'
.calculate_ahr_vectorized <- function(S_control, S_treatment) {
  # --- Assertions ---
  checkmate::assert_matrix(S_control)
  checkmate::assert_matrix(S_treatment)

  if (any(dim(S_control) != dim(S_treatment))) {
    stop("S_control and S_treatment matrices must have identical dimensions.")
  }
  n_draws <- nrow(S_control)
  ahr_draws <- numeric(n_draws)

  for (i in 1:n_draws) {
    sc <- S_control[i, ]
    st <- S_treatment[i, ]

    # Calculate differences (hazard approximation)
    dsc <- -diff(c(1, sc))
    dst <- -diff(c(1, st))

    # Handle edge cases
    if (length(dsc) != length(dst) || length(dsc) == 0) {
      ahr_draws[i] <- NA
      next
    }

    # Pad to same length if needed
    if (length(dsc) < length(sc)) dsc <- c(dsc, 0)
    if (length(dst) < length(st)) dst <- c(dst, 0)

    # Calculate AHR
    numerator <- sum(sc * dst, na.rm = TRUE)
    denominator <- sum(st * dsc, na.rm = TRUE)

    ahr_draws[i] <- if (denominator == 0 || is.na(denominator)) NA else numerator / denominator
  }

  return(ahr_draws)
}

