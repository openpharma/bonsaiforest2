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
#' @param brms_fit A fitted `brmsfit` object from `fit_brms_model()` or
#'    `run_brms_analysis()`.
#' @param original_data A `data.frame` containing the data. While this parameter
#'    is called "original_data" for backward compatibility, the function internally
#'    uses the processed data from `brms_fit$data` to ensure consistency with the
#'    model's factor coding and contrasts. This parameter is mainly used for
#'    validation purposes.
#' @param trt_var A character string specifying the name of the treatment variable.
#' @param subgroup_vars A character vector of subgroup variable names found in
#'    `original_data`. If set to `"auto"` (the default), the function
#'    attempts to automatically identify subgroup variables from the model's
#'    formula.
#' @param response_type The type of outcome variable. One of "binary", "count",
#'    "continuous", or "survival".
#' @param ndraws An integer specifying the number of posterior draws to use.
#'    If `NULL` (default), all available draws are used.
#'
#' @return
#' A `data.frame` (tibble) where each row corresponds to a subgroup
#' (or the "Overall" effect), providing the estimated marginal effect and
#' posterior summaries (e.g., `mean`, `sd`, `q2.5`, `q97.5`).
#'
#' @importFrom checkmate assert_class assert_data_frame assert_string assert_choice
#' @importFrom checkmate assert_character assert_count assert_logical assert_subset
#' @importFrom stringr str_detect str_match_all
#' @importFrom dplyr bind_rows bind_cols
#' @importFrom stats contrasts contrasts<-
#' @importFrom tibble tibble
#' @export
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

  # Check that trt_var exists in the model data
  checkmate::assert_subset(trt_var, names(brms_fit$data))

  # Check response type
  response_type <- checkmate::assert_choice(response_type,
                                            choices = c("continuous", "binary", "count", "survival"))
  checkmate::assert_count(ndraws, null.ok = TRUE, positive = TRUE)

  # CRITICAL FIX: Use brms_fit$data (processed data with contrasts) for all operations
  # This ensures consistency between predictions and subgroup membership
  # The original_data parameter is kept for backward compatibility but we use model data
  model_data <- brms_fit$data

  # Validate that model_data and original_data are compatible (same number of rows)
  if (nrow(model_data) != nrow(original_data)) {
    warning("original_data has different number of rows than the fitted model data. ",
            "Using the data from brms_fit object for consistency.")
  }

  # --- 2. Create counterfactual datasets (UPDATED for Factor Logic) ---
  # Use model data (with contrasts) for everything
  message("Step 1: Creating counterfactual datasets...")

  # Prepare Subgroup Variables - use model_data for consistency
  prep <- .prepare_subgroup_vars(
    brms_fit = brms_fit,
    original_data = model_data,  # Use model data, not user-provided original
    trt_var = trt_var,
    subgroup_vars = subgroup_vars
  )

  counterfactual_data <- .create_counterfactual_datasets(
    model_data = prep$data,  # Use prep$data which may have "Overall" added
    trt_var = trt_var
  )

  # --- 3. Generate posterior predictions ---
  message("Step 2: Generating posterior predictions...")
  posterior_preds <- .get_posterior_predictions(
    brms_fit = brms_fit,
    data_control = counterfactual_data$control,
    data_treatment = counterfactual_data$treatment,
    response_type = response_type,
    original_data = prep$data, # Use prep$data (might have Overall added)
    ndraws = ndraws
  )

  # --- 4. Calculate marginal effects ---
  message("Step 3: Calculating marginal effects...")
  results <- .calculate_and_summarize_effects(
    posterior_preds = posterior_preds,
    original_data = prep$data,  # Same data used for subgroup membership
    subgroup_vars = prep$subgroup_vars,
    is_overall = prep$is_overall,
    response_type = response_type
  )

  message("Done.")
  return(results)
}

#' Prepare and Validate Subgroup Variables (Formula Parsing)
#'
#' Validates subgroup inputs. If "auto", it scans the brms formula object
#' for interaction terms involving the treatment variable (e.g. `trt:subgroup`).
#'
#' @noRd
.prepare_subgroup_vars <- function(brms_fit, original_data, trt_var, subgroup_vars) {

  # Logic to handle "auto" or specific columns
  checkmate::assert(
    checkmate::check_string(subgroup_vars, pattern = "^auto$"),
    checkmate::check_character(subgroup_vars, null.ok = TRUE, min.len = 1, unique = TRUE)
  )

  is_overall <- FALSE

  if (is.null(subgroup_vars)) {
    is_overall <- TRUE
    subgroup_vars <- "Overall"
  } else if (identical(subgroup_vars, "auto")) {
    message("`subgroup_vars` set to 'auto'. Detecting from model data...")

    # With the new explicit dummy approach, interaction dummies follow pattern: trt_VARNAMElevel
    # We look for columns in model_data that match this pattern
    dummy_pattern <- paste0("^", trt_var, "_")
    interaction_cols <- grep(dummy_pattern, names(original_data), value = TRUE)

    if (length(interaction_cols) > 0) {
      # Extract variable names by removing the trt_ prefix and the level suffix
      # We check which factor variables in the data could have produced these dummies
      detected_vars <- character(0)

      for (var_name in names(original_data)) {
        if (is.factor(original_data[[var_name]]) && var_name != trt_var) {
          # Check if any dummy matches this variable
          var_levels <- levels(original_data[[var_name]])
          # Use make.names to match the sanitized names created in prepare_formula_model
          expected_dummies <- make.names(paste0(trt_var, "_", var_name, var_levels), unique = FALSE)
          if (any(expected_dummies %in% interaction_cols)) {
            detected_vars <- c(detected_vars, var_name)
          }
        }
      }
      detected_vars <- unique(detected_vars)
    } else {
      detected_vars <- character(0)
    }

    if (length(detected_vars) == 0) {
      message("...no specific interaction terms detected. Calculating overall effect.")
      is_overall <- TRUE
      subgroup_vars <- "Overall"
    } else {
      subgroup_vars <- detected_vars
      message(paste("...detected subgroup variable(s):", paste(subgroup_vars, collapse = ", ")))
    }
  } else {
    # Validate user provided vars exist
    checkmate::assert_subset(subgroup_vars, names(original_data))
  }

  if (is_overall) {
    original_data$Overall <- "Overall"
  }

  return(list(
    subgroup_vars = subgroup_vars,
    data = original_data,
    is_overall = is_overall
  ))
}

#' Create Counterfactual Datasets (With Explicit Dummy Recreation)
#'
#' Creates "all control" and "all treatment" datasets.
#' Importantly, this function also recreates any interaction dummy variables
#' (e.g., trt_regionA, trt_regionB) that were created by prepare_formula_model.
#'
#' @noRd
.create_counterfactual_datasets <- function(model_data, trt_var) {
  checkmate::assert_data_frame(model_data)
  checkmate::assert_string(trt_var)

  # Ensure treatment is a factor to get levels
  if (!is.factor(model_data[[trt_var]])) {
     # Fallback if it was somehow converted to int, though prepare_formula handles this
     model_data[[trt_var]] <- as.factor(model_data[[trt_var]])
  }

  trt_levels <- levels(model_data[[trt_var]])
  if(length(trt_levels) < 2) stop("Treatment variable must have at least 2 levels.")

  # Store original contrasts to preserve them
  trt_contrasts <- contrasts(model_data[[trt_var]])

  ref_level <- trt_levels[1] # Usually 0
  alt_level <- trt_levels[2] # Usually 1

  # Identify interaction dummy columns (pattern: trt_varname)
  # These are created by .process_predictive_terms as trt_regionA, trt_regionB, etc.
  interaction_dummy_pattern <- paste0("^", trt_var, "_")
  interaction_cols <- grep(interaction_dummy_pattern, names(model_data), value = TRUE)

  # Create Control (Reference) Data
  data_control <- model_data
  data_control[[trt_var]] <- factor(rep(ref_level, nrow(model_data)), levels = trt_levels)
  contrasts(data_control[[trt_var]]) <- trt_contrasts

  # Set all interaction dummies to 0 (since trt = control)
  for (col in interaction_cols) {
    data_control[[col]] <- 0
  }

  # Create Treatment Data
  data_treatment <- model_data
  data_treatment[[trt_var]] <- factor(rep(alt_level, nrow(model_data)), levels = trt_levels)
  contrasts(data_treatment[[trt_var]]) <- trt_contrasts

  # Recreate interaction dummies for treatment arm
  # Pattern: trt_subgroupvarLEVEL, e.g., trt_regionA
  # Value should be 1 if patient is in that subgroup level, 0 otherwise
  for (col in interaction_cols) {
    # Find which variable this corresponds to by checking which factor variable
    # in model_data has this level
    matched <- FALSE
    for (var_name in names(model_data)) {
      if (is.factor(model_data[[var_name]]) && var_name != trt_var) {
        # Check if this dummy name contains a level of this factor
        var_levels <- levels(model_data[[var_name]])
        for (level in var_levels) {
          expected_dummy_name <- paste0(trt_var, "_", var_name, level)
          if (col == expected_dummy_name) {
            # This is the dummy for this variable-level combination
            # Set to 1 for patients in this level, 0 otherwise
            data_treatment[[col]] <- as.numeric(model_data[[var_name]] == level)
            matched <- TRUE
            break
          }
        }
        if (matched) break
      }
    }
    
    # If no match found, the dummy column might be from a variable that's not a factor
    # in the current data. Keep it at 0 for safety.
    if (!matched) {
      warning(paste("Could not match interaction column", col, "to any factor variable. Setting to 0."))
      data_treatment[[col]] <- 0
    }
  }

  # Verify the counterfactuals were created correctly
  if (!all(data_control[[trt_var]] == ref_level)) {
    stop("Control counterfactual data creation failed - not all patients set to control")
  }
  if (!all(data_treatment[[trt_var]] == alt_level)) {
    stop("Treatment counterfactual data creation failed - not all patients set to treatment")
  }

  return(list(control = data_control, treatment = data_treatment))
}

#' Get Posterior Predictions (Robust)
#'
#' Generates predictions. Handles Survival (manual bhaz reconstruction)
#' and Standard (posterior_epred) models.
#'
#' @noRd
.get_posterior_predictions <- function(brms_fit, data_control, data_treatment, response_type, original_data, ndraws = NULL) {

  checkmate::assert_class(brms_fit, "brmsfit")

  # OPTIMIZATION: Combine datasets for a single call to brms
  # This reduces overhead of Stan I/O
  data_control$._sim_group <- "control"
  data_treatment$._sim_group <- "treatment"
  data_combined <- rbind(data_control, data_treatment)

  if (response_type == "survival") {
    message("... (reconstructing baseline hazard and getting linear predictors)...")

    # 1. Get Linear Predictors (eta) for combined data
    # brms handles the trt:subgroup interaction here automatically
    linpred_combined <- if (is.null(ndraws)) {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = NA)
    } else {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = NA, ndraws = ndraws)
    }

    n_control <- nrow(data_control)
    linpred_control <- linpred_combined[, 1:n_control]
    linpred_treatment <- linpred_combined[, (n_control + 1):ncol(linpred_combined)]

    # 2. Reconstruct Baseline Hazard
    # (Reuse the robust logic from your previous code)
    h0_res <- .extract_baseline_hazard(brms_fit, original_data, ndraws)

    return(list(
      H0_posterior = h0_res$H0_posterior,
      linpred_control = linpred_control,
      linpred_treatment = linpred_treatment,
      strat_var = h0_res$strat_var,
      original_data = original_data
    ))

  } else {
    message("... (predicting expected outcomes)...")

    pred_combined <- if (is.null(ndraws)) {
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = NA)
    } else {
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = NA, ndraws = ndraws)
    }

    n_control <- nrow(data_control)
    pred_control <- pred_combined[, 1:n_control]
    pred_treatment <- pred_combined[, (n_control + 1):ncol(pred_combined)]

    return(list(
      pred_control = pred_control,
      pred_treatment = pred_treatment
    ))
  }
}

#' Helper: Extract Baseline Hazard details for Survival
#'
#' Separated for cleanliness. Parses formula to find bhaz() terms.
#' @noRd
.extract_baseline_hazard <- function(brms_fit, original_data, ndraws) {
    # Parse formula for bhaz
    lhs_formula_str_vec <- deparse(brms_fit$formula$formula[[2]])
    bhaz_term <- paste(lhs_formula_str_vec, collapse = " ")

    resp_match <- stringr::str_match(bhaz_term, "(\\w+)\\s*\\|\\s*cens\\(1\\s*-\\s*(\\w+)\\)")
    if (is.na(resp_match[1,1])) stop("Could not parse 'time | cens(1 - status)' structure.")

    time_var <- resp_match[1, 2]
    status_var <- resp_match[1, 3]

    strat_match <- stringr::str_match(bhaz_term, "gr\\s*=\\s*(\\w+)")
    strat_var <- if (!is.na(strat_match[1, 2])) strat_match[1, 2] else NULL

    bknots_str <- stringr::str_extract(bhaz_term, "Boundary\\.knots = c\\(.*?\\)")
    knot_str <- stringr::str_extract(bhaz_term, "(?<!Boundary\\.)knots = c\\(.*?\\)")

    knots <- eval(parse(text = gsub("knots = ", "", knot_str)))
    bknots <- eval(parse(text = gsub("Boundary.knots = ", "", bknots_str)))

    times_to_predict <- sort(unique(original_data[[time_var]][original_data[[status_var]] == 1]))

    sbhaz_draws_df <- if (is.null(ndraws)) {
      posterior::as_draws_df(brms_fit, variable = "^sbhaz", regex = TRUE)
    } else {
      posterior::as_draws_df(brms_fit, variable = "^sbhaz", regex = TRUE, ndraws = ndraws)
    }

    all_sbhaz_names <- names(sbhaz_draws_df)[grep("^sbhaz", names(sbhaz_draws_df))]
    i_spline_basis <- splines2::iSpline(times_to_predict, knots = knots, Boundary.knots = bknots, intercept = FALSE)

    H0_posterior_list <- list()

    if (is.null(strat_var)) {
      sbhaz_matrix <- as.matrix(sbhaz_draws_df[, all_sbhaz_names])
      H0_posterior_list[["_default_"]] <- as.matrix(i_spline_basis %*% t(sbhaz_matrix))
    } else {
      strat_levels <- levels(factor(original_data[[strat_var]]))
      for (level in strat_levels) {
        param_pattern <- paste0("^sbhaz\\[", level, ",")
        level_sbhaz_cols <- grep(param_pattern, all_sbhaz_names, value = TRUE)
        # Sort cols numerically
        numeric_indices <- as.numeric(stringr::str_extract(level_sbhaz_cols, "(?<=,)\\d+(?=\\$$)"))
        sorted_cols <- level_sbhaz_cols[order(numeric_indices)]

        if(length(sorted_cols) > 0) {
            sbhaz_matrix_level <- as.matrix(sbhaz_draws_df[, sorted_cols])
            H0_posterior_list[[level]] <- as.matrix(i_spline_basis %*% t(sbhaz_matrix_level))
        }
      }
    }

    return(list(H0_posterior = H0_posterior_list, strat_var = strat_var))
}


#' Calculate and Summarize Marginal Effects (Identical Logic to Previous)
#' @noRd
.calculate_and_summarize_effects <- function(posterior_preds, original_data, subgroup_vars, is_overall, response_type) {

  all_results_list <- list()
  all_draws_list <- list()

  # Pre-compute factor info for iteration
  subgroup_factors <- list()
  if (!is_overall) {
    for (var in subgroup_vars) {
      factor_var <- as.factor(original_data[[var]])
      subgroup_factors[[var]] <- list(factor = factor_var, levels = levels(factor_var))
    }
  }

  for (current_subgroup_var in subgroup_vars) {
    if (!is_overall) message(paste("... processing", current_subgroup_var))

    if (is_overall) {
      current_data_subgroups <- as.factor(rep("Overall", nrow(original_data)))
      subgroup_levels <- "Overall"
    } else {
      current_data_subgroups <- subgroup_factors[[current_subgroup_var]]$factor
      subgroup_levels <- subgroup_factors[[current_subgroup_var]]$levels
    }

    level_results_list <- list()

    for (level in subgroup_levels) {
      subgroup_indices <- which(current_data_subgroups == level)

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
        # Vectorized means
        # Check dimensions: rows=draws, cols=obs
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
#' (Unchanged from your logic, included for completeness)
#' @noRd
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
#' @noRd
.get_marginal_survival_vectorized <- function(linpred_posterior, H0_post_list, sub_indices, strat_variable, full_data) {
  # ... [Same as your previous implementation] ...
  subgroup_linpred <- linpred_posterior[, sub_indices, drop = FALSE]
  n_draws <- nrow(subgroup_linpred)

  if (is.null(strat_variable)) {
    H0_post <- H0_post_list[["_default_"]]
    exp_eta <- exp(subgroup_linpred)
    S_marginal <- matrix(NA, nrow = n_draws, ncol = nrow(H0_post))
    for (i in 1:n_draws) {
       S_marginal[i, ] <- exp(-H0_post[, i] * mean(exp_eta[i, ]))
    }
    return(S_marginal)
  }

  # Stratified logic
  subgroup_strata <- full_data[[strat_variable]][sub_indices]
  unique_strata <- unique(subgroup_strata)
  n_times <- nrow(H0_post_list[[unique_strata[1]]])
  S_marginal <- matrix(NA, nrow = n_draws, ncol = n_times)

  for (i in 1:n_draws) {
    S_indiv <- matrix(NA, nrow = n_times, ncol = length(sub_indices))
    for (s in unique_strata) {
       idx <- which(subgroup_strata == s)
       h0 <- H0_post_list[[s]][, i]
       eta <- subgroup_linpred[i, idx]
       S_indiv[, idx] <- exp(-outer(h0, exp(eta)))
    }
    S_marginal[i, ] <- rowMeans(S_indiv, na.rm=TRUE)
  }
  return(S_marginal)
}

#' Calculate AHR Vectorized
#' @noRd
.calculate_ahr_vectorized <- function(S_control, S_treatment) {
  # ... [Same as your previous implementation] ...
  n_draws <- nrow(S_control)
  ahr_draws <- numeric(n_draws)
  for (i in 1:n_draws) {
    sc <- S_control[i, ]
    st <- S_treatment[i, ]
    dsc <- -diff(c(1, sc))
    dst <- -diff(c(1, st))

    # Pad
    if(length(dsc) < length(sc)) dsc <- c(dsc, 0)
    if(length(dst) < length(st)) dst <- c(dst, 0)

    num <- sum(sc * dst, na.rm=TRUE)
    den <- sum(st * dsc, na.rm=TRUE)
    ahr_draws[i] <- if (den == 0) NA else num/den
  }
  return(ahr_draws)
}


