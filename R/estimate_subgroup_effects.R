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
#' @param original_data A `data.frame` containing the data. Used for validation
#'    and structure, while internal predictions use the model's processed data.
#' @param trt_var A character string specifying the name of the treatment variable.
#' @param subgroup_vars A character vector of subgroup variable names.
#'    If "auto", detects both fixed interaction terms (Colon syntax) and
#'    random effect grouping factors (Pipe syntax).
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

  # --- 2. Prepare Subgroup Variables (UPDATED for Pipe Support) ---
  message("Step 1: Identifying subgroups and creating counterfactuals...")

  prep <- .prepare_subgroup_vars(
    brms_fit = brms_fit,
    original_data = model_data,
    trt_var = trt_var,
    subgroup_vars = subgroup_vars
  )

  # --- 3. Create Counterfactuals ---
  # Note: For Pipe syntax, this function essentially just sets the trt factor.
  # The 'subgroup' factor remains as is, which is exactly what brms needs
  # to look up the random effect.
  counterfactual_data <- .create_counterfactual_datasets(
    model_data = prep$data,
    trt_var = trt_var
  )

  # --- 4. Generate posterior predictions (UPDATED for RE support) ---
  message("Step 2: Generating posterior predictions...")
  posterior_preds <- .get_posterior_predictions(
    brms_fit = brms_fit,
    data_control = counterfactual_data$control,
    data_treatment = counterfactual_data$treatment,
    response_type = response_type,
    original_data = prep$data,
    ndraws = ndraws,
    trt_var = trt_var
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

#' Prepare and Validate Subgroup Variables
#'
#' Now detects BOTH:
#' 1. Fixed Interactions (dummies starting with trt_var_)
#' 2. Random Effects Grouping Factors (from Pipe syntax)
#' @noRd
.prepare_subgroup_vars <- function(brms_fit, original_data, trt_var, subgroup_vars) {

  checkmate::assert(
    checkmate::check_string(subgroup_vars, pattern = "^auto$"),
    checkmate::check_character(subgroup_vars, null.ok = TRUE, min.len = 1, unique = TRUE)
  )

  is_overall <- FALSE

  if (is.null(subgroup_vars)) {
    is_overall <- TRUE
    subgroup_vars <- "Overall"
  } else if (identical(subgroup_vars, "auto")) {
    message("`subgroup_vars` set to 'auto'. Detecting from model...")
    detected_vars <- character(0)

    # A. Detect FIXED interactions (Colon syntax)
    # -------------------------------------------
    dummy_pattern <- paste0("^", trt_var, "_")
    interaction_cols <- grep(dummy_pattern, names(original_data), value = TRUE)

    if (length(interaction_cols) > 0) {
      for (var_name in names(original_data)) {
        if (is.factor(original_data[[var_name]]) && var_name != trt_var) {
          # Check if dummies match this variable
          var_levels <- levels(original_data[[var_name]])
          expected_dummies <- make.names(paste0(trt_var, "_", var_name, var_levels), unique = FALSE)
          if (any(expected_dummies %in% interaction_cols)) {
            detected_vars <- c(detected_vars, var_name)
          }
        }
      }
    }

    # B. Detect RANDOM effects grouping factors (Pipe syntax)
    # -----------------------------------------------------
    # brms::ranef() returns a list where names are the grouping factors
    # We check if the random effects for a group include the treatment slope
    # Only attempt this if the model actually has random effects
    re_structure <- tryCatch(
      brms::ranef(brms_fit, summary = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(re_structure)) {
      # re_structure is a list of arrays. Names of list are grouping vars (e.g. "subgroup")
      grouping_factors <- names(re_structure)

      for (g_var in grouping_factors) {
        # We only care if this grouping factor is in our data (it should be)
        if (g_var %in% names(original_data)) {
           # Check if 'trt' is a slope for this group
           # The dimnames of the array usually contain the coef names
           # e.g., dimnames(re_structure$subgroup)[[3]] might contain "trt1"
           re_coefs <- dimnames(re_structure[[g_var]])[[3]] # The 3rd dim is the parameter name

           # Robust check: Does any coefficient name START with the treatment variable?
           # This handles trt1, trt_1, trt_control, etc. across different coding schemes
           if (any(grepl(paste0("^", trt_var), re_coefs))) {
              detected_vars <- c(detected_vars, g_var)
           }
        }
      }
    }

    detected_vars <- unique(detected_vars)

    if (length(detected_vars) == 0) {
      message("...no subgroup terms detected. Calculating overall effect.")
      is_overall <- TRUE
      subgroup_vars <- "Overall"
    } else {
      subgroup_vars <- detected_vars
      message(paste("...detected subgroup variable(s):", paste(subgroup_vars, collapse = ", ")))
    }
  } else {
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

#' Create Counterfactual Datasets
#' Unchanged, but commented to explain why it works for Pipe syntax.
#' @noRd
.create_counterfactual_datasets <- function(model_data, trt_var) {
  checkmate::assert_data_frame(model_data)

  if (!is.factor(model_data[[trt_var]])) model_data[[trt_var]] <- as.factor(model_data[[trt_var]])
  trt_levels <- levels(model_data[[trt_var]])
  trt_contrasts <- contrasts(model_data[[trt_var]])

  ref_level <- trt_levels[1]
  alt_level <- trt_levels[2]

  # 1. Identify Fixed Interaction Dummies (Colon syntax)
  interaction_dummy_pattern <- paste0("^", trt_var, "_")
  interaction_cols <- grep(interaction_dummy_pattern, names(model_data), value = TRUE)

  # 2. Create Control Data
  data_control <- model_data
  data_control[[trt_var]] <- factor(rep(ref_level, nrow(model_data)), levels = trt_levels)
  contrasts(data_control[[trt_var]]) <- trt_contrasts

  # For Colon syntax: Zero out dummies
  # For Pipe syntax: interaction_cols is empty, loop does nothing. Correct.
  for (col in interaction_cols) data_control[[col]] <- 0

  # 3. Create Treatment Data
  data_treatment <- model_data
  data_treatment[[trt_var]] <- factor(rep(alt_level, nrow(model_data)), levels = trt_levels)
  contrasts(data_treatment[[trt_var]]) <- trt_contrasts

  # For Colon syntax: Recreate dummies based on subgroup membership
  # For Pipe syntax: We do NOTHING. brms will see 'trt' = 1 and 'subgroup' = 'S1'
  # and automatically look up the random slope for S1.
  for (col in interaction_cols) {
    matched <- FALSE
    for (var_name in names(model_data)) {
      if (is.factor(model_data[[var_name]]) && var_name != trt_var) {
        var_levels <- levels(model_data[[var_name]])
        for (level in var_levels) {
          expected_dummy_name <- make.names(paste0(trt_var, "_", var_name, level), unique = FALSE)
          if (col == expected_dummy_name) {
            data_treatment[[col]] <- as.numeric(model_data[[var_name]] == level)
            matched <- TRUE
            break
          }
        }
        if (matched) break
      }
    }
    if (!matched) data_treatment[[col]] <- 0
  }

  return(list(control = data_control, treatment = data_treatment))
}

#' Get Posterior Predictions (Robust)
#'
#' UPDATED: Intelligently switches `re_formula` based on whether the treatment
#' effect is modeled as Fixed (Colon) or Random (Pipe).
#' @noRd
.get_posterior_predictions <- function(brms_fit, data_control, data_treatment, response_type, original_data, ndraws = NULL, trt_var) {

  checkmate::assert_class(brms_fit, "brmsfit")

  # --- Determine re_formula strategy ---
  # If trt_var is part of a random slope (Pipe model), we MUST include random effects (re_formula = NULL).
  # If trt_var is only fixed effects (Colon model), we SHOULD ignore random effects (re_formula = NA).

  has_random_trt <- FALSE
  re_structure <- tryCatch(
    brms::ranef(brms_fit, summary = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(re_structure)) {
    for (g_var in names(re_structure)) {
      re_coefs <- dimnames(re_structure[[g_var]])[[3]]
      # Robust check matching any coefficient name starting with trt_var
      if (any(grepl(paste0("^", trt_var), re_coefs))) {
        has_random_trt <- TRUE
        break
      }
    }
  }

  # re_formula = NULL includes all random effects
  # re_formula = NA includes NO random effects (Fixed only)
  prediction_re_formula <- if (has_random_trt) NULL else NA

  if (has_random_trt) {
    message("... detected Random Effects (Pipe model). Predicting with re_formula = NULL.")
  } else {
    message("... detected Fixed Effects (Colon model). Predicting with re_formula = NA.")
  }

  data_control$._sim_group <- "control"
  data_treatment$._sim_group <- "treatment"
  data_combined <- rbind(data_control, data_treatment)

  if (response_type == "survival") {
    message("... (reconstructing baseline hazard and getting linear predictors)...")

    linpred_combined <- if (is.null(ndraws)) {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula)
    } else {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula, ndraws = ndraws)
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
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula)
    } else {
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula, ndraws = ndraws)
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


