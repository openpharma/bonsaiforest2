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
#' @param brms_fit A `brmsfit` object. Fitted model object from `fit_brms_model()` or
#'   `run_brms_analysis()`. Must contain the necessary attributes for extracting treatment
#'   variable and response type information.
#' @param trt_var A character string or `NULL`. Treatment variable name. If `NULL`, automatically
#'   extracted from model attributes (set by `fit_brms_model()`). Must be a binary variable
#'   coded as 0/1 in the dataset.
#' @param data A data frame or `NULL`. Dataset used for model fitting. If `NULL`, automatically
#'   extracted from model attributes (set by `fit_brms_model()`). This dataset is used for
#'   generating counterfactual predictions.
#' @param subgroup_vars A character vector or `"auto"`. Subgroup variable names for which
#'   to estimate treatment effects. If `"auto"` (default), automatically detects treatment
#'   interaction terms (colon syntax) and random effect grouping factors (pipe syntax) from
#'   all formula components (`unshrunktermeffect`, `shprogeffect`, `shpredeffect`).
#' @param response_type A character string or `NULL`. Outcome type, one of `"binary"`,
#'   `"count"`, `"continuous"`, or `"survival"`. If `NULL`, automatically extracted from
#'   model attributes (set by `fit_brms_model()`). This determines the appropriate scale
#'   for effect estimation.
#' @param ndraws An integer or `NULL`. Number of posterior draws to use for estimation.
#'   If `NULL` (default), all available posterior draws are used. Reducing this can speed
#'   up computation at the cost of precision.
#' @param conf A numeric scalar. Credible interval confidence level. Default is 0.95
#'   (corresponding to 95% credible intervals). Example: use 0.90 for 90% credible intervals.
#'
#' @return `list` with two named elements:
#'   \describe{
#'     \item{`estimates`}{`tibble` where each row represents a subgroup (or "Overall" effect),
#'       with columns for `Subgroup`, `Median`, `CI_Lower`, and `CI_Upper` from posterior distribution}
#'     \item{`draws`}{`data.frame` containing posterior draws for each subgroup}
#'   }
#'
#' @importFrom checkmate assert_class assert_data_frame assert_string assert_choice
#' @importFrom checkmate assert_character assert_count assert_logical assert_subset
#' @importFrom stringr str_detect str_match_all
#' @importFrom dplyr bind_rows bind_cols
#' @importFrom stats contrasts contrasts<-
#' @importFrom tibble tibble
#' @importFrom utils head
#' @export
estimate_subgroup_effects <- function(brms_fit,
                                      trt_var = NULL,
                                      data = NULL,
                                      subgroup_vars = "auto",
                                      response_type = NULL,
                                      ndraws = NULL,
                                      conf = 0.95) {

  # --- 1. Validate inputs and determine which subgroups to analyze ---
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

  # Extract data from model attributes if not provided
  if (is.null(data)) {
    data <- attr(brms_fit, "model_data")
    if (is.null(data)) {
      # Fallback to brms_fit$data if model_data attribute is not available
      data <- brms_fit$data
      message("Using data from brms_fit$data (model_data attribute not found)")
    } else {
      message("Using data from model attributes")
    }
  }
  checkmate::assert_data_frame(data, min.rows = 1)

  # Check that trt_var exists in the data
  checkmate::assert_subset(trt_var, names(data))

  # Extract response_type from model attributes if not provided
  if (is.null(response_type)) {
    response_type <- attr(brms_fit, "response_type")
    if (is.null(response_type)) {
      stop("response_type must be specified or stored in the model attributes (via fit_brms_model()).")
    }
    message("Using response_type from model attributes: ", response_type)
  }
  # Check response type
  response_type <- checkmate::assert_choice(response_type,
                                            choices = c("continuous", "binary", "count", "survival"))
  checkmate::assert_count(ndraws, null.ok = TRUE, positive = TRUE)
  checkmate::assert_number(conf, lower = 0, upper = 1)

  # Use the data extracted from attributes or provided by user
  # This ensures consistency between predictions and subgroup membership
  model_data <- data

  # --- 2. Prepare Subgroup Variables ---
  message("Step 1: Identifying subgroups and creating counterfactuals...")

  prep <- .prepare_subgroup_vars(
    brms_fit = brms_fit,
    model_data = model_data,
    trt_var = trt_var,
    subgroup_vars = subgroup_vars
  )

  # --- 3. Create Counterfactuals ---
  # Create counterfactual datasets under potential outcomes framework
  # Control dataset: All observations assigned to control arm (treatment = 0)
  # Treatment dataset: All observations assigned to treatment arm (treatment = 1)
  # This enables estimation of average treatment effects via posterior predictive distributions
  counterfactual_data <- .create_counterfactual_datasets(
    model_data = prep$data,
    trt_var = trt_var
  )

  # --- 4. Generate posterior predictions ---
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
    display_names = prep$display_names,
    is_overall = prep$is_overall,
    response_type = response_type,
    conf = conf
  )

  message("Done.")
  return(results)
}

#' Strip _onehot Suffix from Variable Name
#'
#' Removes the _onehot suffix used internally for duplicate variables.
#' This is used to display clean variable names to users.
#' @noRd
.strip_onehot_suffix <- function(var_name) {
  sub("_onehot$", "", var_name)
}


#' Strip _onehot Suffix from Variable Name
#'
#' Removes the _onehot suffix used internally for duplicate variables.
#' This is used to display clean variable names to users.
#' @noRd
.strip_onehot_suffix <- function(var_name) {
  sub("_onehot$", "", var_name)
}


#' Prepare and Validate Subgroup Variables
#'
#' Detects treatment interaction terms from the fitted brms model.
#' Detects both fixed interactions (colon syntax) and random effects grouping
#' factors (pipe syntax) from all formula components.
#' @noRd
.prepare_subgroup_vars <- function(brms_fit, model_data, trt_var, subgroup_vars) {

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
    message(sprintf("Model data has %d rows and %d columns", nrow(model_data), ncol(model_data)))
    message(sprintf("Column names: %s", paste(names(model_data), collapse = ", ")))
    message(sprintf("Treatment variable: '%s'", trt_var))
    
    detected_vars <- character(0)

    # A. Detect FIXED interactions (Colon syntax)
    # -------------------------------------------
    # Detect treatment interactions from ANY formula component:
    # - unshrunktermeffect: all unshrunk terms (both main effects and interactions)
    # - shprogeffect: shrunk prognostic effects
    # - shpredeffect: shrunk predictive effects (treatment interactions)
    fixef_matrix <- brms::fixef(brms_fit)
    all_coefs <- rownames(fixef_matrix)
    
    message("All coefficient names:")
    message(paste(all_coefs, collapse = "\n"))
    
    # Look for any coefficient containing "trt:" pattern
    # This will find ALL treatment interactions regardless of formula component
    interaction_pattern <- paste0(trt_var, ":")
    message(sprintf("Looking for treatment interactions with pattern: '%s'", interaction_pattern))
    interaction_coefs <- grep(interaction_pattern, all_coefs, value = TRUE, fixed = TRUE)
    
    # Also check for reverse order: "var:trt" 
    reverse_pattern <- paste0(":", trt_var)
    reverse_coefs <- grep(reverse_pattern, all_coefs, value = TRUE, fixed = TRUE)
    interaction_coefs <- unique(c(interaction_coefs, reverse_coefs))
    
    message(sprintf("Found %d treatment interaction coefficients", length(interaction_coefs)))
    if (length(interaction_coefs) > 0) {
      message("Treatment interaction coefficients found:")
      message(paste(interaction_coefs, collapse = "\n"))
    }
    
    if (length(interaction_coefs) > 0) {
      # Extract the variable names from the interaction terms
      # Patterns to handle (all formula components):
      # - "unshrunktermeffect_trt:biomarkerLow" or "unshrunktermeffect_biomarker:trt"
      # - "shpredeffect_trt:sexF" or "shpredeffect_sex:trt"
      # - "shprogeffect_trt:regionEast" (less common but possible)
      # General form: "componentname_trt:var" or "componentname_var:trt"
      
      for (coef in interaction_coefs) {
        # Remove the effect prefix (e.g., "unprogeffect_", "shpredeffect_")
        # Pattern: effectname_trt:var or effectname_var:trt
        coef_without_prefix <- sub("^[^_]+_", "", coef)
        
        # Now we have either "trt:varnameLevel" or "varnameLevel:trt"
        # Split by colon to get the two parts
        parts <- strsplit(coef_without_prefix, ":", fixed = TRUE)[[1]]
        
        if (length(parts) == 2) {
          # Determine which part is NOT trt
          if (parts[1] == trt_var) {
            varname_with_level <- parts[2]
          } else if (parts[2] == trt_var) {
            varname_with_level <- parts[1]
          } else {
            # Neither part is trt, skip
            next
          }
          
          # Check which factor variable in the data this corresponds to
          # by seeing if any factor variable name is a prefix of varname_with_level
          # Sort by length descending to match longer names first (e.g., x_10 before x_1)
          factor_vars <- names(model_data)[sapply(model_data, is.factor)]
          factor_vars <- factor_vars[factor_vars != trt_var]
          factor_vars_sorted <- factor_vars[order(nchar(factor_vars), decreasing = TRUE)]
          
          for (var_name in factor_vars_sorted) {
            # Check if varname_with_level starts with var_name
            # e.g., "biomarkerLow" starts with "biomarker", "x_10a" starts with "x_10"
            if (startsWith(varname_with_level, var_name)) {
              detected_vars <- c(detected_vars, var_name)
              message(sprintf("Detected subgroup variable '%s' from coefficient '%s'", var_name, coef))
              break
            }
          }
        }
      }
    }
    
    # Also check for legacy dummy column approach (backward compatibility)
    # This was used in older versions before the consolidated architecture
    dummy_pattern <- paste0("^", trt_var, "_")
    interaction_cols <- grep(dummy_pattern, names(model_data), value = TRUE)
    
    if (length(interaction_cols) > 0) {
      for (var_name in names(model_data)) {
        if (is.factor(model_data[[var_name]]) && var_name != trt_var) {
          # Check if dummies match this variable
          var_levels <- levels(model_data[[var_name]])
          expected_dummies <- make.names(paste0(trt_var, "_", var_name, var_levels), unique = FALSE)
          if (any(expected_dummies %in% interaction_cols)) {
            detected_vars <- c(detected_vars, var_name)
          }
        }
      }
    }

    # B. Detect RANDOM effects grouping factors (Pipe syntax)
    # -----------------------------------------------------
    # Check Stan parameter names for random effects like: r_biomarker__shpredeffect[High,trt]
    # This pattern indicates treatment slopes varying by biomarker levels.
    # Random effects are typically used for shrunk predictive effects (treatment interactions)
    # specified with pipe syntax: ~ (0 + trt || biomarker)
    message("Checking for random effects parameters...")
    all_params <- tryCatch(
      brms::variables(brms_fit),
      error = function(e) {
        message("Error retrieving variables: ", e$message)
        character(0)
      }
    )
    
    message(sprintf("Retrieved %d total parameters from model", length(all_params)))
    
    if (length(all_params) > 0) {
      # Look for random effects parameters with treatment
      # Pattern: r_GROUPVAR__TERM[LEVEL,TRT]
      # Example: r_biomarker__shpredeffect[High,trt] or r_age_group__shpredeffect[<50,trt]
      # Need to handle special characters in level names AND underscores in variable names
      # Use __ (double underscore) as the delimiter between grouping var and term
      re_pattern <- sprintf("^r_(.+)__[^\\[]+\\[[^,]+,%s\\]", trt_var)
      message(sprintf("Using regex pattern: '%s'", re_pattern))
      re_params <- grep(re_pattern, all_params, value = TRUE)
      
      message(sprintf("Found %d matching random effect parameters", length(re_params)))
      
      if (length(re_params) > 0) {
        message(sprintf("Found %d random effect parameters with treatment slopes:", length(re_params)))
        message(paste("  ", head(re_params, 3), collapse = "\n"))
        if (length(re_params) > 3) message(sprintf("  ... and %d more", length(re_params) - 3))
        
        # Extract grouping variable names
        for (param in re_params) {
          # Extract the grouping variable name from r_GROUPVAR__...
          # Use the double underscore __ as delimiter
          grouping_var <- sub("^r_(.+)__.*", "\\1", param)
          message(sprintf("Extracted grouping variable: '%s' from parameter: '%s'", grouping_var, param))
          # Verify this variable exists in the model data
          if (grouping_var %in% names(model_data)) {
            if (is.factor(model_data[[grouping_var]]) || is.character(model_data[[grouping_var]])) {
              detected_vars <- c(detected_vars, grouping_var)
              message(sprintf("Detected subgroup variable '%s' from random effects (treatment slopes)", grouping_var))
            } else {
              message(sprintf("Variable '%s' exists but is not a factor/character", grouping_var))
            }
          } else {
            message(sprintf("Variable '%s' not found in model data", grouping_var))
          }
        }
      } else {
        message("No random effect parameters matching the pattern were found")
      }
    } else {
      message("Could not retrieve model variables")
    }

    detected_vars <- unique(detected_vars)

    if (length(detected_vars) == 0) {
      message("...no treatment interaction terms detected. Calculating overall treatment effect.")
      message("   (averaging over all covariates regardless of formula component)")
      is_overall <- TRUE
      subgroup_vars <- "Overall"
    } else {
      subgroup_vars <- detected_vars
      message(paste("...detected subgroup variable(s):", paste(subgroup_vars, collapse = ", ")))
    }
  } else {
    checkmate::assert_subset(subgroup_vars, names(model_data))
  }

  if (is_overall) {
    model_data$Overall <- "Overall"
  }

  # Create mapping of internal variable names to display names
  # Strip _onehot suffix for user-facing output
  display_names <- setNames(
    sapply(subgroup_vars, .strip_onehot_suffix),
    subgroup_vars
  )
  
  return(list(
    subgroup_vars = subgroup_vars,
    display_names = display_names,
    data = model_data,
    is_overall = is_overall
  ))
}

#' Create Counterfactual Datasets
#'
#' Creates two versions of the data: one with all observations set to control,
#' and one with all set to treatment. Handles both numeric (0/1) and factor
#' treatment variables. Works for treatment effects in any formula component.
#' @noRd
.create_counterfactual_datasets <- function(model_data, trt_var) {
  checkmate::assert_data_frame(model_data)

  # Check if treatment is numeric or factor
  is_numeric_trt <- is.numeric(model_data[[trt_var]])
  
  if (is_numeric_trt) {
    # Treatment is numeric (0/1)
    # Validate that treatment is binary 0/1
    trt_vals <- unique(model_data[[trt_var]])
    if (!all(trt_vals %in% c(0, 1))) {
      stop("Treatment variable '", trt_var, "' must be binary (0/1) when numeric.")
    }
    
    control_val <- 0
    treatment_val <- 1
    
  } else {
    # Treatment is factor (for backward compatibility with older workflows)
    if (!is.factor(model_data[[trt_var]])) model_data[[trt_var]] <- as.factor(model_data[[trt_var]])
    trt_levels <- levels(model_data[[trt_var]])
    trt_contrasts <- contrasts(model_data[[trt_var]])
    
    ref_level <- trt_levels[1]
    alt_level <- trt_levels[2]
  }

  # 2. Create Control Data
  data_control <- model_data
  
  # 3. Create Treatment Data
  data_treatment <- model_data
  
  if (is_numeric_trt) {
    # For numeric treatment, interactions are handled automatically by R/brms
    # When we have trt:sex where trt is 0/1 and sex is a factor with contrasts,
    # R automatically computes trt * sexM, trt * sexF, etc.
    # We just need to set trt to 0 or 1, and keep all other variables as-is
    data_control[[trt_var]] <- rep(control_val, nrow(model_data))
    data_treatment[[trt_var]] <- rep(treatment_val, nrow(model_data))
    
  } else {
    # Factor treatment with explicit interaction dummy columns
    data_control[[trt_var]] <- factor(rep(ref_level, nrow(model_data)), levels = trt_levels)
    contrasts(data_control[[trt_var]]) <- trt_contrasts
    
    data_treatment[[trt_var]] <- factor(rep(alt_level, nrow(model_data)), levels = trt_levels)
    contrasts(data_treatment[[trt_var]]) <- trt_contrasts
    
    # 1. Identify Fixed Interaction Dummies (Colon syntax with factor treatment)
    interaction_dummy_pattern <- paste0("^", trt_var, "_")
    interaction_cols <- grep(interaction_dummy_pattern, names(model_data), value = TRUE)
    
    # For Colon syntax: Zero out dummies in control
    for (col in interaction_cols) data_control[[col]] <- 0
    
    # For Colon syntax: Recreate dummies based on subgroup membership in treatment
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
  }

  return(list(control = data_control, treatment = data_treatment))
}

#' Get Posterior Predictions
#'
#' Intelligently switches `re_formula` based on whether the treatment
#' effect is modeled as fixed effects (colon syntax) or random effects (pipe syntax).
#' @noRd
.get_posterior_predictions <- function(brms_fit, data_control, data_treatment, response_type, original_data, ndraws = NULL, trt_var) {

  checkmate::assert_class(brms_fit, "brmsfit")

  # Determine prediction strategy based on model structure:
  # - Fixed effects model (colon syntax): re_formula = NA (population-level predictions)
  # - Random effects model (pipe syntax): re_formula = NULL (include group-level effects)
  # This distinction ensures appropriate marginalization across model architectures
  has_random_trt <- FALSE
  
  # Check Stan parameter names for random effects with treatment
  # Pattern: r_GROUPVAR__TERM[LEVEL,TRT]
  all_params <- tryCatch(
    brms::variables(brms_fit),
    error = function(e) character(0)
  )
  
  if (length(all_params) > 0) {
    # Look for random effects parameters with treatment slopes
    # Pattern: r_GROUPVAR__TERM[LEVEL,TRT] where GROUPVAR can contain underscores
    # Use .+ instead of [^_]+ to match variable names with underscores (like age_group)
    re_pattern <- sprintf("^r_.+__[^\\[]+\\[[^,]+,%s\\]", trt_var)
    re_params <- grep(re_pattern, all_params, value = TRUE)
    
    if (length(re_params) > 0) {
      has_random_trt <- TRUE
      message(sprintf("... found %d random treatment slope parameters", length(re_params)))
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
    # Survival model estimation requires two-stage process:
    # 1. Extract linear predictors from fitted model
    # 2. Reconstruct baseline hazard via B-spline basis (computed during model fitting)
    # This approach leverages brms' internal spline representation for computational efficiency
    message("... (reconstructing baseline hazard and getting linear predictors)...")

    linpred_combined <- if (is.null(ndraws)) {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula,
                             allow_new_levels = FALSE)
    } else {
      brms::posterior_linpred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula, 
                             ndraws = ndraws, allow_new_levels = FALSE)
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
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula, 
                           allow_new_levels = FALSE)
    } else {
      brms::posterior_epred(brms_fit, newdata = data_combined, re_formula = prediction_re_formula, 
                           ndraws = ndraws, allow_new_levels = FALSE)
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

#' Extract Baseline Hazard Details for Survival Models
#'
#' Parses the brms formula to extract baseline hazard function parameters
#' and reconstructs the hazard using B-spline basis functions.
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


#' Calculate and Summarize Marginal Effects
#' @noRd
.calculate_and_summarize_effects <- function(posterior_preds, original_data, subgroup_vars, display_names, is_overall, response_type, conf = 0.95) {

  all_results_list <- list()
  all_draws_list <- list()

  # Pre-compute factor info for iteration
  subgroup_factors <- list()
  if (!is_overall) {
    for (var in subgroup_vars) {
      factor_var <- as.factor(original_data[[var]])
      # Use the actual levels from the data - brms preserves special characters
      # in factor levels even if they contain <, >, -, etc.
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

      # Use display name (without _onehot suffix) for user-facing output
      display_var_name <- if (is_overall) "Overall" else display_names[current_subgroup_var]
      subgroup_name <- if (is_overall) "Overall" else paste0(display_var_name, ": ", level)
      all_draws_list[[subgroup_name]] <- effect_draws

      point_estimate <- median(effect_draws, na.rm = TRUE)
      alpha <- 1 - conf
      ci <- quantile(effect_draws, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

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
#' Computes the average hazard ratio across all time points by integrating
#' the marginal survival curves for control and treatment arms.
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
#'
#' Computes marginal survival curves by averaging over individual-level
#' survival predictions for a given subgroup. Uses vectorized operations
#' to avoid loops and improve performance.
#' @noRd
.get_marginal_survival_vectorized <- function(linpred_posterior, H0_post_list, sub_indices, strat_variable, full_data) {
  subgroup_linpred <- linpred_posterior[, sub_indices, drop = FALSE]
  n_draws <- nrow(subgroup_linpred)
  n_individuals <- ncol(subgroup_linpred)
  
  # Pre-compute exp(eta) once for all draws
  exp_linpred <- exp(subgroup_linpred)

  if (is.null(strat_variable)) {
    H0_post <- H0_post_list[["_default_"]]
    n_times <- nrow(H0_post)
    
    # Vectorized approach: compute all draws at once using lapply
    # Then rbind the results
    S_list <- lapply(seq_len(n_draws), function(i) {
      # Use tcrossprod instead of outer (slightly faster)
      # tcrossprod(a, b) = a %*% t(b) which is same as outer(a, b) for our case
      S_indiv <- exp(-tcrossprod(H0_post[, i], exp_linpred[i, ]))
      rowMeans(S_indiv)
    })
    
    S_marginal <- do.call(rbind, S_list)
    return(S_marginal)
  }

  # Stratified logic (vectorized across draws)
  subgroup_strata <- full_data[[strat_variable]][sub_indices]
  unique_strata <- unique(subgroup_strata)
  n_times <- nrow(H0_post_list[[unique_strata[1]]])
  
  # Vectorized approach for stratified case
  S_list <- lapply(seq_len(n_draws), function(i) {
    S_indiv <- matrix(NA, nrow = n_times, ncol = n_individuals)
    for (s in unique_strata) {
      idx <- which(subgroup_strata == s)
      h0 <- H0_post_list[[s]][, i]
      eta_s <- exp_linpred[i, idx]
      S_indiv[, idx] <- exp(-tcrossprod(h0, eta_s))
    }
    rowMeans(S_indiv, na.rm = TRUE)
  })
  
  S_marginal <- do.call(rbind, S_list)
  return(S_marginal)
}

#' Calculate AHR from Survival Curves
#'
#' Computes the average hazard ratio using the restricted mean survival time approach.
#' @noRd
.calculate_ahr_vectorized <- function(S_control, S_treatment) {
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


