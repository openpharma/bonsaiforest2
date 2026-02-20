#' Simulate Study Data for Multiple Endpoints and Scenarios
#'
#' Generates replicated datasets for TTE, binary, count, or continuous
#' endpoints based on 6 scenarios from Wolbers et al. (2025).
#'
#' @param endpoint The type of outcome data to generate.
#'   One of: "tte", "binary", "count", "continuous".
#' @param scenario The simulation scenario (1-6).
#' @param n_datasets The number of datasets (replications) to generate.
#' @param ... Additional arguments. For "tte", these are passed to
#'   `simul_scenario` (e.g., n_patients, n_events).
#'
#' @return A list of `n_datasets` data frames, each containing a
#'   simulated dataset.
#'
#' @details
#'   - Assumes `simul_scenario`, `simul_data`, `simul_covariates`, etc.,
#'     are all available in the environment.
#'
simul_study_data <- function(endpoint = c("tte", "binary", "count", "continuous"),
                             scenario = c("1", "2", "3", "4", "5", "6"),
                             n_datasets = 1000,
                             ...) {
  # Match arguments
  endpoint <- match.arg(endpoint)
  scenario <- match.arg(scenario)
  assert_count(n_datasets)

  # --- Branch 1: Time-to-Event (TTE) Endpoint ---
  # Use the user's provided function for TTE.
  if (endpoint == "tte") {
    message("Endpoint 'tte': Calling provided simul_scenario() function...")

    # Call the user's function (must be in the Global Environment)
    # We pass '...' which includes n_patients, n_events, etc.
    return(
      simul_scenario(
        scenario = scenario,
        n_datasets = n_datasets,
        ...
      )
    )
  }

  # --- Branch 2: New Endpoints (Binary, Count, Continuous) ---
  message(paste("Endpoint '", endpoint, "': Calling new simulation generator...", sep = ""))

  # Get model parameters (N, intercept, etc.) for the chosen endpoint
  model_params <- .get_model_parameters(endpoint)

  # Get the (non-zero) coefficients for this scenario
  scenario_coefs <- .get_scenario_coefs(scenario, model_params)

  # Generate the list of n_datasets
  replicate(
    n_datasets,
    .simul_new_endpoint_single(
      n = model_params$N,
      endpoint = endpoint,
      model_params = model_params,
      coefs = scenario_coefs,
      add_interaction = (scenario == "6")
    ),
    simplify = FALSE
  )
}

#' Simulate a Single Dataset for New Endpoints
#' (Internal Helper Function)
#'
#' @description
#'   This function is designed to be a "twin" of your `simul_data` function.
#'   It uses the *exact same* covariate generation (`simul_covariates`)
#'   and *exact same* design matrix logic, but substitutes a new
#'   outcome (y) for the TTE outcome.
#'
#' @param n Sample size.
#' @param endpoint "binary", "count", or "continuous".
#' @param model_params List of parameters from .get_model_parameters().
#' @param coefs Named vector of *non-zero* coefficients from .get_scenario_coefs().
#' @param add_interaction Boolean, for scenario 6.
#'
#' @return A single simulated data.frame.
.simul_new_endpoint_single <- function(n,
                                       endpoint,
                                       model_params,
                                       coefs,
                                       add_interaction = FALSE) {

  # 1. GENERATE COVARIATES
  # Call your 'simul_covariates' function, just like 'simul_data' does.
  covariates <- simul_covariates(n = n, p_catvar = 10, add_contvars = FALSE, arm_factor = TRUE)

  # 2. BUILD DESIGN MATRIX
  # This logic is copied *directly* from your 'simul_data' function
  # to ensure the design matrix is identical.
  subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10

  if (add_interaction) {
    covariates$x_1_2 <- factor(
      with(
        covariates,
        paste(as.character(x_1), as.character(x_2), sep = "")
      )
    )
    subgroup_model <- stats::update(subgroup_model, ~ . + x_1_2)
  }

  design_main <- stats::model.matrix(
    stats::update(subgroup_model, ~ arm + .),
    data = covariates
  )
  subgroup_vars <- all.vars(subgroup_model)
  design_ia <- NULL
  for (j in subgroup_vars) {
    ia_j <- stats::model.matrix(
      stats::as.formula(paste("~", j, "-1")),
      data = covariates
    ) * design_main[, "arm1"] # Assumes arm_factor=TRUE created 'arm1'
    design_ia <- cbind(design_ia, ia_j)
  }
  colnames(design_ia) <- paste(colnames(design_ia), "arm", sep = "_")
  colnames(design_ia) <- gsub(" ", "", colnames(design_ia))
  design_matrix <- cbind(design_main, design_ia)

  if (add_interaction) {
    # Clean up covariates data.frame, just as in simul_data
    covariates$x_1_2 <- NULL
  }

  # 3. BUILD COEFFICIENT VECTOR
  # This logic is also from 'simul_data': create a full vector
  # of zeros and fill in the non-zero scenario coefficients.
  reg_coef <- rep(0, ncol(design_matrix))
  names(reg_coef) <- colnames(design_matrix)

  # Check that our scenario coefs are all valid
  assert_subset(names(coefs), names(reg_coef))
  reg_coef[names(coefs)] <- coefs

  # Validate that reg_coef is numeric (no NAs)
  if (any(is.na(reg_coef))) {
    stop("NA values detected in reg_coef. NA positions: ",
         paste(names(reg_coef)[is.na(reg_coef)], collapse = ", "),
         "\nCoefficients provided: ", paste(names(coefs), "=", coefs, collapse = ", "))
  }
  if (!is.numeric(reg_coef)) {
    stop("reg_coef is not numeric. Class: ", class(reg_coef))
  }

  # 4. CALCULATE LINEAR PREDICTOR (eta)
  lp <- design_matrix %*% reg_coef

  # 5. GENERATE NEW OUTCOME 'y'
  # This is the only part that's truly "new".
  y <- switch(
    endpoint,
    "binary" = {
      p <- plogis(lp) # Inverse-logit to get probabilities
      rbinom(n, 1, p)
    },
    "count" = {
      mu <- exp(lp) # Inverse-log to get rates
      rnbinom(n, mu = mu, size = model_params$overdispersion)
    },
    "continuous" = {
      rnorm(n, mean = lp, sd = model_params$sd)
    }
  )

  # 6. RETURN FINAL DATA.FRAME
  # We return 'covariates' (which has x_1...x_10, arm) and the new outcome 'y'.
  # This matches the structure of the TTE data (which has covariates, tt_pfs, ev_pfs).
  data.frame(id = 1:n, covariates, y)
}

#' Get Model Parameters for New Endpoints (Internal Helper Function)
#'
#' Returns a list of parameters (N, intercept, base_effect, sigma_scale)
#' for the binary, count, and continuous models, calibrated for ~80% power.
#'
#' CALIBRATION TARGETS:
#' - Binary: N = 800 (400 per arm) for 80% power to detect OR=0.66
#' - Count: N = 770 (385 per arm) with overdispersion shape=0.34
#' - Continuous: N = 770 (385 per arm) with residual SD=2.06
#'
#' @param endpoint One of "binary", "count", "continuous".
#' @return A list of model parameters including N, intercept, base_effect,
#'   sigma_scale, and endpoint-specific parameters.
.get_model_parameters <- function(endpoint) {
  switch(
    endpoint,
    "binary" = list(
      endpoint = "binary",
      # Sample size: 800 (400 per arm) for ~80% power to detect OR=0.66
      N = 800,
      # Intercept -0.2 gives baseline event rate ~45%
      intercept = -0.2,
      # Base treatment effect: log(0.66) = -0.416 (beneficial)
      base_effect = log(0.66),
      # No coefficient scaling for binary
      sigma_scale = 1.0
    ),
    "count" = list(
      endpoint = "count",
      # Sample size: 770 (385 per arm) for ~80% power with overdispersion
      N = 770,
      # Intercept 0.0 gives baseline count = exp(0.0) = 1.0
      intercept = 0.0,
      # Base treatment effect: log(0.66) = -0.416 (beneficial: reduced rate)
      base_effect = log(0.66),
      # Negative Binomial shape parameter: 0.34 (variance/mean ~3)
      overdispersion = 0.34,
      # No coefficient scaling for count
      sigma_scale = 1.0
    ),
    "continuous" = list(
      endpoint = "continuous",
      # Sample size: 770 (385 per arm) for ~80% power with residual SD=2.06
      N = 770,
      # Intercept 0.0 (baseline mean = 0)
      intercept = 0.0,
      # Base treatment effect: d = -0.3 (beneficial: lower score)
      base_effect = -0.3,
      # Residual standard deviation: 2.06 (calibrated for power)
      sd = 2.06,
      # No coefficient scaling for continuous
      sigma_scale = 1.0
    )
  )
}

#' Get Scenario Coefficients (Direct Translation of TTE Logic)
#'
#' @param scenario The simulation scenario (1-6).
#' @param model_params A list of parameters from .get_model_parameters().
#' @return A named vector of coefficients.
#' Internal: Get Scenario Coefficients
#' Applies endpoint-specific scaling factors and treatment effect signs.
.get_scenario_coefs <- function(scenario, model_params) {
  RNGkind('Mersenne-Twister')
  set.seed(5)
  # --- SCALING FACTOR ---
  # Extract from model_params (0.85 for TTE, 1.0 for others)
  sigma_scale <- model_params$sigma_scale

  subgroup_covariates <- c(
    "x_1a", "x_1b", "x_2a", "x_2b", "x_3a", "x_3b",
    "x_4a", "x_4b", "x_4c", "x_5a", "x_5b", "x_5c", "x_5d",
    "x_6a", "x_6b", "x_7a", "x_7b", "x_8a", "x_8b", "x_8c",
    "x_9a", "x_9b", "x_10a", "x_10b", "x_10c"
  )
  ia_names <- paste0(subgroup_covariates, "_arm")

  # Initialize
  coefs <- c("(Intercept)" = as.numeric(model_params$intercept))

  # Determine sign for treatment effects based on endpoint
  # TTE: -log(HR) because lower hazard is better (HR<1 is beneficial)
  # Binary: log(OR) because interpretation depends on outcome direction
  # Count: log(RR) because Y is "bad" events, so RR<1 is beneficial
  # Continuous: direct coefficient because Y is "bad" score, so negative is beneficial
  endpoint <- model_params$endpoint

  # For prognostic effects, always use same direction as TTE AFT model
  # In TTE AFT: positive coef = longer survival = beneficial
  # x_4c: HR=0.7 (beneficial) -> AFT coef = -log(0.7) * 0.85 = positive
  # x_6b: HR=1.5 (harmful) -> AFT coef = -log(1.5) * 0.85 = negative
  # For non-TTE: keep same direction of effect
  if (endpoint == "continuous") {
    # For continuous, use direct scaling without log transformation
    coefs["x_4c"] <- -log(0.7) * sigma_scale  # Beneficial (lower bad score)
    coefs["x_6b"] <- -log(1.5) * sigma_scale  # Harmful (higher bad score)
  } else {
    # For binary/count, use log scale
    coefs["x_4c"] <- log(0.7) * sigma_scale   # Beneficial
    coefs["x_6b"] <- log(1.5) * sigma_scale   # Harmful
  }

  # Helper to get treatment effect coefficient based on endpoint
  # Returns the coefficient corresponding to a given hazard ratio
  get_trt_coef <- function(hr) {
    if (endpoint == "continuous") {
      # For continuous: use direct standardized effect
      # HR=0.66 maps to d=-0.3 (beneficial = negative for "bad" outcome)
      # Scale proportionally: log(hr) / log(0.66) * (-0.3)
      (log(hr) / log(0.66)) * (-0.3) * sigma_scale
    } else {
      # For binary/count: use log(effect measure)
      # OR=0.66 or RR=0.66 (beneficial = <1)
      log(hr) * sigma_scale
    }
  }

  if (scenario == "1") {
    coefs["arm1"] <- get_trt_coef(0.66)
  }
  else if (scenario == "2") {
    coefs["arm1"] <- get_trt_coef(0.66)
    # x_4a has NO effect. Cancels out the arm effect.
    coefs["x_4a_arm"] <- -get_trt_coef(0.66)
    coefs["x_4b_arm"] <- get_trt_coef(0.8)
    coefs["x_4c_arm"] <- get_trt_coef(0.8)
  }
  else if (scenario == "3") {
    coefs["arm1"] <- 0
    coefs["x_4a_arm"] <- get_trt_coef(0.5)
    coefs["x_4b_arm"] <- get_trt_coef(1.25)
    coefs["x_4c_arm"] <- get_trt_coef(1.25)
  }
  else if (scenario == "4") {
    coefs["arm1"] <- 0
    # Generate random treatment effects with mild heterogeneity
    if (endpoint == "continuous") {
      # For continuous: random effects around 0
      ia_coefs <- rnorm(25, mean = 0, sd = 0.15) * sigma_scale
    } else {
      # For binary/count: random log effects
      # Negative mean for beneficial direction on average
      ia_coefs <- -rnorm(25, mean = 0, sd = 0.15) * sigma_scale
    }
    names(ia_coefs) <- ia_names
    coefs <- c(coefs, ia_coefs)
  }
  else if (scenario == "5") {
    coefs["arm1"] <- 0
    # Generate random treatment effects with large heterogeneity
    if (endpoint == "continuous") {
      # For continuous: random effects around 0
      ia_coefs <- rnorm(25, mean = 0, sd = 0.3) * sigma_scale
    } else {
      # For binary/count: random log effects
      ia_coefs <- -rnorm(25, mean = 0, sd = 0.3) * sigma_scale
    }
    names(ia_coefs) <- ia_names
    coefs <- c(coefs, ia_coefs)
  }
  else if (scenario == "6") {
    coefs["arm1"] <- get_trt_coef(0.66)
    coefs["x_1_2aa_arm"] <- get_trt_coef(1.5)
    coefs["x_1_2ba_arm"] <- get_trt_coef(0.5)
    coefs["x_1_2ab_arm"] <- get_trt_coef(0.92)
    coefs["x_1_2bb_arm"] <- get_trt_coef(1.07)
  }

  return(coefs)
}

## Helper functions from original Bonsaiforest library

#' Helper for Cutting into Normal Quantiles
#'
#' @param x (`numeric`)\cr continuous covariate values.
#' @param prob (`numeric`)\cr probabilities for the standard normal quantiles.
#' @param labels (`character`)\cr corresponding labels for the resulting factor.
#'
#' @return The factor variable.
#' @keywords internal
cut_norm_quant <- function(x, prob, labels = letters[seq_along(c(prob, 1))]) {
  assert_numeric(x)
  assert_numeric(prob)
  assert_character(labels)
  assert_true(identical(length(prob) + 1L, length(labels)))

  norm_quantiles <- stats::qnorm(p = prob)
  breaks <- c(-Inf, norm_quantiles, Inf)
  cut(x, breaks, labels = labels)
}

#' Generation of a Design Matrix for Simulations
#'
#' This function uses a block diagonal covariance matrix for the underlying
#' multivariate normal data to create the design matrix in blocks of 10, see
#' the details.
#'
#' @details
#' The following pattern is repeated for the covariate blocks:
#'
#' - The first 5 covariates are uncorrelated with everything.
#' - The covariates 6 to 8 have "moderate" correlation (0.25) between each other.
#' - The covariates 9 and 10 have "high" correlation (0.5).
#'
#' By default, only the resulting categorical covariates obtained by thresholding
#' are included. Optionally also the original continuous covariates are included
#' in the returned design matrix.
#'
#' @param n (`count`)\cr number of rows (observations).
#' @param p_catvar (`count`)\cr number of covariates (excluding treatment arm).
#' @param add_contvars (`flag`)\cr whether to add continuous covariates.
#' @param arm_factor (`flag`)\cr whether to make the arm variable a factor.
#'
#' @return The design matrix.
#' @export
#'
#' @examples
#' simul_covariates(n = 10, p_catvar = 3, add_contvars = FALSE)
#' simul_covariates(n = 10, p_catvar = 3, add_contvars = TRUE)
#' simul_covariates(n = 10, p_catvar = 3, add_contvars = TRUE, arm_factor = TRUE)
simul_covariates <- function(n, p_catvar = 10, add_contvars = FALSE, arm_factor = FALSE) {
  assert_count(n, positive = TRUE)
  assert_count(p_catvar, positive = TRUE)
  assert_flag(add_contvars)
  assert_flag(arm_factor)

  sigma <- matrix(0, nrow = 10, ncol = 10)
  first_grp <- 1:5
  second_grp <- 6:8
  third_grp <- 9:10
  sigma[first_grp, first_grp] <- 0
  sigma[second_grp, second_grp] <- 0.25
  sigma[third_grp, third_grp] <- 0.5
  diag(sigma) <- 1

  no_10_blocks <- ceiling(p_catvar / 10)
  x <- NULL
  z <- NULL
  for (j in seq_len(no_10_blocks)) {
    # Continuous version.
    z_j <- data.frame(MASS::mvrnorm(n, mu = rep(0, 10), Sigma = sigma))
    colnames(z_j) <- paste("z", (j - 1) * 10 + 1:10, sep = "_")
    z <- if (j == 1) {
      z_j
    } else {
      cbind(z, z_j)
    }
    # Categorized version.
    x_j <- data.frame(
      v1 = cut_norm_quant(z_j[, 1], prob = 0.5),
      v2 = cut_norm_quant(z_j[, 2], prob = 0.4),
      v3 = cut_norm_quant(z_j[, 3], prob = 0.2),
      v4 = cut_norm_quant(z_j[, 4], prob = c(0.3, 0.6)),
      v5 = cut_norm_quant(z_j[, 5], prob = c(0.15, 0.3, 0.6)),
      v6 = cut_norm_quant(z_j[, 6], prob = 0.4),
      v7 = cut_norm_quant(z_j[, 7], prob = 0.4),
      v8 = cut_norm_quant(z_j[, 8], prob = c(0.2, 0.5)),
      v9 = cut_norm_quant(z_j[, 9], prob = 0.2),
      v10 = cut_norm_quant(z_j[, 10], prob = c(0.2, 0.5))
    )
    colnames(x_j) <- paste("x", (j - 1) * 10 + seq_len(10), sep = "_")
    x <- if (j == 1) {
      x_j
    } else {
      cbind(x, x_j)
    }
  }
  n_ctrl <- n %/% 2
  n_exp <- n - n_ctrl
  trt_arm <- sample(
    rep(
      c(0, 1),
      c(n_ctrl, n_exp)
    )
  )
  index_catvar <- seq_len(p_catvar)
  x <- cbind(
    arm = if (arm_factor) factor(trt_arm) else trt_arm,
    x[, index_catvar, drop = FALSE]
  )
  if (add_contvars) {
    x <- cbind(x, z[, index_catvar, drop = FALSE])
  }
  x
}

#' Simulation of Progression Free Survival Times
#'
#' @param lp_aft (`numeric`)\cr linear predictor values for the accelerate failure time model (AFT).
#' @param sigma_aft (`number`)\cr standard deviation for the AFT model.
#' @param recr_duration (`number`)\cr duration of recruitment.
#' @param rate_cens (`number`)\cr rate for the exponentially distributed censoring process.
#' @param n_events (`count`)\cr number of events to reach for the study end.
#' @param add_uncensored_pfs (`flag`)\cr whether to add the uncensored PFS as well to the resulting
#'   `data.frame`.
#'
#' @return A `data.frame` with columns `tt_pfs` (PFS time) and `ev_pfs` (corresponding
#'   event indicator with 1 for an event and 0 for censored), and optionally
#'   `tt_pfs_uncens`.
#' @export
#'
#' @examples
#' set.seed(123)
#' simul_pfs(
#'   lp_aft = rnorm(100),
#'   sigma_aft = 1,
#'   recr_duration = 0.2,
#'   rate_cens = 2,
#'   n_events = 20
#' )
simul_pfs <- function(lp_aft,
                      sigma_aft,
                      recr_duration,
                      rate_cens,
                      n_events,
                      add_uncensored_pfs = FALSE) {
  assert_numeric(lp_aft)
  assert_number(sigma_aft, lower = .Machine$double.xmin)
  assert_number(recr_duration, lower = .Machine$double.xmin)
  assert_number(rate_cens, lower = .Machine$double.xmin)
  assert_count(n_events, positive = TRUE)
  assert_flag(add_uncensored_pfs)

  n <- length(lp_aft)
  # Uncensored event time.
  log_tt_pfs <- c(lp_aft + sigma_aft * log(stats::rexp(n, rate = 1)))
  tt_pfs_uncens <- exp(log_tt_pfs)

  # Censoring step 1:
  # with rate_cens.
  tt_pfs_cens1 <- stats::rexp(n, rate = rate_cens)
  tt_pfs_cens1 <- pmin(tt_pfs_uncens, tt_pfs_cens1)
  ev_pfs_cens1 <- ifelse(tt_pfs_uncens <= tt_pfs_cens1, 1, 0)
  if (sum(ev_pfs_cens1) < n_events) {
    stop(paste(
      "Impossible to reach", n_events,
      "events with", n, "patients,",
      "a censoring rate of", rate_cens,
      "and the specified linear predictor."
    ))
  }

  # Censoring step 2:
  # due to staggerred recruitment and recruiting only until target_ev
  # events have been observed.
  rec_time <- stats::runif(n, min = 0, max = recr_duration)
  tt_pfs_cens1_calendar <- rec_time + tt_pfs_cens1
  study_stop_time <- sort(tt_pfs_cens1_calendar[ev_pfs_cens1 == 1])[n_events]
  if (study_stop_time < max(rec_time)) {
    warning("Target number of events reached before all subjects were enrolled.")
  }

  tt_pfs <- pmax(0, pmin(tt_pfs_cens1_calendar, study_stop_time) - rec_time)
  ev_pfs <- ifelse(tt_pfs_cens1_calendar <= study_stop_time, ev_pfs_cens1, 0)
  result <- data.frame(tt_pfs = tt_pfs, ev_pfs = ev_pfs)
  if (add_uncensored_pfs) {
    result$tt_pfs_uncens <- tt_pfs_uncens
  }
  result
}

#' Simulate Covariates and Progression Free Survival Data
#'
#' This combines the covariates simulation via [simul_covariates()] with 10
#' categorical covariates, and the PFS simulation via [simul_pfs()].
#'
#' @details
#' Regression coefficients are for an AFT with over-parametrized dummy
#' coding for arm-subgroup interactions.
#'
#' @param n (`count`)\cr number of patients.
#' @param coefs (`numeric`)\cr named vector of coefficients to set.
#' @param add_interaction (`flag`)\cr whether to add interaction terms between covariates
#'   1 and 2.
#' @param \dots additional parameters apart from the linear predictor values
#'   needed for [simul_pfs()].
#'
#' @return A combined `data.frame` with the `id` column, the design matrix and the
#'    PFS outcomes.
#' @export
#'
#' @examples
#' set.seed(321)
#' simul_data(
#'   n = 100,
#'   coefs = c(arm1 = 1),
#'   sigma_aft = 1,
#'   recr_duration = 0.2,
#'   rate_cens = 2,
#'   n_events = 20
#' )
simul_data <- function(n,
                       add_interaction = FALSE,
                       coefs,
                       ...) {
  assert_flag(add_interaction)
  assert_numeric(coefs, min.len = 1L, names = "unique")
  covariates <- simul_covariates(n = n, p_catvar = 10, add_contvars = FALSE, arm_factor = TRUE)
  subgroup_model <- ~ x_1 + x_2 + x_3 + x_4 + x_5 + x_6 + x_7 + x_8 + x_9 + x_10
  if (add_interaction) {
    covariates$x_1_2 <- factor(
      with(
        covariates,
        paste(as.character(x_1), as.character(x_2), sep = "")
      )
    )
    subgroup_model <- stats::update(subgroup_model, ~ . + x_1_2)
  }
  design_main <- stats::model.matrix(
    stats::update(subgroup_model, ~ arm + .),
    data = covariates
  )
  subgroup_vars <- all.vars(subgroup_model)
  design_ia <- NULL
  for (j in subgroup_vars) {
    ia_j <- stats::model.matrix(
      stats::as.formula(paste("~", j, "-1")),
      data = covariates
    ) * design_main[, "arm1"]
    design_ia <- cbind(design_ia, ia_j)
  }
  colnames(design_ia) <- paste(colnames(design_ia), "arm", sep = "_")
  colnames(design_ia) <- gsub(" ", "", colnames(design_ia))
  design_matrix <- cbind(design_main, design_ia)

  if (add_interaction) {
    # Remove created variable again such that final covariates dataset is identical for all scenarios.
    covariates$x_1_2 <- NULL
  }

  reg_coef <- rep(0, ncol(design_matrix))
  names(reg_coef) <- colnames(design_matrix)
  assert_subset(names(coefs), names(reg_coef))
  reg_coef[names(coefs)] <- coefs

  lp_aft <- design_matrix %*% reg_coef
  outcome <- simul_pfs(lp_aft = lp_aft, ...)
  d <- cbind(id = 1:n, covariates, outcome)
  d
}


# Simulate data sets from a given scenario, defined by named vector of coefficients in
# `coefs`.
simul_scenario <- function(scenario = c("1", "2", "3", "4", "5", "6"),
                           n_datasets = 1000,
                           n_patients = 1000 * inflation_factor,
                           n_events = 247 * inflation_factor,
                           sigma_aft = 0.85,
                           # Recruitment duration in years:
                           recr_duration = 3,
                           # Censoring rate at 1 year:
                           rate_cens = 0.02,
                           inflation_factor = 1,
                           add_uncensored_pfs = FALSE) {
  scenario <- match.arg(scenario)
  assert_count(n_datasets)

  # Set constant (across scenarios) intercept and prognostic factors.
  # For TTE AFT model: positive coef = longer survival
  constant_coefs <- c(
    "(Intercept)" = 2,
    "x_4c" = -log(0.7) * sigma_aft,   # Beneficial (HR=0.7)
    "x_6b" = -log(1.5) * sigma_aft    # Harmful (HR=1.5)
  )

  # Names of all the subgroup specific coefficients.
  group_coefs_names <- c(
    "x_1a_arm", "x_1b_arm", "x_2a_arm", "x_2b_arm",
    "x_3a_arm", "x_3b_arm", "x_4a_arm", "x_4b_arm", "x_4c_arm", "x_5a_arm",
    "x_5b_arm", "x_5c_arm", "x_5d_arm", "x_6a_arm", "x_6b_arm", "x_7a_arm",
    "x_7b_arm", "x_8a_arm", "x_8b_arm", "x_8c_arm", "x_9a_arm", "x_9b_arm",
    "x_10a_arm", "x_10b_arm", "x_10c_arm"
  )

  coefs <- switch(scenario,
                  # Positive trial, homogeneous treatment effect.
                  "1" = c(
                    constant_coefs,
                    "arm1" = -log(0.66) * sigma_aft
                  ),
                  # Overall HR~0.66, but no effect in x_4a.
                  "2" = c(
                    constant_coefs,
                    "arm1" = -log(0.66) * sigma_aft,
                    # No effect in x_4a:
                    "x_4a_arm" = log(0.66) * sigma_aft,
                    # Slightly enhanced effect in x_4b and x_4c to "compensate" no effect in x_4a:
                    "x_4b_arm" = -log(0.8) * sigma_aft,
                    "x_4c_arm" = -log(0.8) * sigma_aft
                  ),
                  # Overall HR~1, but HR~0.5 in x_4a.
                  "3" = c(
                    constant_coefs,
                    "arm1" = 0,
                    "x_4a_arm" = -log(0.5) * sigma_aft,
                    # Detrimental effect in x_4b and x_4c to "compensate" effect in x_4a:
                    "x_4b_arm" = -log(1.25) * sigma_aft,
                    "x_4c_arm" = -log(1.25) * sigma_aft
                  ),
                  # Mild heterogeneity.
                  "4" = {
                    set.seed(5)
                    c(
                      constant_coefs,
                      "arm1" = 0,
                      setNames(
                        -rnorm(25, sd = 0.15) * sigma_aft,
                        group_coefs_names
                      )
                    )
                  },
                  # Large heterogeneity.
                  "5" = {
                    set.seed(5)
                    c(
                      constant_coefs,
                      arm1 = 0,
                      setNames(
                        -rnorm(25, sd = 0.3) * sigma_aft,
                        group_coefs_names
                      )
                    )
                  },
                  # Model with interaction.
                  "6" = c(
                    constant_coefs,
                    "arm1" = -log(0.66) * sigma_aft,
                    "x_1_2aa_arm" = -log(1.5) * sigma_aft,
                    "x_1_2ba_arm" = -log(0.5) * sigma_aft,
                    "x_1_2ab_arm" = -log(0.92) * sigma_aft,
                    "x_1_2bb_arm" = -log(1.07) * sigma_aft
                  )
  )
  add_interaction <- scenario == "6"

  replicate(
    n_datasets,
    simul_data(
      n = n_patients,
      coefs = coefs,
      add_interaction = add_interaction,
      sigma_aft = sigma_aft,
      recr_duration = recr_duration,
      rate_cens = rate_cens,
      n_events = n_events,
      add_uncensored_pfs = add_uncensored_pfs
    ),
    simplify = FALSE
  )
}


# Initialize a `data.frame` with given row and column names.
init_data_frame <- function(row_names, col_names) {
  assert_character(row_names)
  assert_character(col_names)
  data.frame(matrix(
    nrow = length(row_names),
    ncol = length(col_names),
    dimnames = list(
      row_names,
      col_names
    )
  ))
}


# Sanitize subgroup string format.
sanitize_subgroups <- function(subgroups) {
  assert_character(subgroups)
  # Simply replace x_ with S_ and keep the rest (including dots)
  gsub(pattern = "^x_", replacement = "S_", x = subgroups)
}

#' Generate Stacked Data for Subgroup Analysis
#'
#' Creates a long-format dataset where each patient appears once for each
#' subgroup they belong to, with a 'subgroup' indicator.
#'
#' @param base_model Formula for the base model
#' @param subgr_model Formula with subgrouping variables
#' @param data Original data frame
#' @param resptype Type of response ("survival", "binary", "count", "continuous")
#'
#' @return A stacked data frame with subgroup indicators
#'
generate_stacked_data <- function(base_model, subgr_model, data, resptype) {

  # Extract variable names from formulas
  base_vars <- all.vars(base_model)
  subgr_vars <- all.vars(subgr_model)

  # Get the treatment variable name to preserve factor levels
  trt_var <- if (resptype == "survival") {
    base_vars[length(base_vars)]
  } else {
    base_vars[2]
  }

  # For survival, extract time and status variable names
  if (resptype == "survival") {
    time_var <- base_vars[1]
    status_var <- base_vars[2]
  }

  # Store the original arm factor levels
  arm_levels <- levels(data[[trt_var]])

  # Initialize list to hold stacked data for each subgroup
  stacked_list <- list()

  # For each subgrouping variable
  for (var in subgr_vars) {
    # Get the levels of this variable
    var_levels <- levels(data[[var]])

    # For each level, create a subset of data for that subgroup
    for (level in var_levels) {
      subgroup_data <- data[data[[var]] == level, ]

      # Skip subgroups with no observations
      if (nrow(subgroup_data) == 0) next

      # Skip subgroups that don't have both treatment arms
      # (needed to fit a model with treatment effect)
      trt_arms_present <- unique(subgroup_data[[trt_var]])
      if (length(trt_arms_present) < 2) {
        warning(sprintf("Subgroup %s.%s has only one treatment arm - skipping", var, level))
        next
      }

      # Add subgroup identifier
      subgroup_data$subgroup <- paste0(var, ".", level)

      # Rename variables for the naive function
      if (resptype == "survival") {
        # For survival: create standardized time, status, arm columns
        subgroup_data$time <- subgroup_data[[time_var]]
        subgroup_data$status <- subgroup_data[[status_var]]
        subgroup_data$arm <- factor(subgroup_data[[trt_var]], levels = arm_levels)
      } else {
        # For binary/count/continuous: standardize to y and arm
        resp_var <- base_vars[1]
        subgroup_data$y <- subgroup_data[[resp_var]]
        subgroup_data$arm <- factor(subgroup_data[[trt_var]], levels = arm_levels)
      }

      stacked_list[[paste0(var, ".", level)]] <- subgroup_data
    }
  }

  # Combine all subgroups into one stacked data frame
  stacked_data <- do.call(rbind, stacked_list)
  rownames(stacked_data) <- NULL

  return(stacked_data)
}

#' Naive Model Estimation (Extended for 4 Endpoints)
#'
#' Fits a separate model to the data for each subgroup.
#'
#' @param resp (`string`)\cr The response variable name.
#' @param trt (`string`)\cr The treatment variable name (factor, 2 levels).
#' @param subgr (`character`)\cr Vector of subgrouping variable names.
#' @param data (`data frame`)\cr The data frame.
#' @param resptype (`string`)\cr Type of data: "survival", "binary",
#'   "count", or "continuous".
#' @param status (`string`)\cr For "survival" only, the status variable name.
#'
#' @return List with `fit`, `estimates`, `model`, `resptype`, `data`.
#'
naive <- function(resp, trt, subgr, data,
                  resptype = c("survival", "binary", "count", "continuous"),
                  status = NULL) {

  # --- 1. Assertions ---
  checkmate::assert_string(resp)
  checkmate::assert_string(trt)
  checkmate::assert_character(subgr)
  checkmate::assert_data_frame(data)
  checkmate::assert_factor(data[[trt]])
  resptype <- match.arg(resptype)
  subgr_model <- stats::as.formula(paste("~", paste0(subgr, collapse = "+")))

  # --- 2. Fit Models ---
  if (resptype == "survival") {
    checkmate::assert_string(status)
    base_model <- stats::as.formula(paste("Surv(", resp, ",", status, ") ~ ", trt))

    stacked_data <- generate_stacked_data(base_model, subgr_model, data, resptype = resptype)
    list_subg <- split(stacked_data, ~subgroup)

    fit <- lapply(list_subg, FUN = function(d) {
      survival::coxph(survival::Surv(time, status) ~ arm, data = d)
    })
    # Keep dots in subgroup names for consistency with truth data
    # names(fit) <- gsub("\\.", "", names(fit))

    # Robust extraction for survival
    naive_estimates <- dplyr::bind_rows(lapply(fit, broom::tidy), .id = "subgroup") %>%
      dplyr::filter(term == "arm1") # Assuming arm is factor 0/1

  } else {
    # --- NON-SURVIVAL BLOCK ---
    base_model <- stats::as.formula(paste(resp, " ~ ", trt))

    # FIX 1: Pass the correct resptype, do not hardcode "binary"
    stacked_data <- generate_stacked_data(base_model, subgr_model, data,
                                          resptype = resptype)

    list_subg <- split(stacked_data, ~subgroup)

    model_fun <- switch(resptype,
                        "binary"     = function(d) stats::glm(y ~ arm, data = d, family = "binomial"),
                        "count"      = function(d) { MASS::glm.nb(y ~ arm, data = d) },
                        "continuous" = function(d) stats::lm(y ~ arm, data = d)
    )

    fit <- lapply(list_subg, FUN = model_fun)

    # FIX 2: Robust Extraction using dplyr
    # We stack all results first, then filter for the treatment term.
    # This handles cases where a subgroup drops the arm term (singular fit).
    naive_estimates <- dplyr::bind_rows(lapply(fit, broom::tidy), .id = "subgroup") %>%
      dplyr::filter(term == paste0("arm", levels(data[[trt]])[2]))
    # Logic: If trt levels are 0,1, term is usually "arm1".
    # If your data uses different factor levels, adjust "term" check accordingly.
    # Alternatively: filter(term != "(Intercept)") if no other covariates exist.
  }

  # --- 3. Return Results ---
  result <- list(
    fit = fit,
    estimates = naive_estimates,
    model = "naive",
    resptype = resptype,
    data = data
  )
  class(result) <- c("bonsaiforest", "naive")
  return(result)
}


subgroup_method_endpoint <- function(df, simul_no) {
  assert_data_frame(df)
  assert_count(simul_no)

  df$arm <- factor(df$arm) # Ensure arm is a factor

  # Call the naive function (loaded from functions.R)
  # Use all 10 variables to match data generation
  model <- naive(
    resp = endpoint_params$resp,
    trt = "arm",
    subgr = c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10"),
    data = df,
    resptype = endpoint_params$resptype,
    status = endpoint_params$status
  )

  # Format the results
  est_df <- model$estimates

  # Calculate Wald CIs (estimate +/- 1.96 * std.error)
  est_df$lower_ci <- est_df$estimate - 1.96 * est_df$std.error
  est_df$upper_ci <- est_df$estimate + 1.96 * est_df$std.error

  # Return a clean data.frame
  data.frame(
    simul_no = simul_no,
    estimator = "subgroup",
    subgroup = sanitize_subgroups(est_df$subgroup),
    estimate = est_df$estimate,
    std.error = est_df$std.error,
    lower_ci = est_df$lower_ci,
    upper_ci = est_df$upper_ci,
    p.value = est_df$p.value
  )
}

# Function constructor to return a function that analyses
# a single scenario with a `fun_method`.
fun_analysis <- function(fun_method) {
  function(scenario, scenario_no) {
    assert_list(scenario)
    assert_character(scenario_no)
    simul_no <- seq_along(scenario)

    # Map the method function over the list of datasets
    results <- Map(f = fun_method, scenario, simul_no)

    results <- do.call(rbind, results)
    results <- cbind(scenario_no = scenario_no, results)
    results
  }
}


# Compute results across all scenarios for a specific method,
# with optional caching.
compute_results <- function(scenarios,
                            analyze,
                            cache = NULL,
                            scenario_no = as.character(seq_along(scenarios))) {
  assert_list(scenarios)
  assert_function(analyze)
  assert_string(cache, null.ok = TRUE)

  # Use 2 cores by default, or more if specified
  cores <- 1

  if (!is.null(cache) && file.exists(cache)) {
    message(paste("Loading cached results from:", cache))
    return(readRDS(cache))
  }

  message(paste("Computing results using", cores, "cores (this may take a while)..."))

  res <- parallel::mcmapply(
    FUN = function(x, y) analyze(x, y),
    x = scenarios,
    y = scenario_no,
    SIMPLIFY = FALSE,
    mc.cores = cores
  )

  res <- do.call(rbind, res)

  if (!is.null(cache)) {
    # Ensure directory exists before saving
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    message(paste("Saving results to:", cache))
    saveRDS(res, file = cache)
  }
  res
}


#' Naive Population Model Estimation (Overall Model)
#'
#' Function to fit a single model to the entire dataset.
#'
#' @param resp (`string`)\cr The response variable name.
#' @param trt (`string`)\cr The treatment variable name (factor, 2 levels).
#' @param data (`data frame`)\cr The data frame.
#' @param resptype (`string`)\cr Type of data: "survival", "binary",
#'   "count", or "continuous".
#' @param status (`string`)\cr For "survival" only, the status variable name.
#'
#' @return List with `fit` (the model object) and `estimates` (tidy output
#'   data frame with a `subgroup` column).
#'
naivepop <- function(resp, trt, data,
                     resptype = c("survival", "binary", "count", "continuous"),
                     status = NULL) {

  # --- 1. Assertions ---
  assert_string(resp)
  assert_string(trt)
  assert_data_frame(data)
  assert_factor(data[[trt]])
  resptype <- match.arg(resptype)

  # --- 2. Build model formula ---
  if (resptype == "survival") {

    assert_string(status)
    model_formula <- stats::as.formula(paste("Surv(", resp, ",", status, ") ~ ", trt))

    # 3. Fit model
    fit <- survival::coxph(model_formula, data = data)

  } else {

    # This handles binary, count, continuous
    model_formula <- stats::as.formula(paste(resp, " ~ ", trt))

    # 3. Fit model
    fit <- switch(resptype,
                  "binary"     = stats::glm(model_formula, data = data, family = "binomial"),
                  "count"      = {
                    # Try glm.nb with better initialization
                    tryCatch({
                      suppressWarnings(
                        MASS::glm.nb(model_formula, data = data,
                                    init.theta = 1.0,  # Use reasonable starting value
                                    control = glm.control(maxit = 100))
                      )
                    }, error = function(e) {
                      # Fallback to Poisson if glm.nb fails (should be rare with full population)
                      stats::glm(model_formula, data = data, family = "poisson")
                    })
                  },
                  "continuous" = stats::lm(model_formula, data = data)
    )
  }

  # --- 4. Extract estimates ---
  # Use broom::tidy to get the coefficient for the treatment arm
  est <- broom::tidy(fit)

  if (resptype == "survival") {
    # coxph only has one term
    arm_est_row <- est[1, ]
  } else {
    # lm/glm have Intercept, so arm is 2nd term
    arm_est_row <- est[est$term == "arm1" | est$term == "arm", ][1, ]
  }

  # Add the 'subgroup' column to match the 'naive' function's output structure
  arm_est <- cbind(
    subgroup = "Overall",
    arm_est_row
  )


  # --- 5. Return a consistent list structure ---
  result <- list(
    fit = fit,         # The raw model object
    estimates = arm_est, # The tidy data frame (now with 'subgroup' column)
    model = "naivepop",
    resptype = resptype,
    data = data
  )
  class(result) <- c("bonsaiforest", "naivepop") # Mimic your class
  return(result)
}


#' Method for a single data set (Robust and Dynamic Version)
#'
#' This version is now fully dynamic and uses the `endpoint_params`
#' defined in the main script.
#'
population_method_endpoint <- function(df, simul_no) {
  assert_data_frame(df)
  assert_count(simul_no)
  df$arm <- factor(df$arm)

  # --- 1. Get subgroup list dynamically from data ---
  # This list is used to "recycle" the single population
  # estimate, matching the format of the subgroup output.
  # Use all 10 variables to match data generation and subgroup_method_endpoint
  subgr_vars <- c("x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10")

  # Generate all subgroup names by iterating through variables and their levels
  all_subgroups <- character(0)
  for (var in subgr_vars) {
    if (var %in% names(df) && is.factor(df[[var]])) {
      var_levels <- levels(df[[var]])
      # Create subgroup names like x_1.a, x_1.b, etc.
      all_subgroups <- c(all_subgroups, paste0(var, ".", var_levels))
    }
  }

  # Sanitize: replace x_ with S_
  all_subgroups <- gsub("^x_", "S_", all_subgroups)

  # --- 2. Run model with robustness check ---
  # This now uses the DYNAMIC parameters from the main script
  model <- try(
    naivepop(
      resp = endpoint_params$resp,     # <-- FIXED (dynamic)
      trt = "arm",
      data = df,
      resptype = endpoint_params$resptype, # <-- FIXED (dynamic)
      status = endpoint_params$status     # <-- FIXED (dynamic)
    ),
    silent = TRUE
  )

  # --- 3. Check for model failure ---
  if (inherits(model, "try-error") || is.null(model) || is.null(model$fit)) {
    message(paste("Model failed for simul_no:", simul_no, "- Returning NA rows."))

    # Return NA rows in the *standardized format*
    return(data.frame(
      simul_no = simul_no,
      estimator = "population",
      subgroup = all_subgroups, # N rows
      estimate = NA_real_,
      std.error = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      p.value = NA_real_
    ))
  }

  # --- 4. Process successful model ---

  # 'naivepop' returns a tidy data frame with the 'Overall' subgroup
  estimates_df <- model$estimates

  # Get CI on the log scale (which is the default 'estimate' scale)
  # Use confint.default() to avoid slow profile likelihood computation for count data
  ci_log <- confint.default(model$fit)

  # Get the correct CI row (1 for coxph, 2 for glm/lm)
  conf_int_log_vector <- if (endpoint_params$resptype == "survival") ci_log[1, ] else ci_log[2, ]

  # --- 5. Create final data frame (using recycling) ---
  # This now returns the standardized log-scale estimates
  res <- data.frame(
    simul_no = simul_no,
    estimator = "population",
    subgroup = all_subgroups, # Replicate this row for all subgroups
    estimate = as.numeric(estimates_df$estimate),
    std.error = as.numeric(estimates_df$std.error),
    lower_ci = as.numeric(conf_int_log_vector[1]),
    upper_ci = as.numeric(conf_int_log_vector[2]),
    p.value = as.numeric(estimates_df$p.value),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(res)
}

#'

# Helper function for global model
run_global_model_task <- function(sim_id,
                                  model_type,
                                  prior_name,
                                  all_data,
                                  subgr_vars,
                                  log_path,
                                  prior_spec,
                                  endpoint_params) {

  # --- 1. Thread Safety ---
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  data.table::setDTthreads(threads = 1)

  # --- 2. Logging Setup ---
  task_description <- sprintf("SimID: %s | Model: %s | Prior: %s | Endpoint: %s",
                              sim_id, model_type, prior_name, endpoint_params$resptype)
  start_time <- Sys.time()

  cat(sprintf("[%s] STARTING: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S"), task_description),
      file = log_path, append = TRUE)

  # --- 3. Execute the Core Task ---
  results <- tryCatch({
    df <- all_data[[sim_id]]

    # --- Model Fitting Logic ---
    prognostic_str <- paste(subgr_vars, collapse = " + ")
    predictive_str <- paste(paste0("arm:", subgr_vars), collapse = " + ")

    # Set unique cmdstanr output directory per task to avoid conflicts when
    # running multiple jobs in parallel on ACE platform. Create job-specific
    # directory using environment variable or fallback to tempdir.
    job_temp_base <- Sys.getenv("JOB_TEMP_DIR", tempdir())
    safe_sim_id <- gsub("[^[:alnum:]_]", "_", sim_id)
    safe_prior_name <- gsub("[^[:alnum:]_]", "_", prior_name)

    # Add random suffix to cmdstan dir to prevent compilation collisions
    # when multiple jobs compile similar models simultaneously
    random_suffix <- paste0(sample(c(0:9, letters), 6, replace = TRUE), collapse = "")
    cmdstan_output_dir <- file.path(
      job_temp_base,
      sprintf("cmdstan_%s_%s_%s", safe_prior_name, safe_sim_id, random_suffix)
    )
    dir.create(cmdstan_output_dir, recursive = TRUE, showWarnings = FALSE)

    fit <- run_brms_analysis(
      data = df,
      response_formula = endpoint_params$resp_formula,  # Dynamic
      response_type = endpoint_params$resptype,             # Dynamic
      unshrunk_terms_formula  = paste("~", prognostic_str),
      shrunk_predictive_formula  = paste("~ 0 +", predictive_str),
      shrunk_predictive_prior =  prior_spec$prior,
      stanvars = prior_spec$stanvars,
      chains = 4, iter = 2000, warmup = 1000, cores = 1,
      backend = "cmdstanr",
      output_dir = cmdstan_output_dir
    )

    # Summarize results
    output <- summary_subgroup_effects(
      fit
    )

    # Clean up cmdstan output directory to save disk space
    unlink(cmdstan_output_dir, recursive = TRUE)

    output$estimates

  }, error = function(e) {
    clean_message <- gsub("\\n", " ", e$message)
    cat(sprintf("[%s] ERROR in %s: %s\n", Sys.time(), task_description, clean_message),
        file = log_path, append = TRUE)
    # Clean up cmdstan directory even on error
    if (exists("cmdstan_output_dir") && dir.exists(cmdstan_output_dir)) {
      unlink(cmdstan_output_dir, recursive = TRUE)
    }
    return(tibble(error = clean_message))
  })

  # --- 4. Final Logging & Formatting ---
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("[%s] FINISHED: %s | Duration: %.2f mins\n", format(end_time, "%Y-%m-%d %H:%M:%S"), task_description, duration),
      file = log_path, append = TRUE)

  if (is.data.frame(results) && nrow(results) > 0) {
    results %>%
      mutate(
        model_type = model_type,
        prior_name = prior_name,
        .before = 1
      )
  } else {
    tibble(
      model_type = model_type,
      prior_name = prior_name,
      error = ifelse(is.data.frame(results), "No results returned", results$error[1])
    )
  }
}

