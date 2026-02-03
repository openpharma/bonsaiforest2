# tests/testthat/test-fit_brms_model.R

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Helper to run minimal brms models for testing
run_quick_brm <- function(prepared_model, sigma_ref = 1, ...) {
  suppressMessages(
    fit_brms_model(
      prepared_model = prepared_model,
      sigma_ref = sigma_ref,
      ...,
      chains = 1, iter = 1000, warmup = 50, refresh = 0,
      backend = "cmdstanr",
      cores = 1
    )
  )
}

# Helper to extract prior for specific nlpar/class/coef
get_prior_string <- function(fit, nlpar_name, class_name = "b", coef_name = "") {
  df <- as.data.frame(fit$prior)
  rows <- which(df$nlpar == nlpar_name & df$class == class_name & df$coef == coef_name)
  if (length(rows) == 0) return(NA_character_)
  prior_val <- df[rows[1], "prior"]
  if (is.na(prior_val) || prior_val == "") return(NA_character_)
  return(prior_val)
}

# ============================================================================
# TEST DATA SETUP
# ============================================================================

set.seed(42)
test_data_fit <- data.frame(
  outcome = rnorm(20, mean = 10, sd = 5),
  binary_outcome = sample(0:1, 20, replace = TRUE),
  count_outcome = rpois(20, lambda = 5),
  time = round(runif(20, 1, 50)),
  status = sample(0:1, 20, replace = TRUE),
  trt = factor(sample(0:1, 20, replace = TRUE), levels = c(0, 1)),
  age = rnorm(20, 50, 10),
  region = factor(sample(c("A", "B"), 20, replace = TRUE)),
  subgroup = factor(sample(c("S1", "S2", "S3"), 20, replace = TRUE))
)

# ============================================================================
# BASIC MODEL FITTING TESTS
# ============================================================================

test_that("Continuous response type works", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = sd(test_data_fit$outcome))
  expect_s3_class(fit, "brmsfit")
})

test_that("Binary response type works", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "binary_outcome ~ trt",
    response_type = "binary",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)
  expect_s3_class(fit, "brmsfit")
})

test_that("Count response type works", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "count_outcome ~ trt",
    response_type = "count",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)
  expect_s3_class(fit, "brmsfit")
})

test_that("Survival response type works", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    unshrunk_terms_formula = "~ age"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)
  expect_s3_class(fit, "brmsfit")
})

# ============================================================================
# PRIOR SPECIFICATION TESTS
# ============================================================================

test_that("Default priors are correct for continuous models", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  sigma_ref_val <- 2.5
  fit <- run_quick_brm(prep, sigma_ref = sigma_ref_val)

  # Check defaults: horseshoe(1) for shrunk terms
  expect_match(get_prior_string(fit, "shprogeffect"), "horseshoe\\(1\\)")
  expect_match(get_prior_string(fit, "shpredeffect"), "horseshoe\\(1\\)")

  # Check unshrunk: normal(0, 2.5 * sigma_ref)
  expected_sd <- 2.5 * sigma_ref_val
  expect_match(get_prior_string(fit, "unshrunktermeffect"),
               paste0("normal\\(0,\\s*", expected_sd, "\\)"))
})

test_that("Default priors are correct for binary models", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "binary_outcome ~ trt",
    response_type = "binary",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)

  # Binary: normal(0, 1.5) for unshrunk
  expect_match(get_prior_string(fit, "unshrunktermeffect"), "normal\\(0,\\s*1\\.5\\)")
  expect_match(get_prior_string(fit, "shpredeffect"), "horseshoe\\(1\\)")
})

test_that("Default priors are correct for count models", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "count_outcome ~ trt",
    response_type = "count",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)

  # Count: normal(0, 2) for unshrunk
  expect_match(get_prior_string(fit, "unshrunktermeffect"), "normal\\(0,\\s*2\\)")
  expect_match(get_prior_string(fit, "shpredeffect"), "horseshoe\\(1\\)")
})

test_that("Default priors are correct for survival models", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)

  # Survival: normal(0, 1.5) for unshrunk
  expect_match(get_prior_string(fit, "unshrunktermeffect"), "normal\\(0,\\s*1\\.5\\)")
  expect_match(get_prior_string(fit, "shpredeffect"), "horseshoe\\(1\\)")
})

test_that("User-specified priors work (string format)", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(
    prep,
    sigma_ref = 2,
    unshrunk_prior = "normal(0, 5)",
    shrunk_prognostic_prior = "normal(0, 1)",
    shrunk_predictive_prior = "normal(0, 0.5)"
  )

  expect_match(get_prior_string(fit, "unshrunktermeffect"), "normal\\(0,\\s*5\\)")
  expect_match(get_prior_string(fit, "shprogeffect"), "normal\\(0,\\s*1\\)")
  expect_match(get_prior_string(fit, "shpredeffect"), "normal\\(0,\\s*0\\.5\\)")
})

test_that("sigma_ref substitution works in prior strings", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  sigma_ref_val <- 3
  fit <- run_quick_brm(
    prep,
    sigma_ref = sigma_ref_val,
    unshrunk_prior = "normal(0, 2 * sigma_ref)",
    shrunk_predictive_prior = "normal(0, sigma_ref)"
  )

  # Check that "sigma_ref" was replaced with actual value (string substitution, not evaluated)
  expect_match(get_prior_string(fit, "shpredeffect"),
               "normal\\(0,\\s*3\\)")

})

test_that("Intercept centering at outcome mean works", {
  # Suppress warning about intercept in shrunk formula - this is intentional for this test
  suppressWarnings({
    prep <- prepare_formula_model(
      data = test_data_fit,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      shrunk_predictive_formula = ~ 0 + trt:region,
      shrunk_prognostic_formula = ~ region  # Intentionally includes intercept to test shrinkage
    )
  })

  outcome_mean <- mean(test_data_fit$outcome)
  sigma_ref_val <- 2

  fit <- run_quick_brm(prep, sigma_ref = sigma_ref_val)

  # Default intercept prior should be centered at outcome mean
  intercept_prior <- get_prior_string(fit, "unshrunktermeffect", "b", "Intercept")
  # Use flexible regex to match the outcome mean (allows for floating point precision)
  expect_match(intercept_prior, "normal\\([0-9]+\\.[0-9]+,\\s*[0-9]+")
  # Also verify it's close to the expected value
  expect_true(grepl(as.character(round(outcome_mean, 1)), intercept_prior))
  custom_prior <- brms::set_prior("student_t(3, 0, 2.5)", class = "b")

  fit <- run_quick_brm(
    prep,
    sigma_ref = 2,
    shrunk_prognostic_prior = custom_prior
  )

  expect_match(get_prior_string(fit, "shprogeffect"), "student_t\\(3,\\s*0,\\s*2\\.5\\)")
})

test_that("Complex priors work (R2D2)", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age"
  )

  fit <- run_quick_brm(
    prep,
    sigma_ref = 2,
    shrunk_prognostic_prior = "R2D2(mean_R2 = 0.5, prec_R2 = 1)"
  )

  expect_match(get_prior_string(fit, "shprogeffect"), "R2D2\\(")
})

# ============================================================================
# INTERCEPT HANDLING TESTS
# ============================================================================

test_that("Models with no intercept work (~ 0 + ...)", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",  # trt in response formula
    response_type = "continuous",
    unshrunk_terms_formula = "~ 0 + age",  # Use different variable, explicitly no intercept
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 2)

  # Should not have intercept prior when unshrunk_terms has ~ 0 +
  intercept_prior <- get_prior_string(fit, "unshrunktermeffect", "b", "Intercept")
  expect_true(is.na(intercept_prior))
})

test_that("Survival models have no intercept", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)

  # Cox models should not have intercept
  intercept_prior <- get_prior_string(fit, "unshrunktermeffect", "b", "Intercept")
  expect_true(is.na(intercept_prior))
})

# ============================================================================
# RANDOM EFFECTS TESTS
# ============================================================================

test_that("Random effects (|| syntax) work", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + (1 + trt || subgroup)"  # Use ~ 0 + to remove intercept
  )

  fit <- run_quick_brm(prep, sigma_ref = 2, shrunk_predictive_prior = NULL)
  expect_s3_class(fit, "brmsfit")

  # Check that random effects get priors on sd class
  df <- as.data.frame(fit$prior)
  sd_priors <- df[df$class == "sd" & df$nlpar == "shpredeffect", ]

  # Should have priors for Intercept and trt
  expect_true(nrow(sd_priors) >= 2)
  # Check that priors were assigned (may be empty string for brms defaults)
  expect_true(nrow(sd_priors) > 0)
})

test_that("Mixed notation (random + fixed effects) works", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + (1 + trt || subgroup) + trt:region"  # Add ~ 0 +
  )

  fit <- run_quick_brm(prep, sigma_ref = 2)
  expect_s3_class(fit, "brmsfit")

  df <- as.data.frame(fit$prior)

  # Should have both sd (random) and b (fixed) priors for shpredeffect
  sd_priors <- df[df$class == "sd" & df$nlpar == "shpredeffect", ]
  b_priors <- df[df$class == "b" & df$nlpar == "shpredeffect", ]

  expect_true(nrow(sd_priors) > 0)
  expect_true(nrow(b_priors) > 0)
})

# ============================================================================
# STRATIFICATION TESTS
# ============================================================================

test_that("Stratified models work", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region",
    stratification_formula = "~ region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 2)
  expect_s3_class(fit, "brmsfit")
})

# ============================================================================
# ATTRIBUTE ATTACHMENT TESTS
# ============================================================================

test_that("Attributes are attached correctly", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 2)

  # Check response_type attribute
  expect_equal(attr(fit, "response_type"), "continuous")

  # Check model_data attribute
  expect_s3_class(attr(fit, "model_data"), "data.frame")

  # Check trt_var attribute
  expect_equal(attr(fit, "trt_var"), "trt")
})

# ============================================================================
# MESSAGE/WARNING TESTS
# ============================================================================

test_that("Warning for default sigma_ref = 1 in continuous models", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  # Capture the message about default sigma_ref
  msgs <- capture.output(
    fit <- fit_brms_model(
      prep,
      sigma_ref = 1,  # Using default
      chains = 1, iter = 1000, warmup = 50, refresh = 0,
      backend = "cmdstanr", cores = 1
    ),
    type = "message"
  )

  # Check that the warning was issued
  expect_true(any(grepl("Using default sigma_ref = 1 for continuous outcome", msgs)))
})

test_that("No warning for non-continuous models with sigma_ref = 1", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "binary_outcome ~ trt",
    response_type = "binary",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  # Should NOT produce warning for binary outcome
  msgs <- capture.output(
    fit <- fit_brms_model(
      prep,
      sigma_ref = 1,
      chains = 1, iter = 1000, warmup = 50, refresh = 0,
      backend = "cmdstanr", cores = 1
    ),
    type = "message"
  )

  expect_false(any(grepl("Using default sigma_ref = 1 for continuous outcome", msgs)))
})

# ============================================================================
# STANVARS TESTS
# ============================================================================

test_that("Stanvars argument is accepted", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  sv <- brms::stanvar(scode = "real dummy_var;", block = "parameters")

  fit <- run_quick_brm(prep, sigma_ref = 2, stanvars = sv)
  expect_s3_class(fit, "brmsfit")
  expect_true(grepl("dummy_var", fit$model))
})

# ============================================================================
# INPUT VALIDATION TESTS
# ============================================================================

test_that("Invalid prepared_model type is rejected", {
  expect_error(
    fit_brms_model(prepared_model = "not a list", sigma_ref = 1),
    regexp = "Must be of type 'list'"
  )
})

test_that("Missing required elements in prepared_model are caught", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  bad_list <- list(formula = prep$formula, data = prep$data)
  expect_error(
    fit_brms_model(prepared_model = bad_list, sigma_ref = 1),
    regexp = "missing elements.*'response_type'"
  )
})

test_that("Invalid formula class is rejected", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  prep$formula <- "outcome ~ trt"
  expect_error(
    fit_brms_model(prepared_model = prep, sigma_ref = 1),
    regexp = "Must inherit from class 'brmsformula'"
  )
})

test_that("Invalid data class is rejected", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  prep$data <- as.matrix(prep$data)
  expect_error(
    fit_brms_model(prepared_model = prep, sigma_ref = 1),
    regexp = "Must be of type 'data.frame'"
  )
})

test_that("Invalid response_type is rejected", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  prep$response_type <- "gaussian"
  expect_error(
    fit_brms_model(prepared_model = prep, sigma_ref = 1),
    regexp = "Must be element of set"
  )
})

test_that("Negative sigma_ref is rejected", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  expect_error(
    fit_brms_model(prepared_model = prep, sigma_ref = -1),
    regexp = "Element 1 is not >= 0"
  )
})

test_that("Non-numeric sigma_ref is rejected", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  expect_error(
    fit_brms_model(prepared_model = prep, sigma_ref = "text"),
    regexp = "Must be of type 'number'"
  )
})

test_that("Invalid stanvars class is rejected", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  expect_error(
    fit_brms_model(prepared_model = prep, sigma_ref = 1, stanvars = list("dummy")),
    regexp = "Must inherit from class 'stanvars'"
  )
})

# ============================================================================
# SIGMA PRIOR TESTS
# ============================================================================

test_that("Sigma prior is added for non-stratified continuous models", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  sigma_ref_val <- 3
  fit <- run_quick_brm(prep, sigma_ref = sigma_ref_val)

  df <- as.data.frame(fit$prior)
  sigma_prior <- df[df$class == "sigma", ]

  expect_true(nrow(sigma_prior) > 0)
  expect_match(sigma_prior$prior[1], paste0("student_t\\(3,\\s*0,\\s*", sigma_ref_val, "\\)"))
})

test_that("No sigma prior for binary models", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "binary_outcome ~ trt",
    response_type = "binary",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)

  df <- as.data.frame(fit$prior)
  sigma_prior <- df[df$class == "sigma", ]

  expect_equal(nrow(sigma_prior), 0)
})

test_that("No sigma prior for count models (uses shape instead)", {
  prep <- prepare_formula_model(
    data = test_data_fit,
    response_formula = "count_outcome ~ trt",
    response_type = "count",
    shrunk_predictive_formula = "~ 0 + trt:region"
  )

  fit <- run_quick_brm(prep, sigma_ref = 1)

  df <- as.data.frame(fit$prior)
  sigma_prior <- df[df$class == "sigma", ]

  expect_equal(nrow(sigma_prior), 0)
})
