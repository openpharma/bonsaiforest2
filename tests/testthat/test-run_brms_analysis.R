# tests/testthat/test-run_brms_analysis.R

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Helper to run minimal brms models for testing
run_quick_full_analysis <- function(data, ...) {
  suppressMessages(
    run_brms_analysis(
      data = data,
      ...,
      chains = 1, iter = 1000, warmup = 50, refresh = 0,
      backend = "cmdstanr",
      cores = 1
    )
  )
}

# ============================================================================
# TEST DATA SETUP
# ============================================================================

set.seed(42)
test_data_run <- data.frame(
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
# BASIC INTEGRATION TESTS (GLOBAL SYNTAX - COLON)
# ============================================================================

test_that("run_brms_analysis works with GLOBAL fixed effects (continuous outcome)", {
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",
    shrunk_predictive_formula = "~ 0 + trt:region",
    sigma_ref = sd(test_data_run$outcome)
  )

  expect_s3_class(fit, "brmsfit")

  # Verify priors were set for shrunk components
  priors_df <- as.data.frame(fit$prior)
  expect_true(any(priors_df$nlpar == "shprogeffect"))
  expect_true(any(priors_df$nlpar == "shpredeffect"))
})

test_that("run_brms_analysis works with GLOBAL fixed effects (binary outcome)", {
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "binary_outcome ~ trt",
    response_type = "binary",
    shrunk_predictive_formula = "~ 0 + trt:region",
    sigma_ref = 1  # Default for binary
  )

  expect_s3_class(fit, "brmsfit")
})

test_that("run_brms_analysis works with GLOBAL fixed effects (survival outcome)", {
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    unshrunk_terms_formula = "~ age",
    shrunk_predictive_formula = "~ 0 + trt:subgroup",
    sigma_ref = 1  # Default for survival
  )

  expect_s3_class(fit, "brmsfit")
})

# ============================================================================
# OVAT SYNTAX TESTS (RANDOM EFFECTS - PIPE-PIPE)
# ============================================================================

test_that("run_brms_analysis works with OVAT random effects (continuous outcome)", {
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age",  # Put age in unshrunk instead
    shrunk_predictive_formula = "~ (0 + trt || subgroup)",
    sigma_ref = sd(test_data_run$outcome)
  )

  expect_s3_class(fit, "brmsfit")

  # Verify random effects priors were set correctly
  priors_df <- as.data.frame(fit$prior)
  sd_priors <- priors_df[priors_df$class == "sd", ]
  expect_true(nrow(sd_priors) > 0, info = "Expected SD priors for random effects")
})

test_that("run_brms_analysis works with OVAT random effects (survival outcome)", {
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    unshrunk_terms_formula = "~ age",
    shrunk_predictive_formula = "~ (0 + trt || subgroup)",
    sigma_ref = 1,
    shrunk_predictive_prior = NULL  # No fixed effects in random effects formula
  )

  expect_s3_class(fit, "brmsfit")
})

# ============================================================================
# DEFAULT BEHAVIOR TESTS
# ============================================================================

test_that("run_brms_analysis uses default sigma_ref = 1", {
  # Should work without specifying sigma_ref (defaults to 1)
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "binary_outcome ~ trt",
    response_type = "binary",
    shrunk_predictive_formula = "~ 0 + trt:region"
    # Note: sigma_ref not specified, should default to 1
  )

  expect_s3_class(fit, "brmsfit")
})

test_that("run_brms_analysis warns about default sigma_ref for continuous outcomes", {
  # Don't use helper function that suppresses messages
  expect_message(
    run_brms_analysis(
      data = test_data_run,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      shrunk_predictive_formula = "~ 0 + trt:region",
      # sigma_ref not specified - should warn
      chains = 1, iter = 1000, warmup = 50, refresh = 0,
      backend = "cmdstanr",
      cores = 1
    ),
    "Using default sigma_ref = 1 for continuous outcome"
  )
})

# ============================================================================
# CUSTOM PRIOR TESTS
# ============================================================================

test_that("run_brms_analysis accepts custom priors with sigma_ref placeholder", {
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",
    shrunk_predictive_formula = "~ 0 + trt:region",
    sigma_ref = 5,
    unshrunk_prior = "normal(0, 2.5 * sigma_ref)",
    shrunk_prognostic_prior = "normal(0, 1)",
    shrunk_predictive_prior = "horseshoe(scale_global = sigma_ref)"
  )

  expect_s3_class(fit, "brmsfit")

  # Verify sigma_ref was substituted in priors
  priors_df <- as.data.frame(fit$prior)
  # Check for horseshoe with scale_global = 5 in shrunk predictive
  shpred_priors <- priors_df[priors_df$nlpar == "shpredeffect", "prior"]
  expect_true(any(grepl("horseshoe.*5", shpred_priors, ignore.case = TRUE)))
})

# ============================================================================
# VALIDATION TESTS
# ============================================================================

test_that("run_brms_analysis validates inputs correctly", {
  # Invalid data (not a data frame)
  expect_error(
    run_brms_analysis(
      data = as.matrix(test_data_run),
      response_formula = "outcome ~ trt",
      response_type = "continuous"
    ),
    regexp = "Must be of type 'data.frame'"
  )

  # Invalid response_formula (string without ~)
  expect_error(
    run_brms_analysis(
      data = test_data_run,
      response_formula = "outcome",
      response_type = "continuous"
    ),
    regexp = "Must comply to pattern '~'"
  )

  # Invalid response_type
  expect_error(
    run_brms_analysis(
      data = test_data_run,
      response_formula = "outcome ~ trt",
      response_type = "gaussian"
    ),
    regexp = "Must be element of set"
  )

  # Invalid sigma_ref (negative)
  expect_error(
    run_brms_analysis(
      data = test_data_run,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      sigma_ref = -1
    ),
    regexp = "Element 1 is not >= 0"
  )

  # Invalid stanvars class
  expect_error(
    run_brms_analysis(
      data = test_data_run,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      stanvars = list("dummy")
    ),
    regexp = "Must inherit from class 'stanvars'"
  )
})

# ============================================================================
# MIXED SYNTAX TESTS
# ============================================================================

test_that("run_brms_analysis works with mixed GLOBAL and OVAT syntax", {
  # Use GLOBAL for prognostic, OVAT for predictive
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",  # GLOBAL (fixed effect)
    shrunk_predictive_formula = "~ (0 + trt || subgroup)",  # OVAT (random effect)
    sigma_ref = sd(test_data_run$outcome),
    shrunk_predictive_prior = NULL  # No fixed effects in random effects formula
  )

  expect_s3_class(fit, "brmsfit")
})

# ============================================================================
# ATTRIBUTE PRESERVATION TESTS
# ============================================================================

test_that("run_brms_analysis preserves attributes in fitted model", {
  fit <- run_quick_full_analysis(
    data = test_data_run,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region",
    sigma_ref = sd(test_data_run$outcome)
  )

  # Check that attributes are stored for downstream functions
  expect_equal(attr(fit, "response_type"), "continuous")
  expect_s3_class(attr(fit, "model_data"), "data.frame")
  expect_type(attr(fit, "trt_var"), "character")
  expect_equal(attr(fit, "trt_var"), "trt")
})
