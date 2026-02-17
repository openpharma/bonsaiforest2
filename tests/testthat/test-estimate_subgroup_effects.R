# tests/testthat/test-estimate_subgroup_effects.R

# --- Test Setup ---

# 1. Create Original Sample Data
original_test_data <- data.frame(
  outcome = rnorm(30),
  time = round(runif(30, 1, 50)),
  status = sample(0:1, 30, replace = TRUE),
  trt = sample(0:1, 30, replace = TRUE),
  age = rnorm(30, 50, 10),
  region = factor(sample(c("A", "B", "C"), 30, replace = TRUE)),
  sex = factor(sample(c("M", "F"), 30, replace = TRUE))
)
original_test_data$trt <- factor(original_test_data$trt, levels = c(0, 1))

# 2. Fit a *minimal* brms model ONCE for use in multiple tests
# This is slow, but simpler than full mocking for now.
# We include interactions so 'auto' detection can be tested.
minimal_brms_fit <- suppressMessages(
  run_brms_analysis(
    data = original_test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",  # Use ~ 0 + to avoid intercept warning
    shrunk_predictive_formula = "~ 0 + trt:region + trt:sex", # Interactions needed for 'auto', use ~ 0 + to avoid intercept warning
    chains = 1, iter = 10, warmup = 5, refresh = 0,
    backend = "cmdstanr", # Use faster backend if available
    cores = 1
  )
)

# --- Tests Start Here ---

test_that("Basic functionality works (overall and specific subgroups)", {
  # Test Overall effect - with automatic extraction from model
  res_overall <- estimate_subgroup_effects(
    brms_fit = minimal_brms_fit,
    subgroup_vars = NULL # Explicitly request overall
  )

  expect_s3_class(res_overall$estimates, "tbl_df")
  expect_s3_class(res_overall$draws, "tbl_df")
  expect_named(res_overall, c("estimates", "draws"))
  expect_equal(nrow(res_overall$estimates), 1)
  expect_equal(res_overall$estimates$Subgroup, "Overall")
  expect_equal(ncol(res_overall$draws), 1)
  expect_true("Overall" %in% names(res_overall$draws))
  expect_true(all(c("Subgroup", "Median", "CI_Lower", "CI_Upper") %in% names(res_overall$estimates)))

  # Test Specific subgroups - with automatic extraction
  res_specific <- estimate_subgroup_effects(
    brms_fit = minimal_brms_fit,
    subgroup_vars = c("region", "sex") # Request specific subgroups
  )

  expect_equal(nrow(res_specific$estimates), 5) # 3 region levels + 2 sex levels
  expect_equal(ncol(res_specific$draws), 5)
  expected_subgroups <- c("region: A", "region: B", "region: C", "sex: F", "sex: M")
  expect_equal(sort(res_specific$estimates$Subgroup), sort(expected_subgroups))
  expect_equal(sort(names(res_specific$draws)), sort(expected_subgroups))

  # Test explicit parameter override still works
  res_explicit <- estimate_subgroup_effects(
    brms_fit = minimal_brms_fit,
    trt_var = "trt",  # Explicitly provide
    subgroup_vars = NULL,
    response_type = "continuous"  # Explicitly provide
  )
  expect_equal(nrow(res_explicit$estimates), 1)
  expect_equal(res_explicit$estimates$Subgroup, "Overall")

})

test_that("'auto' subgroup detection works", {
  # minimal_brms_fit was created with trt:region and trt:sex interactions
  res_auto <- estimate_subgroup_effects(
    brms_fit = minimal_brms_fit,
    subgroup_vars = "auto" # Use auto-detection, other params auto-extracted
  )

  # Should detect 'region' and 'sex' from the interaction terms in minimal_brms_fit$data
  expect_equal(nrow(res_auto$estimates), 5) # 3 region levels + 2 sex levels
  expect_equal(ncol(res_auto$draws), 5)
  expected_subgroups <- c("region: A", "region: B", "region: C", "sex: F", "sex: M")
  expect_equal(sort(res_auto$estimates$Subgroup), sort(expected_subgroups))
  expect_equal(sort(names(res_auto$draws)), sort(expected_subgroups))
})


test_that("ndraws argument works", {
  # Get full draws first
  res_full <- estimate_subgroup_effects(
    brms_fit = minimal_brms_fit,
    subgroup_vars = "sex",
    ndraws = NULL
  )
  expect_equal(nrow(res_full$draws), 5) # Check based on iter/warmup

  # Get subset of draws
  res_subset <- estimate_subgroup_effects(
    brms_fit = minimal_brms_fit,
    subgroup_vars = "sex",
    ndraws = 3
  )
  expect_equal(nrow(res_subset$draws), 3) # Check if subsetting worked
})

# Note: Testing survival requires a survival model fit, which adds complexity/time.
# Add similar tests for 'binary' and 'count' if needed, fitting minimal models first.

test_that("Assertions catch invalid inputs", {
  # Invalid brms_fit
  expect_error(
    estimate_subgroup_effects(brms_fit = list()),
    regexp = "Must inherit from class 'brmsfit'"
  )

  # Invalid trt_var (not a string) when explicitly provided
  expect_error(
    estimate_subgroup_effects(brms_fit = minimal_brms_fit, trt_var = 123),
    regexp = "Must be of type 'string'"
  )

  # trt_var not in data when explicitly provided
  expect_error(
    estimate_subgroup_effects(brms_fit = minimal_brms_fit, trt_var = "treatment"),
    regexp = "Must be a subset of"
  )

  # Invalid subgroup_vars type
  expect_error(
    estimate_subgroup_effects(brms_fit = minimal_brms_fit, subgroup_vars = 123),
    regexp = "Assertion failed"
  )

  # subgroup_vars not in data
  expect_error(
    estimate_subgroup_effects(brms_fit = minimal_brms_fit, subgroup_vars = c("region", "missing_var")),
    regexp = "Must be a subset of"
  )

  # Invalid response_type when explicitly provided
  expect_error(
    estimate_subgroup_effects(brms_fit = minimal_brms_fit, response_type = "gaussian"),
    regexp = "Must be element of set"
  )

  # Invalid ndraws
  expect_error(
    estimate_subgroup_effects(brms_fit = minimal_brms_fit, ndraws = 0),
    regexp = "Must be >= 1"
  )

})
