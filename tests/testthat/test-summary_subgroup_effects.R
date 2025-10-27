# tests/testthat/test-summary_subgroup_effects.R

# --- Test Setup ---
# Reusing setup from estimate_subgroup_effects tests for consistency

# 1. Create Original Sample Data
original_test_data_summary <- data.frame(
  outcome = rnorm(30),
  trt = sample(0:1, 30, replace = TRUE),
  age = rnorm(30, 50, 10),
  region = factor(sample(c("A", "B", "C"), 30, replace = TRUE)),
  sex = factor(sample(c("M", "F"), 30, replace = TRUE))
)
original_test_data_summary$trt <- factor(original_test_data_summary$trt, levels = c(0, 1))

# 2. Fit a *minimal* brms model ONCE
# Includes interactions for 'auto' testing
minimal_brms_fit_summary <- suppressMessages(
  run_brms_analysis(
    data = original_test_data_summary,
    response_formula_str = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula_str = "~ age",
    shrunk_predictive_formula_str = "~ trt:region + trt:sex",
    chains = 1, iter = 10, warmup = 5, refresh = 0,
    backend = "cmdstanr", cores = 1
  )
)

# 3. Create a sample summary object ONCE for plot tests
sample_summary_obj <- suppressMessages(
  summary_subgroup_effects(
    brms_fit = minimal_brms_fit_summary,
    original_data = original_test_data_summary,
    trt_var = "trt",
    response_type = "continuous",
    subgroup_vars = c("region", "sex") # Use specific vars for predictable output
  )
)

# --- Tests for summary_subgroup_effects ---

test_that("summary_subgroup_effects works with subgroup_vars = NULL (Overall only)", {
  res_overall <- suppressMessages(
    summary_subgroup_effects(
      brms_fit = minimal_brms_fit_summary,
      original_data = original_test_data_summary,
      trt_var = "trt",
      response_type = "continuous",
      subgroup_vars = NULL
    )
  )

  expect_s3_class(res_overall, "subgroup_summary")
  expect_named(res_overall, c("estimates", "response_type", "ci_level", "trt_var"))
  expect_s3_class(res_overall$estimates, "tbl_df")
  expect_equal(nrow(res_overall$estimates), 1)
  expect_equal(res_overall$estimates$Subgroup, "Overall")
  expect_equal(res_overall$response_type, "continuous")
  expect_equal(res_overall$ci_level, 0.95)
  expect_equal(res_overall$trt_var, "trt")
})

test_that("summary_subgroup_effects works with specific subgroup_vars", {
  # sample_summary_obj was created with c("region", "sex")
  expect_s3_class(sample_summary_obj, "subgroup_summary")
  expect_s3_class(sample_summary_obj$estimates, "tbl_df")
  # Expect Overall + 3 region levels + 2 sex levels
  expect_equal(nrow(sample_summary_obj$estimates), 1 + 3 + 2)
  expect_true("Overall" %in% sample_summary_obj$estimates$Subgroup)
  expect_true(all(c("region: A", "region: B", "region: C") %in% sample_summary_obj$estimates$Subgroup))
  expect_true(all(c("sex: F", "sex: M") %in% sample_summary_obj$estimates$Subgroup))
})

test_that("summary_subgroup_effects works with subgroup_vars = 'auto'", {
  res_auto <- suppressMessages(
    summary_subgroup_effects(
      brms_fit = minimal_brms_fit_summary,
      original_data = original_test_data_summary,
      trt_var = "trt",
      response_type = "continuous",
      subgroup_vars = "auto"
    )
  )
  expect_s3_class(res_auto, "subgroup_summary")
  # Expect Overall + detected 'region' (3 levels) + detected 'sex' (2 levels)
  expect_equal(nrow(res_auto$estimates), 1 + 3 + 2)
  expect_true("Overall" %in% res_auto$estimates$Subgroup)
  expect_true(all(c("region: A", "region: B", "region: C") %in% res_auto$estimates$Subgroup))
  expect_true(all(c("sex: F", "sex: M") %in% res_auto$estimates$Subgroup))
})


test_that("summary_subgroup_effects assertions catch invalid inputs", {
  # Invalid brms_fit
  expect_error(
    summary_subgroup_effects(brms_fit = list(), original_data = original_test_data_summary, trt_var = "trt", response_type = "continuous"),
    regexp = "Must inherit from class 'brmsfit'"
  )
  # Invalid original_data
  expect_error(
    summary_subgroup_effects(brms_fit = minimal_brms_fit_summary, original_data = as.matrix(original_test_data_summary), trt_var = "trt", response_type = "continuous"),
    regexp = "Must be of type 'data.frame'"
  )
  # trt_var not in original_data
  expect_error(
    summary_subgroup_effects(brms_fit = minimal_brms_fit_summary, original_data = original_test_data_summary, trt_var = "treatment", response_type = "continuous"),
    regexp = "Must be a subset of"
  )
  # Invalid subgroup_vars type
  expect_error(
    summary_subgroup_effects(brms_fit = minimal_brms_fit_summary, original_data = original_test_data_summary, trt_var = "trt", subgroup_vars = 123, response_type = "continuous"),
    regexp = "Assertion failed"
  )
  # subgroup_vars not in original_data
  expect_error(
    summary_subgroup_effects(brms_fit = minimal_brms_fit_summary, original_data = original_test_data_summary, trt_var = "trt", subgroup_vars = c("region", "missing_var"), response_type = "continuous"),
    regexp = "Must be a subset of"
  )
})

# --- Tests for plot.subgroup_summary ---

test_that("plot.subgroup_summary runs and returns ggplot object", {
  p <- plot(sample_summary_obj)
  expect_s3_class(p, "ggplot")

  # Test providing labels
  p_labelled <- plot(sample_summary_obj, x_lab = "Custom X", title = "Custom Title")
  expect_s3_class(p_labelled, "ggplot")
  expect_equal(p_labelled$labels$x, "Custom X")
  expect_equal(p_labelled$labels$title, "Custom Title")
})

test_that("plot.subgroup_summary assertions catch invalid inputs", {
  # Invalid input class
  expect_error(
    plot(list(estimates = data.frame())), # Not a subgroup_summary object
    regexp = "does not have components 'x' and 'y'" # <-- Update this pattern
  )

  # Empty estimates table
  bad_summary_obj <- sample_summary_obj
  bad_summary_obj$estimates <- bad_summary_obj$estimates[0, ] # Make it empty
  expect_error(
    plot(bad_summary_obj),
    regexp = "Must have at least 1 rows"
  )
})
