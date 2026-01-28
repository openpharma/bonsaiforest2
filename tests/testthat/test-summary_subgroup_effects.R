# tests/testthat/test-summary_subgroup_effects.R

# --- Test Setup ---

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
minimal_brms_fit_summary <- suppressMessages(
  run_brms_analysis(
    data = original_test_data_summary,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age",  # Use ~ 0 + to avoid intercept warning
    shrunk_predictive_formula = "~ 0 + trt:region + trt:sex",  # Use ~ 0 + to avoid intercept warning
    sigma_ref = sd(original_test_data_summary$outcome),
    chains = 1, iter = 10, warmup = 5, refresh = 0,
    backend = "cmdstanr", cores = 1
  )
)

# 3. Create a sample summary object ONCE for plot tests
sample_summary_obj <- suppressMessages(
  summary_subgroup_effects(
    brms_fit = minimal_brms_fit_summary,
    subgroup_vars = c("region", "sex")  # Other params auto-extracted
  )
)

# --- Tests for summary_subgroup_effects ---

test_that("summary_subgroup_effects throws error if subgroup_vars = NULL", {
  # We expect an error because NULL is no longer allowed.
  # We use a broader regex "Assertion" to match both:
  # "Assertion on 'subgroup_vars' failed" AND "Assertion failed."
  expect_error(
    summary_subgroup_effects(
      brms_fit = minimal_brms_fit_summary,
      subgroup_vars = NULL  # Other params auto-extracted
    ),
    regexp = "Assertion.*failed"
  )
})

test_that("summary_subgroup_effects works with specific subgroup_vars", {
  # sample_summary_obj was created with c("region", "sex")
  expect_s3_class(sample_summary_obj, "subgroup_summary")
  expect_s3_class(sample_summary_obj$estimates, "tbl_df")

  # Expect 3 region levels + 2 sex levels (Overall is REMOVED)
  expect_equal(nrow(sample_summary_obj$estimates), 3 + 2)

  # Verify Overall is NOT there
  expect_false("Overall" %in% sample_summary_obj$estimates$Subgroup)

  # Verify subgroups are there
  expect_true(all(c("region: A", "region: B", "region: C") %in% sample_summary_obj$estimates$Subgroup))
  expect_true(all(c("sex: F", "sex: M") %in% sample_summary_obj$estimates$Subgroup))

  # Check metadata
  expect_equal(sample_summary_obj$response_type, "continuous")
  expect_equal(sample_summary_obj$ci_level, 0.95)
  expect_equal(sample_summary_obj$trt_var, "trt")
})

test_that("summary_subgroup_effects works with subgroup_vars = 'auto'", {
  res_auto <- suppressMessages(
    summary_subgroup_effects(
      brms_fit = minimal_brms_fit_summary,
      subgroup_vars = "auto"  # Other params auto-extracted
    )
  )
  expect_s3_class(res_auto, "subgroup_summary")

  # Expect detected 'region' (3 levels) + detected 'sex' (2 levels)
  # NO Overall row
  expect_equal(nrow(res_auto$estimates), 3 + 2)
  expect_false("Overall" %in% res_auto$estimates$Subgroup)

  expect_true(all(c("region: A", "region: B", "region: C") %in% res_auto$estimates$Subgroup))
  expect_true(all(c("sex: F", "sex: M") %in% res_auto$estimates$Subgroup))
})


test_that("summary_subgroup_effects assertions catch invalid inputs", {
  # Invalid brms_fit
  expect_error(
    summary_subgroup_effects(brms_fit = list()),
    regexp = "Must inherit from class 'brmsfit'"
  )
  
  # trt_var not in data when explicitly provided
  expect_error(
    summary_subgroup_effects(brms_fit = minimal_brms_fit_summary, trt_var = "treatment"),
    regexp = "Must be a subset of"
  )
  
  # Invalid subgroup_vars type (number instead of char or 'auto')
  expect_error(
    summary_subgroup_effects(brms_fit = minimal_brms_fit_summary, subgroup_vars = 123),
    regexp = "Assertion failed"
  )
  
  # subgroup_vars not in data
  expect_error(
    summary_subgroup_effects(brms_fit = minimal_brms_fit_summary, subgroup_vars = c("region", "missing_var")),
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
  # To test the assertion INSIDE plot.subgroup_summary, we must pass an object
  # that has the class but fails internal checks, OR call the function directly.
  # Standard generic plot() on a list calls base::plot, so we use explicit call here:

  expect_error(
    plot.subgroup_summary(list(estimates = data.frame())),
    regexp = "Must inherit from class 'subgroup_summary'"
  )

  # Empty estimates table
  bad_summary_obj <- sample_summary_obj
  bad_summary_obj$estimates <- bad_summary_obj$estimates[0, ] # Make it empty

  expect_error(
    plot(bad_summary_obj),
    regexp = "Must have at least 1 rows"
  )
})
