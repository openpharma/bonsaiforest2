# tests/testthat/test-run_brms_analysis.R

# --- Sample Data for Testing ---
# (Using a smaller dataset for quicker setup)
test_data_run <- data.frame(
  outcome = rnorm(20),
  time = round(runif(20, 1, 50)),
  status = sample(0:1, 20, replace = TRUE),
  trt = sample(0:1, 20, replace = TRUE),
  age = rnorm(20, 50, 10),
  region = factor(sample(c("A", "B"), 20, replace = TRUE))
)
test_data_run$trt <- factor(test_data_run$trt, levels=c(0, 1))

# --- Tests Start Here ---

test_that("run_brms_analysis executes and returns brmsfit object", {
  # Suppress messages and use minimal settings for speed
  fit <- suppressMessages(
    run_brms_analysis(
      data = test_data_run,
      response_formula_str = "outcome ~ trt",
      response_type = "continuous",
      shrunk_prognostic_formula_str = "~ age",
      shrunk_predictive_formula_str = "~ trt:region",
      prognostic_effect_priors = list(shrunk = "normal(0,1)"),
      predictive_effect_priors = list(shrunk = "normal(0,1)"),
      chains = 1, iter = 10, warmup = 5, refresh = 0,
      backend = "cmdstanr", # Use faster backend if available
      cores = 1
    )
  )

  # Basic check: Did it return the right object class?
  expect_s3_class(fit, "brmsfit")

  # Optional: Deeper check - verify a prior was set as expected
  priors_df <- as.data.frame(fit$prior)

  # Find rows matching the nlpar='shprogeffect'
  shrunk_prog_rows <- priors_df$nlpar == "shprogeffect"

  # Check that *at least one* row was found for this nlpar
  expect_true(any(shrunk_prog_rows),
              info = "Expected at least one prior row for nlpar='shprogeffect'")

  # Extract all prior strings for this nlpar
  shrunk_prog_priors_found <- priors_df[shrunk_prog_rows, "prior"]

  # Check if the expected prior string exists anywhere in the found priors
  expect_true(any(grepl("normal\\(0,\\s*1\\)", shrunk_prog_priors_found)),
              info = "Expected prior 'normal(0,1)' not found among priors for nlpar='shprogeffect'")

})

test_that("run_brms_analysis assertions catch invalid inputs", {

  # Invalid data (matrix)
  expect_error(
    run_brms_analysis(data = as.matrix(test_data_run), response_formula_str = "outcome ~ trt", response_type = "continuous"),
    regexp = "Must be of type 'data.frame'"
  )

  # Invalid response_formula_str (missing ~)
  expect_error(
    run_brms_analysis(data = test_data_run, response_formula_str = "outcome", response_type = "continuous"),
    regexp = "Must comply to pattern '~'"
  )

  # Invalid response_type
  expect_error(
    run_brms_analysis(data = test_data_run, response_formula_str = "outcome ~ trt", response_type = "gaussian"),
    regexp = "Must be element of set"
  )

  # Invalid optional formula string (doesn't start with ~)
  expect_error(
    run_brms_analysis(data = test_data_run, response_formula_str = "outcome ~ trt", response_type = "continuous",
                      shrunk_prognostic_formula_str = "age"), # Missing ~
    regexp = "Must comply to pattern"
  )

  # Invalid prior list (not named)
  expect_error(
    run_brms_analysis(data = test_data_run, response_formula_str = "outcome ~ trt", response_type = "continuous",
                      prognostic_effect_priors = list("normal(0,1)")), # Unnamed list
    regexp = "Must have names" # <-- Update this pattern
  )
  # Invalid stanvars class
  expect_error(
    run_brms_analysis(data = test_data_run, response_formula_str = "outcome ~ trt", response_type = "continuous",
                      stanvars = list("dummy")), # Not a stanvars object
    regexp = "Must inherit from class 'stanvars'"
  )

})
