# tests/testthat/test-fit_brms_model.R

# --- Setup: Prepare sample data and formula ---
# (We need this setup before the tests run)

# Load brms explicitly for helper functions if needed in tests
# library(brms)

# Sample data
test_data_fit <- data.frame(
  outcome = rnorm(20),
  time = round(runif(20, 1, 50)),
  status = sample(0:1, 20, replace = TRUE),
  trt = sample(0:1, 20, replace = TRUE),
  age = rnorm(20, 50, 10),
  region = factor(sample(c("A", "B"), 20, replace = TRUE))
)
# Ensure trt is factor initially for prepare_formula_model
test_data_fit$trt <- factor(test_data_fit$trt, levels=c(0, 1))

# Prepare a basic formula object (continuous)
prep_cont <- prepare_formula_model(
  data = test_data_fit,
  response_formula_str = "outcome ~ trt",
  response_type = "continuous",
  shrunk_prognostic_formula_str = "~ age",
  shrunk_predictive_formula_str = "~ trt:region"
)

# Prepare a basic formula object (survival) - needed for some prior tests
prep_surv <- prepare_formula_model(
  data = test_data_fit,
  response_formula_str = "Surv(time, status) ~ trt",
  response_type = "survival",
  unshrunk_prognostic_formula_str = "~ age"
)


# --- Helper to run minimal brms models ---
run_quick_brm <- function(...) {
  # Suppress messages and use minimal settings for speed
  suppressMessages(
    fit_brms_model(..., chains = 1, iter = 10, warmup = 5, refresh = 0,
                   # Use a minimal backend if available (or keep default)
                   backend = "cmdstanr", # Faster startup usually
                   cores=1)
  )
}

# --- Tests Start Here ---

test_that("Basic model fitting works for different response types", {
  # Test continuous model
  fit_cont <- run_quick_brm(
    formula = prep_cont$formula,
    data = prep_cont$data,
    response_type = "continuous"
  )
  expect_s3_class(fit_cont, "brmsfit")

  # Test survival model
  fit_surv <- run_quick_brm(
    formula = prep_surv$formula,
    data = prep_surv$data,
    response_type = "survival"
  )
  expect_s3_class(fit_surv, "brmsfit")

  # Add tests for "binary" and "count" if needed, preparing data/formula first
})

test_that("Prior assignment works correctly", {
  # 1. Default priors (empty lists passed)
  fit_default <- run_quick_brm(
    formula = prep_cont$formula,
    data = prep_cont$data,
    response_type = "continuous",
    predictive_effect_priors = list(), # Explicitly empty
    prognostic_effect_priors = list()  # Explicitly empty
  )
  expect_s3_class(fit_default, "brmsfit")

  # --- Revised Prior Check using the prior data frame ---
  default_priors_df <- as.data.frame(fit_default$prior) # Get priors as a data frame

  # Function to get prior for a specific nlpar (REVISED)
  get_nlpar_prior <- function(df, nlpar_name) {
    # Explicitly find rows matching the nlpar
    rows <- which(df$nlpar == nlpar_name)
    # If no rows found, return NA
    if(length(rows) == 0) return(NA_character_)
    # Get the prior value from the *first* matching row
    prior_val <- df[rows[1], "prior"]
    # Now prior_val is guaranteed to be length 1 (or NA if the cell was empty)
    # Check if the single value is NA or empty string
    if(is.na(prior_val) || prior_val == "") return(NA_character_)
    return(prior_val)
  }

  # Check shprogeffect default
  expect_match(get_nlpar_prior(default_priors_df, "shprogeffect"), "horseshoe\\(1\\)")

  # Check shpredeffect default
  expect_match(get_nlpar_prior(default_priors_df, "shpredeffect"), "horseshoe\\(1\\)")

  # Check unprogeffect default (This was the failing one)
  expect_match(get_nlpar_prior(default_priors_df, "unprogeffect"), "normal\\(0,\\s*10\\)")


  # 2. User-specified priors (as strings)
  fit_user_string <- run_quick_brm(
    formula = prep_cont$formula,
    data = prep_cont$data,
    response_type = "continuous",
    predictive_effect_priors = list(shrunk = "normal(0, 1)"), # Custom shrunk prior
    prognostic_effect_priors = list(shrunk = "normal(0, 2)", unshrunk = "normal(0, 3)")
  )
  expect_s3_class(fit_user_string, "brmsfit")

  user_priors_df <- as.data.frame(fit_user_string$prior)

  # Check shrunk predictive prior
  expect_match(get_nlpar_prior(user_priors_df, "shpredeffect"), "normal\\(0,\\s*1\\)")
  # Check shrunk prognostic prior
  expect_match(get_nlpar_prior(user_priors_df, "shprogeffect"), "normal\\(0,\\s*2\\)")
  # Check unshrunk prognostic prior
  expect_match(get_nlpar_prior(user_priors_df, "unprogeffect"), "normal\\(0,\\s*3\\)")


  # 3. User-specified priors (as brmsprior object - re-targeting check)
  complex_prior <- "R2D2(mean_R2 = 0.5, prec_R2 = 1)"
  fit_user_brmsprior <- run_quick_brm(
    formula = prep_cont$formula,
    data = prep_cont$data,
    response_type = "continuous",
    prognostic_effect_priors = list(shrunk = complex_prior) # Apply complex prior to shrunk prog
  )
  expect_s3_class(fit_user_brmsprior, "brmsfit")
  prior_summary_brmsprior <- capture.output(print(fit_user_brmsprior$prior))
  # --- Revised Prior Check using the prior data frame ---
  brmsprior_priors_df <- as.data.frame(fit_user_brmsprior$prior)

  # Check if the R2D2 prior was correctly applied to the nlpar
  expect_match(get_nlpar_prior(brmsprior_priors_df, "shprogeffect"), "R2D2\\(") # Check it starts with R2D2(

  # Check that the default horseshoe wasn't used for shprogeffect
  expect_false(grepl("horseshoe", get_nlpar_prior(brmsprior_priors_df, "shprogeffect")))
})

test_that("Stanvars argument is accepted", {
  # Create a dummy stanvar
  sv <- brms::stanvar(scode = "real dummy_var;", block = "parameters")
  fit_stanvar <- run_quick_brm(
    formula = prep_cont$formula,
    data = prep_cont$data,
    response_type = "continuous",
    stanvars = sv
  )
  expect_s3_class(fit_stanvar, "brmsfit")
  # Check if stanvar code is included (requires looking at generated Stan code)
  # This is a basic check; deeper checks are complex.
  expect_true(grepl("dummy_var", fit_stanvar$model))
})


test_that("Assertions work for invalid inputs", {
  # Invalid formula class
  expect_error(
    fit_brms_model(formula = lm(outcome ~ 1, data=test_data_fit), data = prep_cont$data, response_type = "continuous"),
    regexp = "Must inherit from class 'brmsformula'"
  )
  # Invalid data
  expect_error(
    fit_brms_model(formula = prep_cont$formula, data = as.matrix(prep_cont$data), response_type = "continuous"),
    regexp = "Must be of type 'data.frame'" # <-- Update this pattern
  )
  # Invalid response_type
  expect_error(
    fit_brms_model(formula = prep_cont$formula, data = prep_cont$data, response_type = "gaussian"),
    regexp = "Must be element of set"
  )
  # Invalid prior list (not named)
  # Around line 183
  expect_error(
    fit_brms_model(formula = prep_cont$formula, data = prep_cont$data, response_type = "continuous",
                   prognostic_effect_priors = list("normal(0,1)")), # Unnamed list
    regexp = "Must have names" # <-- Update this pattern
  )
  # Invalid stanvars class
  expect_error(
    fit_brms_model(formula = prep_cont$formula, data = prep_cont$data, response_type = "continuous",
                   stanvars = list("dummy")), # Not a stanvars object
    regexp = "Must inherit from class 'stanvars'" # <-- Update this pattern
  )
})
