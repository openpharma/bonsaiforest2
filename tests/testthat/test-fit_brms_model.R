# tests/testthat/test-fit_brms_model.R

# --- Setup: Prepare sample data and formula ---
# (We need this setup before the tests run)

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
test_data_fit$trt <- factor(test_data_fit$trt, levels = c(0, 1))

# Prepare a basic formula object (continuous)
# This now correctly returns a list with formula, data, and response_type
prep_cont <- prepare_formula_model(
  data = test_data_fit,
  response_formula_str = "outcome ~ trt",
  response_type = "continuous",
  shrunk_prognostic_formula_str = "~ age",
  shrunk_predictive_formula_str = "~ trt:region"
)

# Prepare a basic formula object (survival)
prep_surv <- prepare_formula_model(
  data = test_data_fit,
  response_formula_str = "Surv(time, status) ~ trt",
  response_type = "survival",
  unshrunk_prognostic_formula_str = "~ age"
)


# --- Helper to run minimal brms models ---
# This helper is perfect as-is, it passes `...` to fit_brms_model
run_quick_brm <- function(...) {
  # Suppress messages and use minimal settings for speed
  suppressMessages(
    fit_brms_model(..., chains = 1, iter = 1000, warmup = 50, refresh = 0,
                   # Use a minimal backend if available (or keep default)
                   backend = "cmdstanr", # Faster startup usually
                   cores = 1)
  )
}

# --- Tests Start Here ---

test_that("Basic model fitting works for different response types", {
  # Test continuous model
  fit_cont <- run_quick_brm(
    prepared_model = prep_cont
  )
  expect_s3_class(fit_cont, "brmsfit")

  # Test survival model

  fit_surv <- run_quick_brm(
    prepared_model = prep_surv
  )
  expect_s3_class(fit_surv, "brmsfit")

  # Add tests for "binary" and "count" if needed, preparing data/formula first
})

test_that("Prior assignment works correctly", {
  # 1. Default priors (empty lists passed)
  fit_default <- run_quick_brm(
    prepared_model = prep_cont,
    predictive_effect_priors = list(), # Explicitly empty
    prognostic_effect_priors = list()  # Explicitly empty
  )
  expect_s3_class(fit_default, "brmsfit")

  default_priors_df <- as.data.frame(fit_default$prior)

  # Helper function (no change needed)
  get_nlpar_prior <- function(df, nlpar_name) {
    rows <- which(df$nlpar == nlpar_name & df$class == "b" & df$coef == "")
    if (nlpar_name == "unprogeffect") {
      # Need to find the non-intercept 'b' prior
      rows <- which(df$nlpar == "unprogeffect" & df$class == "b" & df$coef == "")
    }
    if (length(rows) == 0) return(NA_character_)
    prior_val <- df[rows[1], "prior"]
    if (is.na(prior_val) || prior_val == "") return(NA_character_)
    return(prior_val)
  }

  # Check shprogeffect default
  expect_match(get_nlpar_prior(default_priors_df, "shprogeffect"), "horseshoe\\(1\\)")

  # Check shpredeffect default
  expect_match(get_nlpar_prior(default_priors_df, "shpredeffect"), "horseshoe\\(1\\)")

  # Check unprogeffect default (for class 'b' non-intercept)
  expect_match(get_nlpar_prior(default_priors_df, "unprogeffect"), "normal\\(0,\\s*5\\)")


  # 2. User-specified priors (as strings)
  fit_user_string <- run_quick_brm(
    prepared_model = prep_cont,
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
  # --- UPDATED CALL ---
  fit_user_brmsprior <- run_quick_brm(
    prepared_model = prep_cont,
    prognostic_effect_priors = list(shrunk = complex_prior) # Apply complex prior
  )
  expect_s3_class(fit_user_brmsprior, "brmsfit")

  brmsprior_priors_df <- as.data.frame(fit_user_brmsprior$prior)

  # Check if the R2D2 prior was correctly applied
  expect_match(get_nlpar_prior(brmsprior_priors_df, "shprogeffect"), "R2D2\\(")
  expect_false(grepl("horseshoe", get_nlpar_prior(brmsprior_priors_df, "shprogeffect")))
})

test_that("Stanvars argument is accepted", {
  # Create a dummy stanvar
  sv <- brms::stanvar(scode = "real dummy_var;", block = "parameters")

  # --- UPDATED CALL ---
  fit_stanvar <- run_quick_brm(
    prepared_model = prep_cont,
    stanvars = sv
  )
  expect_s3_class(fit_stanvar, "brmsfit")
  expect_true(grepl("dummy_var", fit_stanvar$model))
})


test_that("Assertions work for invalid inputs", {
  # This block calls fit_brms_model directly to test assertions

  # --- NEW: Test prepared_model object itself ---
  # Invalid type (not a list)
  expect_error(
    fit_brms_model(prepared_model = "not a list"),
    regexp = "Must be of type 'list'"
  )

  # Missing required names
  bad_list_1 <- list(formula = prep_cont$formula, data = prep_cont$data) # Missing response_type
  expect_error(
    fit_brms_model(prepared_model = bad_list_1),
    # --- THIS IS THE FIX ---
    # The error message starts with "Names must include" or contains "missing elements"
    regexp = "missing elements.*'response_type'"
  )

  # --- NEW: Test contents of prepared_model ---
  # Invalid formula class inside list
  bad_list_2 <- prep_cont
  bad_list_2$formula <- "outcome ~ trt" # Is a string, not brmsformula
  expect_error(
    fit_brms_model(prepared_model = bad_list_2),
    regexp = "Must inherit from class 'brmsformula'"
  )

  # Invalid data class inside list
  bad_list_3 <- prep_cont
  bad_list_3$data <- as.matrix(prep_cont$data) # Is a matrix
  expect_error(
    fit_brms_model(prepared_model = bad_list_3),
    regexp = "Must be of type 'data.frame'"
  )

  # Invalid response_type inside list
  bad_list_4 <- prep_cont
  bad_list_4$response_type <- "gaussian" # Not one of the allowed choices
  expect_error(
    fit_brms_model(prepared_model = bad_list_4),
    regexp = "Must be element of set"
  )

  # --- Test other arguments (priors, stanvars) ---
  # Invalid prior list (not named)
  expect_error(
    fit_brms_model(prepared_model = prep_cont,
                   prognostic_effect_priors = list("normal(0,1)")), # Unnamed list
    regexp = "Must have names"
  )

  # Invalid stanvars class
  expect_error(
    fit_brms_model(prepared_model = prep_cont,
                   stanvars = list("dummy")), # Not a stanvars object
    regexp = "Must inherit from class 'stanvars'"
  )
})
