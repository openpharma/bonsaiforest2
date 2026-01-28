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
  response_formula = "outcome ~ trt",
  response_type = "continuous",
  shrunk_prognostic_formula = "~ 0 + age",  # Use ~ 0 + to avoid intercept warning
  shrunk_predictive_formula = "~ 0 + trt:region"  # Use ~ 0 + to avoid intercept warning
)

# Prepare a basic formula object (survival)
prep_surv <- prepare_formula_model(
  data = test_data_fit,
  response_formula = "Surv(time, status) ~ trt",
  response_type = "survival",
  unshrunk_terms_formula = "~ age"
)


# --- Helper to run minimal brms models ---
# This helper adds sigma_ref and minimal settings for speed
run_quick_brm <- function(prepared_model, sigma_ref = 1, ...) {
  # Suppress messages and use minimal settings for speed
  suppressMessages(
    fit_brms_model(prepared_model = prepared_model, 
                   sigma_ref = sigma_ref,
                   ..., 
                   chains = 1, iter = 1000, warmup = 50, refresh = 0,
                   # Use a minimal backend if available (or keep default)
                   backend = "cmdstanr", # Faster startup usually
                   cores = 1)
  )
}

# --- Tests Start Here ---

test_that("Basic model fitting works for different response types", {
  # Test continuous model
  fit_cont <- run_quick_brm(
    prepared_model = prep_cont,
    sigma_ref = sd(test_data_fit$outcome)
  )
  expect_s3_class(fit_cont, "brmsfit")

  # Test survival model (use sigma_ref = 1 for survival)
  fit_surv <- run_quick_brm(
    prepared_model = prep_surv,
    sigma_ref = 1
  )
  expect_s3_class(fit_surv, "brmsfit")

  # Add tests for "binary" and "count" if needed, preparing data/formula first
})

test_that("Prior assignment works correctly", {
  sigma_ref_val <- sd(test_data_fit$outcome)
  
  # 1. Default priors (NULL for all optional priors)
  fit_default <- run_quick_brm(
    prepared_model = prep_cont,
    sigma_ref = sigma_ref_val
    # All priors use defaults (NULL)
  )
  expect_s3_class(fit_default, "brmsfit")

  default_priors_df <- as.data.frame(fit_default$prior)

  # Helper function (updated for new nlpar name)
  get_nlpar_prior <- function(df, nlpar_name) {
    rows <- which(df$nlpar == nlpar_name & df$class == "b" & df$coef == "")
    if (nlpar_name == "unshrunktermeffect") {
      # Need to find the non-intercept 'b' prior
      rows <- which(df$nlpar == "unshrunktermeffect" & df$class == "b" & df$coef == "")
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

  # Check unshrunktermeffect default (for class 'b' non-intercept)
  # Note: Prior is calculated as normal(0, 5 * sigma_ref) which varies based on data
  expect_match(get_nlpar_prior(default_priors_df, "unshrunktermeffect"), "normal\\(0,\\s*[0-9.]+")


  # 2. User-specified priors (as strings)
  fit_user_string <- run_quick_brm(
    prepared_model = prep_cont,
    sigma_ref = sigma_ref_val,
    shrunk_predictive_prior = "normal(0, 1)", # Custom shrunk predictive prior
    shrunk_prognostic_prior = "normal(0, 2)", # Custom shrunk prognostic prior
    unshrunk_prior = "normal(0, 3)" # Custom unshrunk prior
  )
  expect_s3_class(fit_user_string, "brmsfit")

  user_priors_df <- as.data.frame(fit_user_string$prior)

  # Check shrunk predictive prior
  expect_match(get_nlpar_prior(user_priors_df, "shpredeffect"), "normal\\(0,\\s*1\\)")
  # Check shrunk prognostic prior
  expect_match(get_nlpar_prior(user_priors_df, "shprogeffect"), "normal\\(0,\\s*2\\)")
  # Check unshrunk prior
  expect_match(get_nlpar_prior(user_priors_df, "unshrunktermeffect"), "normal\\(0,\\s*3\\)")


  # 3. User-specified priors (complex prior string)
  complex_prior <- "R2D2(mean_R2 = 0.5, prec_R2 = 1)"
  fit_user_brmsprior <- run_quick_brm(
    prepared_model = prep_cont,
    sigma_ref = sigma_ref_val,
    shrunk_prognostic_prior = complex_prior # Apply complex prior
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

  fit_stanvar <- run_quick_brm(
    prepared_model = prep_cont,
    sigma_ref = sd(test_data_fit$outcome),
    stanvars = sv
  )
  expect_s3_class(fit_stanvar, "brmsfit")
  expect_true(grepl("dummy_var", fit_stanvar$model))
})


test_that("Assertions work for invalid inputs", {
  # This block calls fit_brms_model directly to test assertions

  # --- Test prepared_model object itself ---
  # Invalid type (not a list)
  expect_error(
    fit_brms_model(prepared_model = "not a list", sigma_ref = 1),
    regexp = "Must be of type 'list'"
  )

  # Missing required names
  bad_list_1 <- list(formula = prep_cont$formula, data = prep_cont$data) # Missing response_type
  expect_error(
    fit_brms_model(prepared_model = bad_list_1, sigma_ref = 1),
    regexp = "missing elements.*'response_type'"
  )

  # --- Test contents of prepared_model ---
  # Invalid formula class inside list
  bad_list_2 <- prep_cont
  bad_list_2$formula <- "outcome ~ trt" # Is a string, not brmsformula
  expect_error(
    fit_brms_model(prepared_model = bad_list_2, sigma_ref = 1),
    regexp = "Must inherit from class 'brmsformula'"
  )

  # Invalid data class inside list
  bad_list_3 <- prep_cont
  bad_list_3$data <- as.matrix(prep_cont$data) # Is a matrix
  expect_error(
    fit_brms_model(prepared_model = bad_list_3, sigma_ref = 1),
    regexp = "Must be of type 'data.frame'"
  )

  # Invalid response_type inside list
  bad_list_4 <- prep_cont
  bad_list_4$response_type <- "gaussian" # Not one of the allowed choices
  expect_error(
    fit_brms_model(prepared_model = bad_list_4, sigma_ref = 1),
    regexp = "Must be element of set"
  )

  # --- Test sigma_ref validation ---
  # Missing sigma_ref
  expect_error(
    fit_brms_model(prepared_model = prep_cont),
    regexp = "Must be of type 'number'"
  )
  
  # Invalid sigma_ref (negative)
  expect_error(
    fit_brms_model(prepared_model = prep_cont, sigma_ref = -1),
    regexp = "Element 1 is not >= 0"
  )
  
  # Invalid sigma_ref (not numeric)
  expect_error(
    fit_brms_model(prepared_model = prep_cont, sigma_ref = "text"),
    regexp = "Must be of type 'number'"
  )

  # --- Test other arguments (stanvars) ---
  # Invalid stanvars class
  expect_error(
    fit_brms_model(prepared_model = prep_cont,
                   sigma_ref = 1,
                   stanvars = list("dummy")), # Not a stanvars object
    regexp = "Must inherit from class 'stanvars'"
  )
})
