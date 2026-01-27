# tests/testthat/test-prepare_formula_model.R

# Helper function to extract formula components for easier testing
get_formula_rhs <- function(bf_obj, nlpar = NULL) {
  if (is.null(nlpar)) {
    form_str_vec <- deparse(bf_obj$formula)
  } else {
    if (!nlpar %in% names(bf_obj$pforms)) return(character(0))
    form_str_vec <- deparse(bf_obj$pforms[[nlpar]])
  }
  form_str <- paste(form_str_vec, collapse = "")
  rhs <- gsub("^.*~\\s*", "", form_str)
  if (rhs == "" || rhs == "0") return(character(0))
  terms <- trimws(strsplit(rhs, "\\+")[[1]])

  # **FIX: Remove both "1" and "0" terms**
  terms <- terms[!terms %in% c("1", "0")]

  terms <- terms[nzchar(terms)]
  return(sort(terms))
}

# --- Sample Data for Testing ---
test_data <- data.frame(
  outcome = rnorm(10),
  n_events = rpois(10, 5),
  log_days = log(runif(10, 50, 150)),
  time = round(runif(10, 1, 100)),
  status = sample(0:1, 10, replace = TRUE),
  trt = sample(c("Control", "Treatment"), 10, replace = TRUE),
  age = rnorm(10, 50, 10),
  region = factor(sample(c("A", "B"), 10, replace = TRUE)),
  subgroup1 = factor(sample(c("S1", "S2"), 10, replace = TRUE)),
  subgroup2 = factor(sample(c("X", "Y"), 10, replace = TRUE))
)
test_data$trt <- factor(test_data$trt) # Make sure it's a factor initially

# --- Tests Start Here ---

test_that("Basic functionality works (binary/continuous)", {
  # --- Binary example ---
  res_bin <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "binary",
    unshrunk_terms_formula = "~ age"  # Updated parameter name
  )

  # Check output structure
  expect_s3_class(res_bin$formula, "brmsformula")
  expect_s3_class(res_bin$data, "data.frame")
  expect_named(res_bin, c("formula", "data", "response_type"))
  expect_equal(res_bin$response_type, "binary")

  # Check data modifications - treatment should be numeric binary (0/1)
  expect_true(is.numeric(res_bin$data$trt))
  expect_true(all(res_bin$data$trt %in% c(0, 1)))

  # Check formula components (trt automatically added to unshrunktermeffect)
  expect_equal(get_formula_rhs(res_bin$formula), "unshrunktermeffect") # Main formula has placeholder
  expect_equal(get_formula_rhs(res_bin$formula, "unshrunktermeffect"), c("age", "trt"))

  # --- Continuous example ---
  res_cont <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ age" # Put age in shrunk this time
  )
  expect_equal(res_cont$response_type, "continuous")
  expect_equal(get_formula_rhs(res_cont$formula), c("shprogeffect", "unshrunktermeffect"))
  expect_equal(get_formula_rhs(res_cont$formula, "unshrunktermeffect"), "trt")
  expect_equal(get_formula_rhs(res_cont$formula, "shprogeffect"), "age")
})

# --- Test Interaction Handling ---
test_that("Predictive terms (interactions) are processed correctly", {
  res_int <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ trt:subgroup1",
    unshrunk_terms_formula = "~ trt:region"  # Updated: now goes in unshrunk_terms
  )

  # Check that subgroup variables are factors with contrasts
  expect_true(is.factor(res_int$data$subgroup1))
  expect_true(is.factor(res_int$data$region))
  expect_true(!is.null(contrasts(res_int$data$subgroup1)))
  expect_true(!is.null(contrasts(res_int$data$region)))

  # Check formula components - all unshrunk terms consolidated into unshrunktermeffect
  expect_equal(get_formula_rhs(res_int$formula), c("shpredeffect", "unshrunktermeffect"))
  expect_equal(get_formula_rhs(res_int$formula, "shpredeffect"), "trt:subgroup1")
  # Both trt and trt:region go into unshrunktermeffect
  expect_equal(get_formula_rhs(res_int$formula, "unshrunktermeffect"), c("trt", "trt:region"))

  # Check returned response type
  expect_equal(res_int$response_type, "continuous")
})

# --- Test Overlaps and Defaults ---
test_that("Term overlaps and defaults are handled", {
  # Expect warning when terms overlap between unshrunk and shrunk
  expect_warning(
    res_overlap_prog <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      unshrunk_terms_formula = "~ age + region",  # Updated parameter name
      shrunk_prognostic_formula = "~ age + subgroup1" # age overlaps
    ),
    regexp = "Prioritizing as unshrunk: age"
  )
  expect_equal(res_overlap_prog$response_type, "continuous")

  # Check that 'age' ended up in unshrunk (all unshrunk terms consolidated)
  expect_equal(get_formula_rhs(res_overlap_prog$formula, "unshrunktermeffect"), c("age", "region", "trt"))
  expect_equal(get_formula_rhs(res_overlap_prog$formula, "shprogeffect"), "subgroup1")

  # Test overlap between unshrunk and shrunk predictive
  # Note: Now we don't have separate unshrunk_predictive_formula
  # All unshrunk terms go in unshrunk_terms_formula
  expect_warning(
    res_overlap_pred <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      unshrunk_terms_formula = "~ region",  # Main effect in unshrunk
      shrunk_predictive_formula = "~ trt:region"  # Interaction in shrunk
    ),
    regexp = NA  # No overlap warning expected - different terms
  )
  expect_equal(res_overlap_pred$response_type, "continuous")

  # Check formula components
  expect_true("region" %in% get_formula_rhs(res_overlap_pred$formula, "unshrunktermeffect"))
  expect_equal(get_formula_rhs(res_overlap_pred$formula, "shpredeffect"), "trt:region")
})

# --- Test Different Response Types ---
test_that("Count models handle offset", {
  res_count <- prepare_formula_model(
    data = test_data,
    response_formula = "n_events + offset(log_days) ~ trt",
    response_type = "count"
  )
  # Check output structure
  expect_named(res_count, c("formula", "data", "response_type"))
  expect_equal(res_count$response_type, "count")

  # Check that offset is included correctly in the main formula string
  main_form_str <- deparse(res_count$formula$formula)
  expect_true(grepl("n_events ~ unshrunktermeffect \\+ log_days", main_form_str))
})

test_that("Stratification works for continuous and count", {
  # Continuous - sigma stratified
  res_cont_strat <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    stratification_formula = "~ region"
  )
  expect_equal(res_cont_strat$response_type, "continuous")
  expect_true("sigma" %in% names(res_cont_strat$formula$pforms))
  expect_true(grepl("sigma ~ region", deparse(res_cont_strat$formula$pforms$sigma)))

  # Count - shape stratified
  res_count_strat <- prepare_formula_model(
    data = test_data,
    response_formula = "n_events ~ trt",
    response_type = "count",
    stratification_formula = "~ region"
  )
  expect_equal(res_count_strat$response_type, "count")
  expect_true("shape" %in% names(res_count_strat$formula$pforms))
  expect_true(grepl("shape ~ region", deparse(res_count_strat$formula$pforms$shape)))
})

test_that("Star notation is supported and preserved", {
  res_star <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age + trt*region"  # Star notation
  )

  # Check that region has contrasts applied
  expect_true(is.factor(res_star$data$region))
  expect_true(!is.null(contrasts(res_star$data$region)))

  # Check that star notation is preserved in the formula
  unshrunk_terms <- get_formula_rhs(res_star$formula, "unshrunktermeffect")
  # Should contain age, trt, and trt*region (star notation preserved)
  expect_true("age" %in% unshrunk_terms)
  expect_true("trt" %in% unshrunk_terms)
  # Star notation should be in the formula
  expect_true(any(grepl("\\*", unshrunk_terms)) || "trt:region" %in% unshrunk_terms)
})

test_that("Star notation excludes from marginality check", {
  # Star notation includes main effects, so no warning expected
  expect_message(
    res_star_no_warn <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      unshrunk_terms_formula = "~ trt*region",  # Includes main effects
      shrunk_predictive_formula = "~ trt:subgroup1"  # Missing subgroup1 main effect
    ),
    regexp = "Marginality principle not followed.*subgroup1",
    all = TRUE
  )
  
  # Should NOT get warning about region (it's in star notation)
  # Should get warning about subgroup1 (it's in colon notation)
  expect_equal(res_star_no_warn$response_type, "continuous")
})

test_that("Survival models work (basic, stratified, fallback)", {
  # Basic survival
  res_surv <- prepare_formula_model(
    data = test_data,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival"
  )
  # Check output structure
  expect_named(res_surv, c("formula", "data", "response_type"))
  expect_equal(res_surv$response_type, "survival")

  main_form_str <- paste(deparse(res_surv$formula$formula), collapse = "")
  expect_true(grepl("time \\| cens\\(1 - status\\)", main_form_str))
  expect_true(grepl("\\+ bhaz\\(", main_form_str)) # bhaz term added
  expect_false(grepl("gr =", main_form_str)) # No stratification

  # Stratified survival
  res_surv_strat <- prepare_formula_model(
    data = test_data,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    stratification_formula = "~ region"
  )
  expect_equal(res_surv_strat$response_type, "survival")

  main_form_str_strat <- paste(deparse(res_surv_strat$formula$formula), collapse = "")
  expect_true(grepl("\\+ bhaz\\(", main_form_str_strat))
  expect_true(grepl("gr = region", main_form_str_strat)) # Stratification added

  # Knot calculation fallback (create data with few unique times)
  sparse_data <- test_data[1:5, ]
  sparse_data$time <- c(10, 10, 20, 20, 30) # Only 3 unique times
  expect_warning(
    res_surv_sparse <- prepare_formula_model(
      data = sparse_data,
      response_formula = "Surv(time, status) ~ trt",
      response_type = "survival"
    ),
    regexp = "Not enough unique event times to compute quantile knots"
  )
  expect_equal(res_surv_sparse$response_type, "survival")

  # Check that it fell back to Cox model (no bhaz term)
  main_form_str_sparse <- deparse(res_surv_sparse$formula$formula)
  expect_false(grepl("\\+ bhaz\\(", main_form_str_sparse))
})


# --- Test Error Conditions ---
test_that("Assertions and error checks work", {
  # This entire block tests for errors, so the function never returns.
  # No changes are needed here.

  # Invalid response type
  expect_error(
    prepare_formula_model(test_data, "outcome ~ trt", response_type = "wrong"),
    regexp = "Must be element of set"
  )
  # Formula string without ~
  expect_error(
    prepare_formula_model(test_data, "outcome", response_type = "continuous"),
    regexp = "Must comply to pattern '~'"
  )
  # trt_var not in data
  expect_error(
    prepare_formula_model(test_data, "outcome ~ treatment_var", response_type = "continuous"),
    regexp = "subset of"
  )
  # Strat var not in data (survival)
  expect_error(
    prepare_formula_model(test_data, "Surv(time, status) ~ trt", response_type = "survival", stratification_formula = "~ missing_var"),
    regexp = "Must be a subset of"
  )
  # Strat var not in data (continuous)
  expect_error(
    prepare_formula_model(test_data, "outcome ~ trt", response_type = "continuous", stratification_formula = "~ missing_var"),
    regexp = "Must be a subset of"
  )
  # Survival formula doesn't start with Surv
  expect_error(
    prepare_formula_model(test_data, "time ~ trt", response_type = "survival"),
    regexp = "Must comply to pattern"
  )
})
