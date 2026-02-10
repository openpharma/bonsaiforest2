# tests/testthat/test-prepare_formula_model.R


# HELPER FUNCTIONS


# Extract formula components for easier testing
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
  terms <- terms[!terms %in% c("1", "0")]
  terms <- terms[nzchar(terms)]
  return(sort(terms))
}

# Check standard output structure
expect_valid_output <- function(result, expected_response_type, expected_trt_var = "trt") {
  expect_named(result, c("formula", "data", "response_type", "trt_var",
                         "stan_variable_names", "has_intercept", "has_random_effects"))
  expect_s3_class(result$formula, "brmsformula")
  expect_s3_class(result$data, "data.frame")
  expect_equal(result$response_type, expected_response_type)
  expect_equal(result$trt_var, expected_trt_var)
  expect_type(result$stan_variable_names, "list")
  expect_type(result$has_intercept, "logical")
  expect_type(result$has_random_effects, "logical")
}

# Check that variables are factors with contrasts
expect_factor_with_contrasts <- function(data, var_names) {
  for (var in var_names) {
    expect_true(is.factor(data[[var]]),
                info = paste(var, "should be a factor"))
    expect_true(!is.null(contrasts(data[[var]])),
                info = paste(var, "should have contrasts"))
  }
}

# Check contrast type (one-hot vs dummy)
expect_onehot_contrasts <- function(data, var_name) {
  contrasts_mat <- contrasts(data[[var_name]])
  expect_equal(ncol(contrasts_mat), nlevels(data[[var_name]]),
               info = paste(var_name, "should have one-hot encoding (k columns)"))
}

expect_dummy_contrasts <- function(data, var_name) {
  contrasts_mat <- contrasts(data[[var_name]])
  expect_equal(ncol(contrasts_mat), nlevels(data[[var_name]]) - 1,
               info = paste(var_name, "should have dummy encoding (k-1 columns)"))
}

# Check treatment conversion to binary
expect_binary_treatment <- function(data, trt_var = "trt") {
  expect_true(is.numeric(data[[trt_var]]))
  expect_true(all(data[[trt_var]] %in% c(0, 1)))
}


# TEST DATA


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
test_data$trt <- factor(test_data$trt)



# PART 1: BASIC FUNCTIONALITY


test_that("Binary response type works", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "binary",
    unshrunk_terms_formula = "~ age"
  )

  expect_valid_output(res, "binary")
  expect_binary_treatment(res$data)
  expect_equal(get_formula_rhs(res$formula), "unshrunktermeffect")
  expect_equal(get_formula_rhs(res$formula, "unshrunktermeffect"), c("age", "trt"))
})

test_that("Continuous response type works", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_prognostic_formula = "~ 0 + age"
  )

  expect_valid_output(res, "continuous")
  expect_equal(get_formula_rhs(res$formula), c("shprogeffect", "unshrunktermeffect"))
  expect_equal(get_formula_rhs(res$formula, "unshrunktermeffect"), "trt")
  expect_equal(get_formula_rhs(res$formula, "shprogeffect"), "age")
})

test_that("Count models handle offset", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "n_events + offset(log_days) ~ trt",
    response_type = "count"
  )

  expect_valid_output(res, "count")
  main_form_str <- deparse(res$formula$formula)
  expect_true(grepl("n_events ~ unshrunktermeffect \\+ log_days", main_form_str))
})

test_that("Treatment automatically added to unshrunk when missing", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age + region"
  )

  expect_true("trt" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
  expect_true("age" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
  expect_true("region" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
})

test_that("Minimal model (treatment only) works", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous"
  )

  expect_equal(get_formula_rhs(res$formula, "unshrunktermeffect"), "trt")
  expect_equal(length(get_formula_rhs(res$formula, "shprogeffect")), 0)
  expect_equal(length(get_formula_rhs(res$formula, "shpredeffect")), 0)
})



# PART 2: SURVIVAL MODELS


test_that("Basic survival model works", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival"
  )

  expect_valid_output(res, "survival")
  main_form_str <- paste(deparse(res$formula$formula), collapse = "")
  expect_true(grepl("time \\| cens\\(1 - status\\)", main_form_str))
  expect_true(grepl("\\+ bhaz\\(", main_form_str))
  expect_false(grepl("gr =", main_form_str))
})

test_that("Stratified survival model works", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    stratification_formula = "~ region"
  )

  expect_valid_output(res, "survival")
  main_form_str <- paste(deparse(res$formula$formula), collapse = "")
  expect_true(grepl("\\+ bhaz\\(", main_form_str))
  expect_true(grepl("gr = region", main_form_str))
})

test_that("Survival model with sparse event times falls back to Cox", {
  sparse_data <- test_data[1:5, ]
  sparse_data$time <- c(10, 10, 20, 20, 30)
  sparse_data$trt <- factor(c("Control", "Treatment", "Control", "Treatment", "Control"),
                             levels = c("Control", "Treatment"))

  expect_warning(
    res <- prepare_formula_model(
      data = sparse_data,
      response_formula = "Surv(time, status) ~ trt",
      response_type = "survival"
    ),
    regexp = "Not enough unique event times to compute quantile knots"
  )

  expect_valid_output(res, "survival")
  main_form_str <- deparse(res$formula$formula)
  expect_false(grepl("\\+ bhaz\\(", main_form_str))
})

test_that("Complex survival model with all features", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    unshrunk_terms_formula = "~ age + trt:region",
    shrunk_prognostic_formula = "~ 0 + subgroup1",
    shrunk_predictive_formula = "~ 0 + trt:subgroup2",
    stratification_formula = "~ region"
  )

  expect_valid_output(res, "survival")
  expect_true("unshrunktermeffect" %in% get_formula_rhs(res$formula))
  expect_true("shprogeffect" %in% get_formula_rhs(res$formula))
  expect_true("shpredeffect" %in% get_formula_rhs(res$formula))

  main_form_str <- paste(deparse(res$formula$formula), collapse = "")
  expect_true(grepl("gr = region", main_form_str))
})



# PART 3: STRATIFICATION


test_that("Continuous stratification (sigma)", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    stratification_formula = "~ region"
  )

  expect_valid_output(res, "continuous")
  expect_true("sigma" %in% names(res$formula$pforms))
  expect_true(grepl("sigma ~ region", deparse(res$formula$pforms$sigma)))
})

test_that("Count stratification (shape)", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "n_events ~ trt",
    response_type = "count",
    stratification_formula = "~ region"
  )

  expect_valid_output(res, "count")
  expect_true("shape" %in% names(res$formula$pforms))
  expect_true(grepl("shape ~ region", deparse(res$formula$pforms$shape)))
})



# PART 4: INTERACTION HANDLING


test_that("Colon notation interactions work", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:subgroup1",
    unshrunk_terms_formula = "~ trt:region"
  )

  expect_factor_with_contrasts(res$data, c("subgroup1", "region"))
  expect_equal(get_formula_rhs(res$formula), c("shpredeffect", "unshrunktermeffect"))
  expect_equal(get_formula_rhs(res$formula, "shpredeffect"), "trt:subgroup1")
  expect_equal(get_formula_rhs(res$formula, "unshrunktermeffect"), c("trt", "trt:region"))
})

test_that("Star notation interactions work and are preserved", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age + trt*region"
  )

  expect_factor_with_contrasts(res$data, "region")
  unshrunk_terms <- get_formula_rhs(res$formula, "unshrunktermeffect")
  expect_true("age" %in% unshrunk_terms)
  expect_true("trt" %in% unshrunk_terms)

})

test_that("Pipe-pipe notation creates random effects", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ (0 + trt || subgroup1)"
  )

  expect_true(res$has_random_effects)
  expect_true(is.factor(res$data$subgroup1))

  shrunk_pred_str <- paste(deparse(res$formula$pforms$shpredeffect), collapse = "")
  expect_true(grepl("\\|\\|", shrunk_pred_str))
})

test_that("Multiple random effects work", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ (0 + trt || subgroup1) + (0 + trt || subgroup2)"
  )

  expect_true(res$has_random_effects)
  expect_factor_with_contrasts(res$data, c("subgroup1", "subgroup2"))
})

test_that("Multiple interaction terms in same formula", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region + trt:subgroup1 + trt:subgroup2"
  )

  shrunk_pred_terms <- get_formula_rhs(res$formula, "shpredeffect")
  expect_true(all(c("trt:region", "trt:subgroup1", "trt:subgroup2") %in% shrunk_pred_terms))
  expect_factor_with_contrasts(res$data, c("region", "subgroup1", "subgroup2"))
})

test_that("Mixed syntaxes (colon + star + pipe) work together", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age + trt:region",
    shrunk_predictive_formula = "~ 0 + trt*subgroup1 + (0 + trt || subgroup2)"
  )

  expect_factor_with_contrasts(res$data, c("region", "subgroup1", "subgroup2"))
  expect_true(res$has_random_effects)
  expect_true("age" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
  expect_true("trt:region" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
  expect_true(length(get_formula_rhs(res$formula, "shpredeffect")) > 0)
})

test_that("Duplicate interaction terms are deduplicated", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age + trt:region"
  )

  unshrunk_terms <- get_formula_rhs(res$formula, "unshrunktermeffect")
  interaction_terms <- unshrunk_terms[grepl(":", unshrunk_terms)]
  expect_equal(length(interaction_terms), 1)
})



# PART 5: FORMULA COMPONENTS


test_that("All three formula types can be used together", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age + trt:region",
    shrunk_prognostic_formula = "~ 0 + subgroup1",
    shrunk_predictive_formula = "~ 0 + trt:subgroup2"
  )

  expect_true(all(c("unshrunktermeffect", "shprogeffect", "shpredeffect") %in%
                    get_formula_rhs(res$formula)))
  expect_true("age" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
  expect_true("trt:region" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
  expect_equal(get_formula_rhs(res$formula, "shprogeffect"), "subgroup1")
  expect_equal(get_formula_rhs(res$formula, "shpredeffect"), "trt:subgroup2")
})

test_that("Prognostic term overlaps throw error", {
  expect_error(
    prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      unshrunk_terms_formula = "~ age + region",
      shrunk_prognostic_formula = "~ 0 + age + subgroup1"
    ),
    regexp = "Variables cannot appear in both unshrunk_terms_formula and shrunk_prognostic_formula"
  )
})

test_that("Predictive interaction overlaps create duplicates", {
  expect_warning(
    res <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      unshrunk_terms_formula = "~ region",
      shrunk_predictive_formula = "~ 0 + trt:region"
    ),
    regexp = NA
  )

  # Original variable in unshrunk (dummy encoding)
  expect_true("region" %in% get_formula_rhs(res$formula, "unshrunktermeffect"))
  expect_dummy_contrasts(res$data, "region")
  
  # Duplicate variable with _onehot suffix in shrunk predictive (one-hot encoding)
  expect_equal(get_formula_rhs(res$formula, "shpredeffect"), "trt:region_onehot")
  expect_onehot_contrasts(res$data, "region_onehot")
  
  # Both variables should exist in data
  expect_true("region" %in% names(res$data))
  expect_true("region_onehot" %in% names(res$data))
})



# PART 6: INTERCEPT HANDLING


test_that("Intercept tracking with default (has intercept)", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age + region"
  )
  expect_true(res$has_intercept)
})

test_that("Intercept tracking with ~ 0 + (no intercept)", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ 0 + age + region"
  )
  expect_false(res$has_intercept)
})

test_that("Intercept tracking with ~ -1 + (no intercept)", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ -1 + age + region"
  )
  expect_false(res$has_intercept)
})

test_that("Survival models handle intercept correctly for Cox", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "Surv(time, status) ~ trt",
    response_type = "survival",
    unshrunk_terms_formula = "~ age"
  )

  expect_valid_output(res, "survival")
})



# PART 7: CONTRAST CODING


test_that("Shrunk terms get one-hot encoding, unshrunk get dummy", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ region",
    shrunk_prognostic_formula = "~ 0 + subgroup1"
  )

  expect_dummy_contrasts(res$data, "region")
  expect_onehot_contrasts(res$data, "subgroup1")
})

test_that("User-specified contrasts are preserved", {
  test_data_custom <- test_data
  custom_contrast <- contr.sum(levels(test_data_custom$region))
  contrasts(test_data_custom$region) <- custom_contrast

  res <- prepare_formula_model(
    data = test_data_custom,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ region"
  )

  expect_true(is.factor(res$data$region))
  expect_equal(contrasts(res$data$region), custom_contrast)
})

test_that("Non-factor variables converted to factors for interactions", {
  test_data_nofactor <- test_data
  test_data_nofactor$region <- as.character(test_data$region)
  test_data_nofactor$subgroup1 <- as.character(test_data$subgroup1)

  res <- prepare_formula_model(
    data = test_data_nofactor,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    shrunk_predictive_formula = "~ 0 + trt:region",
    unshrunk_terms_formula = "~ trt:subgroup1"
  )

  expect_factor_with_contrasts(res$data, c("region", "subgroup1"))
})



# PART 8: TREATMENT VARIABLE HANDLING


test_that("Numeric binary treatment handled correctly", {
  test_data_num <- test_data
  test_data_num$trt <- as.numeric(test_data$trt) - 1

  res <- prepare_formula_model(
    data = test_data_num,
    response_formula = "outcome ~ trt",
    response_type = "continuous"
  )

  expect_binary_treatment(res$data)
})

test_that("Character treatment converted correctly", {
  test_data_char <- test_data
  test_data_char$trt <- as.character(test_data$trt)

  res <- prepare_formula_model(
    data = test_data_char,
    response_formula = "outcome ~ trt",
    response_type = "continuous"
  )

  expect_binary_treatment(res$data)
})

test_that("Non-binary treatment throws error", {
  test_data_multi <- test_data
  test_data_multi$trt <- factor(c("A", "B", "C")[sample(1:3, nrow(test_data), replace = TRUE)])

  expect_error(
    prepare_formula_model(
      data = test_data_multi,
      response_formula = "outcome ~ trt",
      response_type = "continuous"
    ),
    regexp = "must have exactly 2 levels"
  )
})



# PART 9: MARGINALITY PRINCIPLE


test_that("Missing main effects generate warnings", {
  expect_message(
    res <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      shrunk_predictive_formula = "~ 0 + trt:subgroup1"
    ),
    regexp = "Marginality principle not followed.*subgroup1"
  )
})

test_that("Star notation excludes from marginality check", {
  expect_message(
    res <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      unshrunk_terms_formula = "~ trt*region",
      shrunk_predictive_formula = "~ 0 + trt:subgroup1"
    ),
    regexp = "Marginality principle not followed.*subgroup1",
    all = TRUE
  )

  expect_valid_output(res, "continuous")
})

test_that("No marginality warning when main effect present", {
  output1 <- capture.output({
    res1 <- suppressMessages(
      prepare_formula_model(
        data = test_data,
        response_formula = "outcome ~ trt",
        response_type = "continuous",
        unshrunk_terms_formula = "~ subgroup1",
        shrunk_predictive_formula = "~ 0 + trt:subgroup1"
      )
    )
  })
  expect_false(any(grepl("Marginality", output1)))

  output2 <- capture.output({
    res2 <- suppressMessages(
      prepare_formula_model(
        data = test_data,
        response_formula = "outcome ~ trt",
        response_type = "continuous",
        shrunk_prognostic_formula = "~ 0 + subgroup1",
        shrunk_predictive_formula = "~ 0 + trt:subgroup1"
      )
    )
  })
  expect_false(any(grepl("Marginality", output2)))
})

test_that("Star notation in mixed model doesn't warn about star variables", {
  output <- capture.output({
    res <- suppressMessages(
      prepare_formula_model(
        data = test_data,
        response_formula = "outcome ~ trt",
        response_type = "continuous",
        shrunk_predictive_formula = "~ 0 + trt*subgroup1"
      )
    )
  })

  expect_false(any(grepl("Marginality", output)))
})



# PART 10: WARNINGS


test_that("Intercept warnings for shrunk prognostic without ~ 0 +", {
  expect_warning(
    res <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      shrunk_prognostic_formula = "~ subgroup1"
    ),
    regexp = "contains an intercept.*consider removing"
  )
})

test_that("Intercept warnings for shrunk predictive without ~ 0 +", {
  expect_warning(
    res <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      shrunk_predictive_formula = "~ trt:subgroup1"
    ),
    regexp = "contains an intercept.*consider removing"
  )
})

test_that("Random effects without ~ 0 + generates warning", {
  expect_warning(
    res <- prepare_formula_model(
      data = test_data,
      response_formula = "outcome ~ trt",
      response_type = "continuous",
      shrunk_predictive_formula = "~ (trt || subgroup1)"
    ),
    regexp = "contains an intercept"
  )

  expect_true(res$has_random_effects)
  shrunk_pred_str <- paste(deparse(res$formula$pforms$shpredeffect), collapse = "")
  expect_true(grepl("\\|\\|", shrunk_pred_str))
})



# PART 11: ADDITIONAL OUTPUTS


test_that("Stan variable names are extracted", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age"
  )

  expect_true("stan_variable_names" %in% names(res))
  expect_type(res$stan_variable_names, "list")
})

test_that("All return values present and correct type", {
  res <- prepare_formula_model(
    data = test_data,
    response_formula = "outcome ~ trt",
    response_type = "continuous",
    unshrunk_terms_formula = "~ age"
  )

  expect_valid_output(res, "continuous")
  expect_true(res$has_intercept)
  expect_false(res$has_random_effects)
})



# PART 12: ERROR VALIDATION


test_that("Invalid response type", {
  expect_error(
    prepare_formula_model(test_data, "outcome ~ trt", response_type = "wrong"),
    regexp = "Must be element of set"
  )
})

test_that("Formula without tilde", {
  expect_error(
    prepare_formula_model(test_data, "outcome", response_type = "continuous"),
    regexp = "Must comply to pattern '~'"
  )
})

test_that("Treatment variable not in data", {
  expect_error(
    prepare_formula_model(test_data, "outcome ~ treatment_var", response_type = "continuous"),
    regexp = "subset of"
  )
})

test_that("Stratification variable not in data (survival)", {
  expect_error(
    prepare_formula_model(test_data, "Surv(time, status) ~ trt",
                         response_type = "survival",
                         stratification_formula = "~ missing_var"),
    regexp = "Must be a subset of"
  )
})

test_that("Stratification variable not in data (continuous)", {
  expect_error(
    prepare_formula_model(test_data, "outcome ~ trt",
                         response_type = "continuous",
                         stratification_formula = "~ missing_var"),
    regexp = "Must be a subset of"
  )
})

test_that("Survival formula doesn't start with Surv", {
  expect_error(
    prepare_formula_model(test_data, "time ~ trt", response_type = "survival"),
    regexp = "Must comply to pattern"
  )
})
