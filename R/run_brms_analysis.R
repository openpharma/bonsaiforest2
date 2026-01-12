#' Run a Full Bayesian Hierarchical Model Analysis
#'
#' This high-level wrapper function streamlines the entire process of preparing
#' and fitting a complex Bayesian hierarchical model with `brms`.
#'
#' This function is the main user-facing entry point. It first calls
#' `prepare_formula_model` to build the `brmsformula` and process the data,
#' then passes the results to `fit_brms_model` to run the analysis.
#'
#' @param data A data.frame containing all the necessary variables.
#' @param response_formula_str A character string for the response part, e.g.,
#'   "outcome ~ trt" or "Surv(time, status) ~ trt".
#' @param response_type The type of outcome variable. One of "binary", "count",
#'   "continuous", or "survival".
#' @param shrunk_predictive_formula_str Predictive terms to be shrunk ('shpredeffect').
#'   E.g., "~ trt:subgroup1" or "~ (trt || subgroup1)".
#' @param unshrunk_prognostic_formula_str Prognostic terms not to be shrunk
#'   ('unprogeffect'). E.g., "~ age + sex".
#' @param unshrunk_predictive_formula_str Predictive terms not to be shrunk
#'   ('unpredeffect'). E.g., "~ trt:important_subgroup" or "~ (trt || important_subgroup)".
#'   Note: Only one interaction syntax (: or ||) can be used across all predictive formulas.
#' @param shrunk_prognostic_formula_str Prognostic terms to be shrunk
#'   ('shprogeffect'). E.g., "~ region + center".
#' @param stratification_formula_str A formula string specifying a stratification
#'   variable, e.g., "~ strata_var".
#' @param predictive_effect_priors A named list with elements `shrunk` and/or `unshrunk`
#'   containing the priors for predictive effects. Can be strings or `brmsprior` objects.
#'   E.g., `list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 5)")`.
#' @param prognostic_effect_priors A named list with elements `shrunk`, `unshrunk` and/or `intercept`
#'   containing the priors for prognostic effects.
#'   E.g., `list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 10)")`.
#' @param stanvars An object created by `brms::stanvar()` to add custom Stan code,
#'   necessary for some hierarchical priors.
#' @param ... Additional arguments passed directly to `brms::brm()` (e.g.,
#'   `chains`, `iter`, `cores`, `backend`).
#'
#' @return A fitted `brmsfit` object.
#' @export
#'
#' @examples
#' if (require("brms") && require("survival")) {
#'   # 1. Create Sample Data
#'   set.seed(123)
#'   n <- 100
#'   sim_data <- data.frame(
#'     time = round(runif(n, 1, 100)),
#'     status = sample(0:1, n, replace = TRUE),
#'     trt = sample(0:1, n, replace = TRUE),
#'     age = rnorm(n, 50, 10),
#'     region = sample(c("A", "B"), n, replace = TRUE),
#'     subgroup = sample(c("S1", "S2", "S3"), n, replace = TRUE)
#'   )
#'   sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
#'   sim_data$region <- as.factor(sim_data$region)
#'   sim_data$subgroup <- as.factor(sim_data$subgroup)
#'
#'   # 2. Run the full analysis
#'   # We use \dontrun{} because fitting a model requires Stan compilation
#'   # which may fail in automated CI/CD environments.
#'   \dontrun{
#'   full_fit <- run_brms_analysis(
#'     data = sim_data,
#'     response_formula_str = "Surv(time, status) ~ trt",
#'     response_type = "survival",
#'     shrunk_predictive_formula_str = "~ trt:subgroup",
#'     unshrunk_prognostic_formula_str = "~ age",
#'     shrunk_prognostic_formula_str = "~ region",
#'     stratification_formula_str = "~ region",
#'     prognostic_effect_priors = list(
#'       shrunk = "normal(0, 1)",
#'       unshrunk = "normal(0, 5)"
#'     ),
#'     predictive_effect_priors = list(
#'       shrunk = "horseshoe(1)"
#'     ),
#'     chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
#'   )
#'
#'   print(full_fit)
#'   }
#' }
run_brms_analysis <- function(data,
                              response_formula_str,
                              response_type,
                              shrunk_predictive_formula_str = NULL,
                              unshrunk_prognostic_formula_str = NULL,
                              unshrunk_predictive_formula_str = NULL,
                              shrunk_prognostic_formula_str = NULL,
                              stratification_formula_str = NULL,
                              predictive_effect_priors = list(), # Use list() as default
                              prognostic_effect_priors = list(), # Use list() as default
                              stanvars = NULL,
                              ...) {

  # --- 1. Prepare the Formula and Data ---
  message("Step 1: Preparing formula and data...")
  prepared_model <- prepare_formula_model(
    data = data,
    response_formula_str = response_formula_str,
    shrunk_predictive_formula_str = shrunk_predictive_formula_str,
    unshrunk_prognostic_formula_str = unshrunk_prognostic_formula_str,
    unshrunk_predictive_formula_str = unshrunk_predictive_formula_str,
    shrunk_prognostic_formula_str = shrunk_prognostic_formula_str,
    response_type = response_type,
    stratification_formula_str = stratification_formula_str
  )

  # --- 2. Fit the Bayesian Model ---
  message("\nStep 2: Fitting the brms model...")

  # --- THIS IS THE MODIFIED PART ---
  model_fit <- fit_brms_model(
    prepared_model = prepared_model, # Pass the entire list
    predictive_effect_priors = predictive_effect_priors,
    prognostic_effect_priors = prognostic_effect_priors,
    stanvars = stanvars,
    ...
  )
  # --- END OF MODIFICATION ---

  # --- 3. Return the Final Model ---
  message("\nAnalysis complete.")
  return(model_fit)
}

