#' Run a Full Bayesian Hierarchical Model Analysis
#'
#' This high-level wrapper function streamlines the entire process of preparing
#' and fitting a complex Bayesian hierarchical model with `brms`.
#'
#' This function is the main user-facing entry point. It first calls
#' `prepare_formula_model` to build the `brmsformula` and process the data,
#' then passes the results to `fit_brms_model` to run the analysis.
#'
#' @param data `data.frame`. Dataset containing all necessary variables.
#' @param response_formula `formula`. Response specification (e.g., `outcome ~ trt` or `Surv(time, status) ~ trt`).
#' @param response_type `character(1)`. Outcome type: `"binary"`, `"count"`, `"continuous"`, or `"survival"`.
#' @param unshrunk_terms_formula `formula` or `NULL`. Unshrunk terms specification (`unshrunktermeffect`).
#'   May include main effects and treatment interactions without regularization (e.g., `~ age + sex + trt*region`).
#' @param shrunk_prognostic_formula `formula` or `NULL`. Prognostic effects for strong regularization
#'   (`shprogeffect`). Must use `~ 0 + ...` syntax (e.g., `~ 0 + biomarker1 + biomarker2`).
#' @param shrunk_predictive_formula `formula` or `NULL`. Treatment interactions for strong regularization
#'   (`shpredeffect`). Must use `~ 0 + ...` syntax (e.g., `~ 0 + trt:subgroup` or `~ (0 + trt || subgroup)`).
#' @param stratification_formula `formula` or `NULL`. Stratification variable specification (e.g., `~ strata_var`).
#' @param sigma_ref `numeric(1)`. Reference scale for priors (REQUIRED). For continuous/count outcomes,
#'   typically `sd(outcome_variable)`. For binary/survival, typically 1. Referenced in prior
#'   expressions (e.g., `"normal(0, 2.5 * sigma_ref)"`).
#' @param intercept_prior `character(1)` or `brmsprior` or `NULL`. Intercept prior for `unshrunktermeffect`.
#'   Not used for survival models. Example: `"normal(0, 10 * sigma_ref)"`.
#' @param unshrunk_prior `character(1)` or `brmsprior` or `NULL`. Prior for unshrunk terms
#'   (non-intercept coefficients in `unshrunktermeffect`). Example: `"normal(0, 2.5 * sigma_ref)"`.
#' @param shrunk_prognostic_prior `character(1)` or `brmsprior` or `NULL`. Prior for regularized
#'   prognostic effects. Typically strong regularization. Example: `"horseshoe(scale_global = sigma_ref)"`.
#' @param shrunk_predictive_prior `character(1)` or `brmsprior` or `NULL`. Prior for regularized
#'   predictive effects (treatment interactions). Example: `"horseshoe(scale_global = 0.5 * sigma_ref)"`.
#' @param stanvars `stanvars` or `NULL`. Custom Stan code via `brms::stanvar()` for hierarchical priors.
#' @param ... Additional arguments for `brms::brm()` (e.g., `chains`, `iter`, `cores`, `backend`).
#'
#' @return `brmsfit`. Fitted Bayesian model object with stored attributes for downstream analysis.
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
#'     response_formula = Surv(time, status) ~ trt,
#'     response_type = "survival",
#'     unshrunk_terms_formula = ~ age,
#'     shrunk_prognostic_formula = ~ region,
#'     shrunk_predictive_formula = ~ trt:subgroup,
#'     stratification_formula = ~ region,
#'     sigma_ref = 1,
#'     unshrunk_prior = "normal(0, 2 * sigma_ref)",
#'     shrunk_prognostic_prior = "normal(0, sigma_ref)",
#'     shrunk_predictive_prior = "horseshoe(scale_global = sigma_ref)",
#'     chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
#'   )
#'
#'   print(full_fit)
#'   }
#' }
run_brms_analysis <- function(data,
                              response_formula,
                              response_type = c("binary", "count", "continuous", "survival"),
                              unshrunk_terms_formula = NULL,
                              shrunk_prognostic_formula = NULL,
                              shrunk_predictive_formula = NULL,
                              stratification_formula = NULL,
                              sigma_ref = NULL,
                              intercept_prior = NULL,
                              unshrunk_prior = NULL,
                              shrunk_prognostic_prior = NULL,
                              shrunk_predictive_prior = NULL,
                              stanvars = NULL,
                              ...) {

  # --- 1. Prepare the Formula and Data ---
  # Builds the brmsformula with three components:
  # - unshrunktermeffect: all unshrunk terms (weakly informative priors)
  # - shprogeffect: shrunk prognostic effects (strong regularization)
  # - shpredeffect: shrunk predictive effects (strong regularization)
  message("Step 1: Preparing formula and data...")
  prepared_model <- prepare_formula_model(
    data = data,
    response_formula = response_formula,
    unshrunk_terms_formula = unshrunk_terms_formula,
    shrunk_prognostic_formula = shrunk_prognostic_formula,
    shrunk_predictive_formula = shrunk_predictive_formula,
    response_type = response_type,
    stratification_formula = stratification_formula
  )

  # --- 2. Fit the Bayesian Model ---
  # Assigns appropriate priors to each formula component and calls brms::brm()
  message("\nStep 2: Fitting the brms model...")

  model_fit <- fit_brms_model(
    prepared_model = prepared_model,
    sigma_ref = sigma_ref,
    intercept_prior = intercept_prior,
    unshrunk_prior = unshrunk_prior,
    shrunk_prognostic_prior = shrunk_prognostic_prior,
    shrunk_predictive_prior = shrunk_predictive_prior,
    stanvars = stanvars,
    ...
  )

  # --- 3. Return the Final Model ---
  message("\nAnalysis complete.")
  return(model_fit)
}

