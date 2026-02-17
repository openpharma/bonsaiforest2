#' Run a Full Bayesian Hierarchical Model Analysis
#'
#' This high-level wrapper function streamlines the entire process of preparing
#' and fitting a complex Bayesian hierarchical model with `brms`. It serves as a
#' convenient single entry point for the complete workflow.
#'
#' This function is the main user-facing entry point. It first calls
#' `prepare_formula_model()` to build the `brmsformula` and process the data,
#' then passes the results to `fit_brms_model()` to run the analysis.
#'
#' @param data A data frame. Dataset containing all necessary variables including response,
#'   treatment, and covariates. Passed to `prepare_formula_model()` for preprocessing.
#' @param response_formula A formula object. Response specification defining the outcome and
#'   treatment. Examples: `outcome ~ trt` for continuous, or `Surv(time, status) ~ trt` for
#'   survival. Passed to `prepare_formula_model()`.
#' @param response_type A character string. Outcome type, one of `"binary"`, `"count"`,
#'   `"continuous"`, or `"survival"`. Determines the likelihood function and link used.
#'   Passed to `prepare_formula_model()`.
#' @param unshrunk_terms_formula A formula object or `NULL`. Unshrunk terms specification for
#'   the `unshrunktermeffect` component. May include main effects and treatment interactions
#'   without regularization (e.g., `~ age + sex + trt*region`). Passed to `prepare_formula_model()`.
#' @param shrunk_prognostic_formula A formula object or `NULL`. Prognostic effects for strong
#'   regularization in the `shprogeffect` component.
#'   - **Fixed effects (colon syntax):** Use `~ 0 + ...` syntax for one-hot encoding
#'     (e.g., `~ 0 + biomarker1 + biomarker2`).
#'   - **Random effects (one-way model, pipe-pipe syntax):** Use random effects notation
#'     (e.g., `~ (0 + 1 || biomarker)`) with automatic `normal(0, 1)` priors on SD scale.
#'   Passed to `prepare_formula_model()`.
#' @param shrunk_predictive_formula A formula object or `NULL`. Treatment interactions for strong
#'   regularization in the `shpredeffect` component.
#'   - **Fixed effects (colon syntax):** Use `~ 0 + ...` syntax for one-hot encoding
#'     (e.g., `~ 0 + trt:subgroup`).
#'   - **Random effects (one-way model, pipe-pipe syntax):** Use random effects notation
#'     (e.g., `~ (0 + trt || subgroup)`) with automatic `normal(0, 1)` priors on SD scale.
#'   Passed to `prepare_formula_model()`.
#' @param stratification_formula A formula object or `NULL`. Stratification variable specification
#'   (e.g., `~ strata_var`). For survival models, estimates separate baseline hazards per stratum.
#'   For other models, models varying distributional parameters. Passed to `prepare_formula_model()`.
#' @param intercept_prior A character string, `brmsprior` object, or `NULL`. Intercept prior for
#'   `unshrunktermeffect` component. Not used for survival models (Cox models have no intercept).
#'   Example: `"normal(0, 10)"`. Passed to `fit_brms_model()`.
#' @param unshrunk_prior A character string, `brmsprior` object, or `NULL`. Prior for unshrunk
#'   terms (non-intercept coefficients in `unshrunktermeffect` component).
#'   Example: `"normal(0, 2.5)"`. Passed to `fit_brms_model()`.
#' @param shrunk_prognostic_prior A character string, `brmsprior` object, or `NULL`. Prior for
#'   regularized prognostic effects in `shprogeffect` component.
#'   - **Fixed effects:** Typically strong regularization like `"horseshoe(scale_global = 1)"`.
#'   - **Random effects (one-way model):** Automatically uses `normal(0, 1)` priors on SD scale.
#'   Passed to `fit_brms_model()`.
#' @param shrunk_predictive_prior A character string, `brmsprior` object, or `NULL`. Prior for
#'   regularized predictive effects (treatment interactions) in `shpredeffect` component.
#'   - **Fixed effects:** Typically strong regularization like
#'     `"horseshoe(scale_global = 0.5)"`.
#'   - **Random effects (one-way model):** Automatically uses `normal(0, 1)` priors on SD scale.
#'   Passed to `fit_brms_model()`.
#' @param stanvars A `stanvars` object or `NULL`. Custom Stan code created via `brms::stanvar()`
#'   for implementing hierarchical priors or other advanced Stan functionality. Passed to `fit_brms_model()`.
#' @param ... Additional arguments passed to `brms::brm()` via `fit_brms_model()`. Common
#'   arguments include `chains`, `iter`, `cores`, `backend`, and `refresh`.
#'
#' @return `brmsfit`. Fitted Bayesian model object with attached attributes: `response_type`,
#'   `model_data` (the processed data used for fitting), and `trt_var` (treatment variable name)
#'   for use by downstream analysis functions.
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
#'   # 2. Run the full analysis using fixed effects (colon syntax)
#'   \dontrun{
#'   full_fit_fixed <- run_brms_analysis(
#'     data = sim_data,
#'     response_formula = Surv(time, status) ~ trt,
#'     response_type = "survival",
#'     unshrunk_terms_formula = ~ age,
#'     shrunk_prognostic_formula = ~ 0 + region,
#'     shrunk_predictive_formula = ~ 0 + trt:subgroup,
#'     stratification_formula = ~ region,
#'     unshrunk_prior = "normal(0, 2)",
#'     shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
#'     shrunk_predictive_prior = "horseshoe(scale_global = 1)",
#'     backend = "cmdstanr",
#'     chains = 2, iter = 1000, warmup = 500, refresh = 0
#'   )
#'
#'   print(full_fit_fixed)
#'
#'   # 3. Alternative: One-way model with random effects (pipe-pipe syntax)
#'   # Random effects automatically get normal(0, 1) priors on SD scale
#'   full_fit_oneway <- run_brms_analysis(
#'     data = sim_data,
#'     response_formula = Surv(time, status) ~ trt,
#'     response_type = "survival",
#'     unshrunk_terms_formula = ~ age,
#'     shrunk_prognostic_formula = ~ (0 + 1 || region),
#'     shrunk_predictive_formula = ~ (0 + trt || subgroup),
#'     stratification_formula = ~ region,
#'     backend = "cmdstanr",
#'     chains = 2, iter = 1000, warmup = 500, refresh = 0
#'   )
#'
#'   print(full_fit_oneway)
#'   }
#' }
run_brms_analysis <- function(data,
                              response_formula,
                              response_type = c("binary", "count", "continuous", "survival"),
                              unshrunk_terms_formula = NULL,
                              shrunk_prognostic_formula = NULL,
                              shrunk_predictive_formula = NULL,
                              stratification_formula = NULL,
                              intercept_prior = NULL,
                              unshrunk_prior = NULL,
                              shrunk_prognostic_prior = NULL,
                              shrunk_predictive_prior = NULL,
                              stanvars = NULL,
                              ...) {

  # --- 1. Prepare the Formula and Data ---
  # Builds the brmsformula with up to three components:
  # - unshrunktermeffect: all unshrunk terms (main effects and interactions without regularization)
  # - shprogeffect: shrunk prognostic effects (main effects with strong regularization)
  # - shpredeffect: shrunk predictive effects (treatment interactions with strong regularization)
  # Supports both fixed effects (colon syntax) and one-way models (random effects via pipe-pipe syntax).
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
  # For fixed effects: uses user-specified or default priors (e.g., horseshoe)
  # For random effects: automatically uses normal(0, 1) priors on SD scale
  message("\nStep 2: Fitting the brms model...")

  model_fit <- fit_brms_model(
    prepared_model = prepared_model,
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

