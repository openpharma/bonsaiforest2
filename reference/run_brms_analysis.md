# Run a Full Bayesian Hierarchical Model Analysis

This high-level wrapper function streamlines the entire process of
preparing and fitting a complex Bayesian hierarchical model with `brms`.

## Usage

``` r
run_brms_analysis(
  data,
  response_formula,
  response_type = c("binary", "count", "continuous", "survival"),
  unshrunk_terms_formula = NULL,
  shrunk_prognostic_formula = NULL,
  shrunk_predictive_formula = NULL,
  stratification_formula = NULL,
  sigma_ref = 1,
  intercept_prior = NULL,
  unshrunk_prior = NULL,
  shrunk_prognostic_prior = NULL,
  shrunk_predictive_prior = NULL,
  stanvars = NULL,
  ...
)
```

## Arguments

- data:

  A data frame. Dataset containing all necessary variables including
  response, treatment, and covariates. Passed to
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md)
  for preprocessing.

- response_formula:

  A formula object. Response specification defining the outcome and
  treatment. Examples: `outcome ~ trt` for continuous, or
  `Surv(time, status) ~ trt` for survival. Passed to
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md).

- response_type:

  A character string. Outcome type, one of `"binary"`, `"count"`,
  `"continuous"`, or `"survival"`. Determines the likelihood function
  and link used. Passed to
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md).

- unshrunk_terms_formula:

  A formula object or `NULL`. Unshrunk terms specification for the
  `unshrunktermeffect` component. May include main effects and treatment
  interactions without regularization (e.g.,
  `~ age + sex + trt*region`). Passed to
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md).

- shrunk_prognostic_formula:

  A formula object or `NULL`. Prognostic effects for strong
  regularization in the `shprogeffect` component. For colon syntax
  (GLOBAL), use `~ 0 + ...` syntax for one-hot encoding (e.g.,
  `~ 0 + biomarker1 + biomarker2`). For pipe-pipe syntax (OVAT), use
  random effects notation (e.g., `~ (1 || biomarker)`). Passed to
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md).

- shrunk_predictive_formula:

  A formula object or `NULL`. Treatment interactions for strong
  regularization in the `shpredeffect` component. For colon syntax
  (GLOBAL), use `~ 0 + ...` syntax for one-hot encoding (e.g.,
  `~ 0 + trt:subgroup`). For pipe-pipe syntax (OVAT), use random effects
  notation (e.g., `~ (trt || subgroup)`). Passed to
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md).

- stratification_formula:

  A formula object or `NULL`. Stratification variable specification
  (e.g., `~ strata_var`). For survival models, estimates separate
  baseline hazards per stratum. For other models, models varying
  distributional parameters. Passed to
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md).

- sigma_ref:

  A numeric scalar. Reference scale for prior specification. Default
  is 1. For continuous outcomes, recommended values are: (1) Assumed
  standard deviation from trial protocol (preferred), or (2)
  `sd(outcome_variable)` if protocol value unavailable. For binary,
  survival and count outcomes, the default value of 1 is typically
  appropriate. Can be referenced in prior strings using the placeholder
  `sigma_ref` (e.g., `"normal(0, 2.5 * sigma_ref)"`). Passed to
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md).

- intercept_prior:

  A character string, `brmsprior` object, or `NULL`. Intercept prior for
  `unshrunktermeffect` component. Supports `sigma_ref` placeholder
  substitution. Not used for survival models (Cox models have no
  intercept). Example: `"normal(0, 10 * sigma_ref)"`. Passed to
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md).

- unshrunk_prior:

  A character string, `brmsprior` object, or `NULL`. Prior for unshrunk
  terms (non-intercept coefficients in `unshrunktermeffect` component).
  Supports `sigma_ref` placeholder substitution. Example:
  `"normal(0, 2.5 * sigma_ref)"`. Passed to
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md).

- shrunk_prognostic_prior:

  A character string, `brmsprior` object, or `NULL`. Prior for
  regularized prognostic effects in `shprogeffect` component. For fixed
  effects (colon syntax), typically strong regularization like
  `"horseshoe(scale_global = sigma_ref)"`. For random effects (pipe-pipe
  syntax), automatically uses `normal(0, sigma_ref)` priors on the
  standard deviation scale. Passed to
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md).

- shrunk_predictive_prior:

  A character string, `brmsprior` object, or `NULL`. Prior for
  regularized predictive effects (treatment interactions) in
  `shpredeffect` component. For fixed effects (colon syntax), typically
  strong regularization like
  `"horseshoe(scale_global = 0.5 * sigma_ref)"`. For random effects
  (pipe-pipe syntax), automatically uses `normal(0, sigma_ref)` priors
  on the standard deviation scale. Passed to
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md).

- stanvars:

  A `stanvars` object or `NULL`. Custom Stan code created via
  [`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
  for implementing hierarchical priors or other advanced Stan
  functionality. Passed to
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md).

- ...:

  Additional arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) via
  [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md).
  Common arguments include `chains`, `iter`, `cores`, `backend`, and
  `refresh`.

## Value

`brmsfit`. Fitted Bayesian model object with stored attributes for
downstream analysis.

## Details

This function is the main user-facing entry point. It first calls
[`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md)
to build the `brmsformula` and process the data, then passes the results
to
[`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)
to run the analysis.

## Examples

``` r
if (require("brms") && require("survival")) {
  # 1. Create Sample Data
  set.seed(123)
  n <- 100
  sim_data <- data.frame(
    time = round(runif(n, 1, 100)),
    status = sample(0:1, n, replace = TRUE),
    trt = sample(0:1, n, replace = TRUE),
    age = rnorm(n, 50, 10),
    region = sample(c("A", "B"), n, replace = TRUE),
    subgroup = sample(c("S1", "S2", "S3"), n, replace = TRUE)
  )
  sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
  sim_data$region <- as.factor(sim_data$region)
  sim_data$subgroup <- as.factor(sim_data$subgroup)

  # 2. Run the full analysis using GLOBAL fixed effects (colon syntax)
  if (FALSE) { # \dontrun{
  # For survival models, typically use sigma_ref = 1 (default)
  # For continuous outcomes, use protocol sigma (preferred) or sd(outcome) (fallback)
  # Example: sigma_ref <- 12.5  # from protocol
  # OR: sigma_ref <- sd(sim_data$outcome)  # if protocol value unavailable

  full_fit_global <- run_brms_analysis(
    data = sim_data,
    response_formula = Surv(time, status) ~ trt,
    response_type = "survival",
    unshrunk_terms_formula = ~ age,
    shrunk_prognostic_formula = ~ 0 + region,
    shrunk_predictive_formula = ~ 0 + trt:subgroup,
    stratification_formula = ~ region,
    sigma_ref = 1,  # Default for survival
    unshrunk_prior = "normal(0, 2 * sigma_ref)",
    shrunk_prognostic_prior = "horseshoe(scale_global = sigma_ref)",
    shrunk_predictive_prior = "horseshoe(scale_global = sigma_ref)",
    chains = 1, iter = 50, warmup = 10, refresh = 0 # For quick example
  )

  print(full_fit_global)

  # 3. Alternative: OVAT random effects (pipe-pipe syntax)
  # Random effects automatically get normal(0, sigma_ref) priors on SD scale
  full_fit_ovat <- run_brms_analysis(
    data = sim_data,
    response_formula = Surv(time, status) ~ trt,
    response_type = "survival",
    unshrunk_terms_formula = ~ age,
    shrunk_prognostic_formula = ~ (1 || region),
    shrunk_predictive_formula = ~ (trt || subgroup),
    stratification_formula = ~ region,
    sigma_ref = 1,
    chains = 1, iter = 50, warmup = 10, refresh = 0
  )

  print(full_fit_ovat)
  } # }
}
```
