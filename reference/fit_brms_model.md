# Fit a Bayesian Hierarchical Model using brms

This function serves as a wrapper to fit a Bayesian model using the
`brms` package. It uses a formula and data prepared by
`prepare_formula_model` and simplifies the process of assigning complex
priors to the different non-linear components of the model.

## Usage

``` r
fit_brms_model(
  prepared_model,
  sigma_ref = NULL,
  intercept_prior = NULL,
  unshrunk_prior = NULL,
  shrunk_prognostic_prior = NULL,
  shrunk_predictive_prior = NULL,
  stanvars = NULL,
  ...
)
```

## Arguments

- prepared_model:

  A list object. Object returned from
  [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md)
  containing the elements: `formula` (a `brmsformula` object), `data` (a
  modified `data.frame` with appropriate contrast coding),
  `response_type` (a character string), and `trt_var` (a character
  string). The `trt_var` element is stored as an attribute on the fitted
  model for use by downstream functions like
  [`estimate_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/estimate_subgroup_effects.md).

- sigma_ref:

  A numeric scalar (REQUIRED). Reference scale for prior specification.
  Recommended values: (1) Standard deviation from trial protocol
  (preferred), or (2) `sd(outcome_variable)` if protocol value
  unavailable (fallback). For binary and survival outcomes, typically
  use 1. Can be referenced in prior strings using the placeholder
  `sigma_ref` (e.g., `"normal(0, 2.5 * sigma_ref)"`).

- intercept_prior:

  A character string, `brmsprior` object, or `NULL`. Prior specification
  for the model intercept in the `unshrunktermeffect` component.
  Supports `sigma_ref` placeholder substitution. Not used for survival
  models (Cox models have no intercept). Example: `"normal(0, 10)"`.

- unshrunk_prior:

  A character string, `brmsprior` object, or `NULL`. Prior specification
  for unshrunk terms in the `unshrunktermeffect` component (excludes
  intercept). These are main effects and interactions for which no
  regularization is desired. Supports `sigma_ref` placeholder
  substitution. Example: `"normal(0, 2.5 * sigma_ref)"`.

- shrunk_prognostic_prior:

  A character string, `brmsprior` object, or `NULL`. Prior specification
  for regularized prognostic effects in the `shprogeffect` component.
  These are main effects for which strong shrinkage/regularization is
  desired. Supports `sigma_ref` placeholder substitution. Example:
  `"horseshoe(1)"` or `"normal(0, sigma_ref)"`.

- shrunk_predictive_prior:

  A character string, `brmsprior` object, or `NULL`. Prior specification
  for regularized predictive effects (treatment interactions) in the
  `shpredeffect` component. These are treatment interactions for which
  strong shrinkage/regularization is desired. Supports `sigma_ref`
  placeholder substitution. Example:
  `"horseshoe(scale_global = 0.5 * sigma_ref)"`.

- stanvars:

  A `stanvars` object or `NULL`. Custom Stan code created via
  [`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
  for implementing hierarchical priors or other advanced Stan
  functionality. See the brms documentation for details on creating
  `stanvars` objects.

- ...:

  Additional arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Common arguments include `chains` (number of MCMC chains), `iter`
  (number of iterations per chain), `cores` (number of CPU cores to
  use), `backend` (Stan backend), and `refresh` (progress update
  frequency).

## Value

`brmsfit`. Fitted Bayesian model object with attributes: `response_type`
(`character`), `model_data` (`data.frame`), and `trt_var` (`character`)
for downstream analysis functions.

## Prior Specification

Priors are specified separately for each model component using dedicated
parameters. This provides explicit control and prevents errors. The
function accepts prior definitions as character strings (e.g.,
`"normal(0, 2.5 * sigma_ref)"`) or as `brmsprior` objects for complex
hierarchical or parameter-specific priors.

**Available prior parameters:**

- `intercept_prior`: Prior for the intercept in the unshrunk component

- `unshrunk_prior`: Prior for all unshrunk terms (main effects you
  trust)

- `shrunk_prognostic_prior`: Prior for shrunk prognostic effects
  (regularized)

- `shrunk_predictive_prior`: Prior for shrunk predictive effects
  (regularized)

**Using sigma_ref in priors:** You can reference `sigma_ref` in prior
strings, and it will be automatically substituted. For example:
`"normal(0, 2.5 * sigma_ref)"` becomes `"normal(0, 8)"` when sigma_ref =
3.2.

**Default Priors by Response Type:**

If you don't specify priors, the function uses sensible defaults:

*Continuous outcomes:*

- `intercept_prior`: `normal(outcome_mean, 5 * sigma_ref)`

- `unshrunk_prior`: `normal(0, 5 * sigma_ref)`

- `shrunk_prognostic_prior`: `horseshoe(1)`

- `shrunk_predictive_prior`: `horseshoe(1)`

- `sigma`: `normal(0, sigma_ref)`

*Binary outcomes:*

- `intercept_prior`: `normal(0, 5 * sigma_ref)`

- `unshrunk_prior`: `normal(0, 5 * sigma_ref)`

- `shrunk_prognostic_prior`: `horseshoe(1)`

- `shrunk_predictive_prior`: `horseshoe(1)`

*Count outcomes:*

- `intercept_prior`: `normal(0, 5 * sigma_ref)`

- `unshrunk_prior`: `normal(0, 5 * sigma_ref)`

- `shrunk_prognostic_prior`: `horseshoe(1)`

- `shrunk_predictive_prior`: `horseshoe(1)`

- `shape`: Uses brms default prior for negative binomial shape parameter

*Survival outcomes:*

- No intercept (Cox models don't have intercepts)

- `unshrunk_prior`: `normal(0, 5 * sigma_ref)`

- `shrunk_prognostic_prior`: `horseshoe(1)`

- `shrunk_predictive_prior`: `horseshoe(1)`

**Example:**

    fit_brms_model(
      prepared_model = prepared,
      sigma_ref = 12.5,  # From trial protocol (preferred)
      # OR sigma_ref = sd(data$outcome),  # If protocol value unavailable
      unshrunk_prior = "normal(0, 2.5 * sigma_ref)",
      shrunk_prognostic_prior = "horseshoe(scale_global = 0.5 * sigma_ref)",
      shrunk_predictive_prior = "horseshoe(scale_global = 0.1 * sigma_ref)"
    )

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

  # 2. Prepare the formula and data using colon syntax
  prepared_model <- prepare_formula_model(
    data = sim_data,
    response_formula = Surv(time, status) ~ trt,
    shrunk_predictive_formula = ~ trt:subgroup,
    unshrunk_terms_formula = ~ age,
    shrunk_prognostic_formula = ~ region,
    response_type = "survival",
    stratification_formula = ~ region
  )

  # 3. Fit the model
  if (FALSE) { # \dontrun{
  # For survival models, typically use sigma_ref = 1
  # For continuous outcomes, use protocol sigma (preferred) or sd(outcome) (fallback)
  # Example: sigma_ref <- 12.5  # from protocol
  # OR: sigma_ref <- sd(sim_data$outcome)  # if protocol value unavailable

  fit <- fit_brms_model(
    prepared_model = prepared_model,
    sigma_ref = 1,
    unshrunk_prior = "normal(0, 2 * sigma_ref)",
    shrunk_prognostic_prior = "horseshoe(scale_global = sigma_ref)",
    shrunk_predictive_prior = "horseshoe(scale_global = sigma_ref)",
    chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
  )

  print(fit)
  } # }
}
#> Loading required package: brms
#> Loading required package: Rcpp
#> Loading 'brms' package (version 2.23.0). Useful instructions
#> can be found by typing help('brms'). A more detailed introduction
#> to the package is available through vignette('brms_overview').
#> 
#> Attaching package: ‘brms’
#> The following object is masked from ‘package:stats’:
#> 
#>     ar
#> Loading required package: survival
#> 
#> Attaching package: ‘survival’
#> The following object is masked from ‘package:brms’:
#> 
#>     kidney
#> Converting treatment variable 'trt' to numeric binary (0/1). '0' = 0, '1' = 1
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'region'.
#> Note: Treatment 'trt' automatically added to unshrunk terms.
#> Note: Applied one-hot encoding to shrunken factor 'subgroup' (will be used with ~ 0 + ...)
#> Note: Marginality principle not followed - interaction term 'subgroup' is used without its main effect. Consider adding 'subgroup' to prognostic terms for proper model hierarchy.
#> Note: Applied one-hot encoding to shrunken factor 'region' (will be used with ~ 0 + ...)
#> DEBUG: Creating sub-formulas...
#>   - all_unshrunk_terms: age, trt
#>   - shrunk_prog_terms: region
#>   - shrunk_pred_formula: trt:subgroup
#> Warning: Formula 'shprogeffect' contains an intercept. For proper regularization/interpretation, consider removing it by adding '~ 0 + ...' or '~ -1 + ...' to your input formula.
#> Warning: Formula 'shpredeffect' contains an intercept. For proper regularization/interpretation, consider removing it by adding '~ 0 + ...' or '~ -1 + ...' to your input formula.
#> DEBUG: Final formula object:
#> time | cens(1 - status) + bhaz(Boundary.knots = c(0.02, 99.98), knots = c(24, 46, 69), intercept = FALSE, gr = region) ~ unshrunktermeffect + shprogeffect + shpredeffect 
#> unshrunktermeffect ~ 0 + age + trt
#> shprogeffect ~ region
#> shpredeffect ~ trt:subgroup
```
