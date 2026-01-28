# Fit a Bayesian Hierarchical Model using brms

This function serves as a wrapper to fit a Bayesian model using the
\`brms\` package. It uses a formula and data prepared by
\`prepare_formula_model\` and simplifies the process of assigning
complex priors to the different non-linear components of the model.

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

  \`list\`. Object returned from \`prepare_formula_model()\` containing
  elements \`formula\` (\`brmsformula\`), \`data\` (\`data.frame\`),
  \`response_type\` (\`character\`), and \`trt_var\` (\`character\`).
  The \`trt_var\` element is stored as an attribute on the fitted model
  for downstream functions.

- sigma_ref:

  \`numeric(1)\`. Reference scale for prior specification (REQUIRED).
  Recommended: (1) Standard deviation from trial protocol (preferred),
  or (2) \`sd(outcome_variable)\` if protocol value unavailable
  (fallback). For binary/survival outcomes, typically use 1. Can be
  referenced in prior strings (e.g., \`"normal(0, 2.5 \* sigma_ref)"\`).

- intercept_prior:

  \`character(1)\` or \`brmsprior\` or \`NULL\`. Prior specification for
  model intercept (e.g., \`"normal(0, 10)"\`).

- unshrunk_prior:

  \`character(1)\` or \`brmsprior\` or \`NULL\`. Prior specification for
  unshrunk terms (\`unshrunktermeffect\` component).

- shrunk_prognostic_prior:

  \`character(1)\` or \`brmsprior\` or \`NULL\`. Prior specification for
  regularized prognostic effects.

- shrunk_predictive_prior:

  \`character(1)\` or \`brmsprior\` or \`NULL\`. Prior specification for
  regularized predictive effects (treatment interactions).

- stanvars:

  \`stanvars\` or \`NULL\`. Custom Stan code created via
  \`brms::stanvar()\`.

- ...:

  Additional arguments passed to \`brms::brm()\` (e.g., \`chains\`,
  \`iter\`, \`cores\`).

## Value

\`brmsfit\`. Fitted Bayesian model object with attributes:
\`response_type\` (\`character\`), \`model_data\` (\`data.frame\`), and
\`trt_var\` (\`character\`) for downstream analysis functions.

## Prior Specification

Priors are specified separately for each model component using dedicated
parameters. This provides explicit control and prevents errors. The
function accepts prior definitions as character strings (e.g.,
\`"normal(0, 2.5 \* sigma_ref)"\`) or as \`brmsprior\` objects for
complex hierarchical or parameter-specific priors.

\*\*Available prior parameters:\*\* - \`intercept_prior\`: Prior for the
intercept in the unshrunk component - \`unshrunk_prior\`: Prior for all
unshrunk terms (main effects you trust) - \`shrunk_prognostic_prior\`:
Prior for shrunk prognostic effects (regularized) -
\`shrunk_predictive_prior\`: Prior for shrunk predictive effects
(regularized)

\*\*Using sigma_ref in priors:\*\* You can reference \`sigma_ref\` in
prior strings, and it will be automatically substituted. For example:
\`"normal(0, 2.5 \* sigma_ref)"\` becomes \`"normal(0, 8)"\` when
sigma_ref = 3.2.

\*\*Default Priors by Response Type:\*\*

If you don't specify priors, the function uses sensible defaults:

\*Continuous outcomes:\* - \`intercept_prior\`: \`normal(outcome_mean, 5
\* sigma_ref)\` - \`unshrunk_prior\`: \`normal(0, 5 \* sigma_ref)\` -
\`shrunk_prognostic_prior\`: \`horseshoe(1)\` -
\`shrunk_predictive_prior\`: \`horseshoe(1)\` - \`sigma\`: \`normal(0,
sigma_ref)\`

\*Binary outcomes:\* - \`intercept_prior\`: \`normal(0, 5 \*
sigma_ref)\` - \`unshrunk_prior\`: \`normal(0, 5 \* sigma_ref)\` -
\`shrunk_prognostic_prior\`: \`horseshoe(1)\` -
\`shrunk_predictive_prior\`: \`horseshoe(1)\`

\*Count outcomes:\* - \`intercept_prior\`: \`normal(0, 5 \*
sigma_ref)\` - \`unshrunk_prior\`: \`normal(0, 5 \* sigma_ref)\` -
\`shrunk_prognostic_prior\`: \`horseshoe(1)\` -
\`shrunk_predictive_prior\`: \`horseshoe(1)\` - \`shape\`: Uses brms
default prior for negative binomial shape parameter

\*Survival outcomes:\* - No intercept (Cox models don't have
intercepts) - \`unshrunk_prior\`: \`normal(0, 5 \* sigma_ref)\` -
\`shrunk_prognostic_prior\`: \`horseshoe(1)\` -
\`shrunk_predictive_prior\`: \`horseshoe(1)\`

\*\*Example:\*\* “\` fit_brms_model( prepared_model = prepared,
sigma_ref = 12.5, \# From trial protocol (preferred) \# OR sigma_ref =
sd(data\$outcome), \# If protocol value unavailable unshrunk_prior =
"normal(0, 2.5 \* sigma_ref)", shrunk_prognostic_prior =
"horseshoe(scale_global = 0.5 \* sigma_ref)", shrunk_predictive_prior =
"horseshoe(scale_global = 0.1 \* sigma_ref)" ) “\`

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
    unshrunk_prognostic_formula = ~ age,
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
#> Error in prepare_formula_model(data = sim_data, response_formula = Surv(time,     status) ~ trt, shrunk_predictive_formula = ~trt:subgroup,     unshrunk_prognostic_formula = ~age, shrunk_prognostic_formula = ~region,     response_type = "survival", stratification_formula = ~region): unused argument (unshrunk_prognostic_formula = ~age)
```
