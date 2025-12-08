# Fit a Bayesian Hierarchical Model using brms

This function serves as a wrapper to fit a Bayesian model using the
\`brms\` package. It uses a formula and data prepared by
\`prepare_formula_model\` and simplifies the process of assigning
complex priors to the different non-linear components of the model.

## Usage

``` r
fit_brms_model(
  prepared_model,
  predictive_effect_priors = list(),
  prognostic_effect_priors = list(),
  stanvars = NULL,
  ...
)
```

## Arguments

- prepared_model:

  A list object returned from \`prepare_formula_model()\`. It must
  contain the elements \`formula\`, \`data\`, and \`response_type\`.

- predictive_effect_priors:

  A named list with elements \`shrunk\` and/or \`unshrunk\` containing
  the priors for predictive effects.

- prognostic_effect_priors:

  A named list with elements \`shrunk\`, \`unshrunk\`, and/or
  \`intercept\` containing the priors for prognostic effects.

- stanvars:

  An object created by \`brms::stanvar()\` to add custom Stan code.

- ...:

  Additional arguments passed directly to \`brms::brm()\`.

## Value

A fitted \`brmsfit\` object.

## Prior Specification

Priors for prognostic and predictive effects should be provided as named
lists. This makes the code explicit and prevents errors. For example:
\`prognostic_effect_priors = list(shrunk = "horseshoe(1)", unshrunk =
"normal(0, 5)", intercept = "normal(0, 10)")\`. The function accepts
prior definitions as character strings (e.g., \`"normal(0, 1)"\`) or as
full \`brmsprior\` objects, which is useful for defining complex
hierarchical or parameter-specific priors.

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

  # 2. Prepare the formula and data
  prepared_model <- prepare_formula_model(
    data = sim_data,
    response_formula_str = "Surv(time, status) ~ trt",
    shrunk_predictive_formula_str = "~ trt:subgroup",
    unshrunk_prognostic_formula_str = "~ age",
    shrunk_prognostic_formula_str = "~ region",
    response_type = "survival",
    stratification_formula_str = "~ region"
  )

  # 3. Fit the model
  if (FALSE) { # \dontrun{
  fit <- fit_brms_model(
    prepared_model = prepared_model,
    # Note: intercept prior is not needed for survival models
    prognostic_effect_priors = list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 2)"),
    predictive_effect_priors = list(shrunk = "horseshoe(1)"),
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
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'region'.
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: subgroup
```
