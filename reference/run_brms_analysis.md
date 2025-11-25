# Run a Full Bayesian Hierarchical Model Analysis

This high-level wrapper function streamlines the entire process of
preparing and fitting a complex Bayesian hierarchical model with
\`brms\`.

## Usage

``` r
run_brms_analysis(
  data,
  response_formula_str,
  response_type,
  shrunk_predictive_formula_str = NULL,
  unshrunk_prognostic_formula_str = NULL,
  unshrunk_predictive_formula_str = NULL,
  shrunk_prognostic_formula_str = NULL,
  stratification_formula_str = NULL,
  predictive_effect_priors = list(),
  prognostic_effect_priors = list(),
  stanvars = NULL,
  ...
)
```

## Arguments

- data:

  A data.frame containing all the necessary variables.

- response_formula_str:

  A character string for the response part, e.g., "outcome ~ trt" or
  "Surv(time, status) ~ trt".

- response_type:

  The type of outcome variable. One of "binary", "count", "continuous",
  or "survival".

- shrunk_predictive_formula_str:

  Predictive terms to be shrunk ('shpredeffect'). E.g., "~
  trt:subgroup1".

- unshrunk_prognostic_formula_str:

  Prognostic terms not to be shrunk ('unprogeffect'). E.g., "~ age +
  sex".

- unshrunk_predictive_formula_str:

  Predictive terms not to be shrunk ('unpredeffect'). E.g., "~
  trt:important_subgroup".

- shrunk_prognostic_formula_str:

  Prognostic terms to be shrunk ('shprogeffect'). E.g., "~ region +
  center".

- stratification_formula_str:

  A formula string specifying a stratification variable, e.g., "~
  strata_var".

- predictive_effect_priors:

  A named list with elements \`shrunk\` and/or \`unshrunk\` containing
  the priors for predictive effects. Can be strings or \`brmsprior\`
  objects. E.g., \`list(shrunk = "horseshoe(1)", unshrunk = "normal(0,
  5)")\`.

- prognostic_effect_priors:

  A named list with elements \`shrunk\`, \`unshrunk\` and/or
  \`intercept\` containing the priors for prognostic effects. E.g.,
  \`list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 10)")\`.

- stanvars:

  An object created by \`brms::stanvar()\` to add custom Stan code,
  necessary for some hierarchical priors.

- ...:

  Additional arguments passed directly to \`brms::brm()\` (e.g.,
  \`chains\`, \`iter\`, \`cores\`, \`backend\`).

## Value

A fitted \`brmsfit\` object.

## Details

This function is the main user-facing entry point. It first calls
\`prepare_formula_model\` to build the \`brmsformula\` and process the
data, then passes the results to \`fit_brms_model\` to run the analysis.

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

  # 2. Run the full analysis
  # We use \donttest{} because fitting a model can be slow
  # and is not suitable for automated CRAN checks.
  # \donttest{
  full_fit <- run_brms_analysis(
    data = sim_data,
    response_formula_str = "Surv(time, status) ~ trt",
    response_type = "survival",
    shrunk_predictive_formula_str = "~ trt:subgroup",
    unshrunk_prognostic_formula_str = "~ age",
    shrunk_prognostic_formula_str = "~ region",
    stratification_formula_str = "~ region",
    prognostic_effect_priors = list(
      shrunk = "normal(0, 1)",
      unshrunk = "normal(0, 5)"
    ),
    predictive_effect_priors = list(
      shrunk = "horseshoe(1)"
    ),
    chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
  )

  print(full_fit)
  # }
}
#> Step 1: Preparing formula and data...
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'region'.
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: subgroup
#> 
#> Step 2: Fitting the brms model...
#> Fitting brms model...
#> Compiling Stan program...
#> Start sampling
#> Warning: There were 16 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 2.11, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> 
#> Analysis complete.
#> Warning: Parts of the model have not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations and/or setting stronger priors.
#> Warning: There were 16 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: cox 
#>   Links: mu = log 
#> Formula: time | cens(1 - status) + bhaz(Boundary.knots = c(0.02, 99.98), knots = c(24, 46, 69), intercept = FALSE, gr = region) ~ unprogeffect + shprogeffect + shpredeffect 
#>          unprogeffect ~ age + trt + subgroup + 0
#>          shprogeffect ~ region + 0
#>          shpredeffect ~ trt_subgroupS1 + trt_subgroupS2 + trt_subgroupS3 + 0
#>    Data: data (Number of observations: 100) 
#>   Draws: 1 chains, each with iter = 50; warmup = 10; thin = 1;
#>          total post-warmup draws = 40
#> 
#> Regression Coefficients:
#>                             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
#> unprogeffect_age                0.01      0.00     0.00     0.02 1.45        2
#> unprogeffect_trt0               0.23      0.10     0.07     0.40 1.09        6
#> unprogeffect_trt1               0.29      0.31    -0.20     0.73 2.13        2
#> unprogeffect_subgroupS2        -0.44      0.19    -0.71    -0.05 2.00        7
#> unprogeffect_subgroupS3         0.28      0.18     0.01     0.57 1.47        2
#> shprogeffect_regionA           -0.13      0.06    -0.22    -0.03 1.33        3
#> shprogeffect_regionB           -0.62      0.16    -0.76    -0.22 1.34       18
#> shpredeffect_trt_subgroupS1    -0.39      0.17    -0.73    -0.15 1.42       11
#> shpredeffect_trt_subgroupS2    -0.41      0.41    -1.12     0.45 1.00       11
#> shpredeffect_trt_subgroupS3    -0.39      0.18    -0.62    -0.02 1.12        8
#>                             Tail_ESS
#> unprogeffect_age                  17
#> unprogeffect_trt0                 17
#> unprogeffect_trt1                 19
#> unprogeffect_subgroupS2           19
#> unprogeffect_subgroupS3           21
#> shprogeffect_regionA              18
#> shprogeffect_regionB              21
#> shpredeffect_trt_subgroupS1       17
#> shpredeffect_trt_subgroupS2       19
#> shpredeffect_trt_subgroupS3       45
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```
