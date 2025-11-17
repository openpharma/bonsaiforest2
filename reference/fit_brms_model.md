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
  # \donttest{
  fit <- fit_brms_model(
    prepared_model = prepared_model,
    # Note: intercept prior is not needed for survival models
    prognostic_effect_priors = list(shrunk = "horseshoe(1)", unshrunk = "normal(0, 2)"),
    predictive_effect_priors = list(shrunk = "horseshoe(1)"),
    chains = 1, iter = 50, warmup = 10, refresh = 0 # For a quick example
  )

  print(fit)
 # }
}
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'region'.
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: subgroup
#> Fitting brms model...
#> Compiling Stan program...
#> Start sampling
#> Warning: There were 33 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 2.06, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Parts of the model have not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations and/or setting stronger priors.
#> Warning: There were 33 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#>  Family: cox 
#>   Links: mu = log 
#> Formula: time | cens(1 - status) + bhaz(Boundary.knots = c(0.02, 99.98), knots = c(24, 46, 69), intercept = FALSE, gr = region) ~ unprogeffect + shprogeffect + shpredeffect 
#>          unprogeffect ~ age + trt + subgroup + 0
#>          shprogeffect ~ region + 0
#>          shpredeffect ~ subgroup_S1_x_trt + subgroup_S2_x_trt + subgroup_S3_x_trt + 0
#>    Data: data (Number of observations: 100) 
#>   Draws: 1 chains, each with iter = 50; warmup = 10; thin = 1;
#>          total post-warmup draws = 40
#> 
#> Regression Coefficients:
#>                                Estimate Est.Error l-95% CI u-95% CI Rhat
#> unprogeffect_age                   0.01      0.00     0.01     0.02 1.89
#> unprogeffect_trt                  -0.07      0.01    -0.09    -0.06 1.50
#> unprogeffect_subgroupS1            0.75      0.00     0.75     0.76 1.94
#> unprogeffect_subgroupS2           -1.03      0.01    -1.04    -1.01 1.94
#> unprogeffect_subgroupS3           -0.53      0.00    -0.53    -0.53 1.50
#> shprogeffect_regionA              -0.01      0.00    -0.01     0.00 1.50
#> shprogeffect_regionB              -1.55      0.01    -1.57    -1.53 1.50
#> shpredeffect_subgroup_S1_x_trt     0.06      0.00     0.06     0.06 1.94
#> shpredeffect_subgroup_S2_x_trt    -0.03      0.00    -0.04    -0.03 1.32
#> shpredeffect_subgroup_S3_x_trt     0.13      0.01     0.11     0.13 1.94
#>                                Bulk_ESS Tail_ESS
#> unprogeffect_age                      9       NA
#> unprogeffect_trt                      2       NA
#> unprogeffect_subgroupS1               9       NA
#> unprogeffect_subgroupS2               2       NA
#> unprogeffect_subgroupS3               2       NA
#> shprogeffect_regionA                  2       NA
#> shprogeffect_regionB                 14        7
#> shpredeffect_subgroup_S1_x_trt        9       NA
#> shpredeffect_subgroup_S2_x_trt       14        7
#> shpredeffect_subgroup_S3_x_trt        2       NA
#> 
#> Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```
