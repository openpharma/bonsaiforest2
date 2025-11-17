# Estimate Marginal Subgroup Treatment Effects

The function uses a counterfactual, marginal approach based on the
posterior predictive distribution. It averages over all other covariates
to provide robust estimates of subgroup-specific effects.

The estimation follows these steps:

1.  Creates two counterfactual datasets based on the model data: one
    where all subjects receive "treatment" and one where all receive
    "control".

2.  Generates posterior predictions (e.g., \`brms::posterior_epred\`)
    for both scenarios.

3.  Calculates the treatment effect (e.g., difference or ratio) for each
    subject and posterior draw.

4.  Averages these individual-level effects within each subgroup to
    produce the final marginal subgroup estimates.

## Usage

``` r
estimate_subgroup_effects(
  brms_fit,
  original_data,
  trt_var,
  subgroup_vars = "auto",
  response_type = c("continuous", "binary", "count", "survival"),
  ndraws = NULL
)
```

## Arguments

- brms_fit:

  A fitted \`brmsfit\` object from \`fit_brms_model()\` or
  \`run_brms_analysis()\`.

- original_data:

  A \`data.frame\` containing the original data before it was processed
  by \`prepare_formula_model()\`. This is essential for correctly
  identifying the original subgroup factor levels.

- trt_var:

  A character string specifying the name of the treatment variable.

- subgroup_vars:

  A character vector of subgroup variable names found in
  \`original_data\`. If set to \`"auto"\` (the default), the function
  attempts to automatically identify subgroup variables from the model's
  interaction terms.

- response_type:

  The type of outcome variable. One of "binary", "count", "continuous",
  or "survival". This determines the scale of the effect (e.g., risk
  difference, rate ratio).

- ndraws:

  An integer specifying the number of posterior draws to use. If
  \`NULL\` (default), all available draws are used.

## Value

A \`data.frame\` (tibble) where each row corresponds to a subgroup (or
the "Overall" effect), providing the estimated marginal effect and
posterior summaries (e.g., \`mean\`, \`sd\`, \`q2.5\`, \`q97.5\`).

## Details

This post-processing function estimates marginal treatment effects for
specified subgroups from a fitted \`brmsfit\` object.

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
  # \donttest{
  full_fit <- run_brms_analysis(
    data = sim_data,
    response_formula_str = "Surv(time, status) ~ trt",
    response_type = "survival",
    shrunk_predictive_formula_str = "~ trt:subgroup",
    unshrunk_prognostic_formula_str = "~ age",
    shrunk_prognostic_formula_str = "~ region",
    chains = 1, iter = 50, warmup = 10, refresh = 0
  )

  # 3. Estimate subgroup effects
  # Note: We pass the original `sim_data`, not the processed model data.
  effects <- estimate_subgroup_effects(
    brms_fit = full_fit,
    original_data = sim_data,
    trt_var = "trt",
    subgroup_vars = c("subgroup", "region"),
    response_type = "survival"
  )

  print(effects)
  # }
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
#> Step 1: Preparing formula and data...
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Treatment 'trt' added to unshrunk prognostic terms by default.
#> Auto-adding missing prognostic effect for interaction: subgroup
#> 
#> Step 2: Fitting the brms model...
#> Using default priors for unspecified effects:
#>   - shrunk prognostic (b): horseshoe(1)
#>   - unshrunk prognostic (b): normal(0, 5)
#>   - shrunk predictive (b): horseshoe(1)
#> Fitting brms model...
#> Compiling Stan program...
#> Start sampling
#> Warning: There were 30 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 2.3, indicating chains have not mixed.
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
#> Step 1: Creating counterfactual datasets...
#> ...setting interaction dummy variables for the 'all treatment' scenario.
#> Step 2: Generating posterior predictions...
#> ... (reconstructing baseline hazard and getting linear predictors)...
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Step 3: Calculating marginal effects...
#> ... processing subgroup
#> ... processing region
#> Done.
#> $estimates
#> # A tibble: 5 × 4
#>   Subgroup     Median CI_Lower CI_Upper
#>   <chr>         <dbl>    <dbl>    <dbl>
#> 1 subgroup: S1  1.01     0.976    1.04 
#> 2 subgroup: S2  0.927    0.893    0.944
#> 3 subgroup: S3  1.07     1.03     1.10 
#> 4 region: A     1.02     0.978    1.04 
#> 5 region: B     1.01     0.976    1.03 
#> 
#> $draws
#> # A tibble: 40 × 5
#>    `subgroup: S1` `subgroup: S2` `subgroup: S3` `region: A` `region: B`
#>             <dbl>          <dbl>          <dbl>       <dbl>       <dbl>
#>  1          0.978          0.895           1.04       0.984       0.980
#>  2          0.978          0.895           1.04       0.984       0.980
#>  3          0.987          0.903           1.05       0.990       0.987
#>  4          0.987          0.903           1.05       0.990       0.987
#>  5          0.987          0.903           1.05       0.990       0.987
#>  6          0.976          0.893           1.03       0.978       0.976
#>  7          0.976          0.893           1.03       0.978       0.976
#>  8          0.976          0.893           1.03       0.978       0.976
#>  9          0.976          0.893           1.03       0.978       0.976
#> 10          0.976          0.893           1.03       0.978       0.976
#> # ℹ 30 more rows
#> 
```
