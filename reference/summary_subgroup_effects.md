# Create a Summary of Marginal Subgroup Treatment Effects

This function orchestrates calls to \`estimate_subgroup_effects\` to
generate a complete summary table. It calculates both the overall
marginal effect and the effects for specific subgroups.

## Usage

``` r
summary_subgroup_effects(
  brms_fit,
  original_data,
  trt_var,
  response_type = c("binary", "count", "continuous", "survival"),
  subgroup_vars = "auto"
)
```

## Arguments

- brms_fit:

  A fitted \`brmsfit\` object.

- original_data:

  The original \`data.frame\` used for model fitting (before processing
  by \`prepare_formula_model\`).

- trt_var:

  The name of the treatment variable (coded 0/1).

- response_type:

  The type of outcome variable. One of "binary", "count", "continuous",
  or "survival".

- subgroup_vars:

  A character vector of subgrouping variable names. Defaults to "auto",
  which detects subgroups from the model. Set to \`NULL\` to get only
  the overall effect.

## Value

An object of class \`subgroup_summary\`, which is a list containing:

- estimates:

  A \`tibble\` combining the overall and subgroup-specific effect
  estimates.

- response_type:

  The specified response type.

- ci_level:

  The credible interval level (hard-coded to 0.95).

- trt_var:

  The name of the treatment variable.

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

  # 3. Get the summary
  # This function calls estimate_subgroup_effects internally
  eff_summary <- summary_subgroup_effects(
    brms_fit = full_fit,
    original_data = sim_data,
    trt_var = "trt",
    response_type = "survival",
    subgroup_vars = c("subgroup", "region")
  )

  print(eff_summary)
  # }
}
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
#> --- Calculating overall marginal effect... ---
#> Step 1: Creating counterfactual datasets...
#> ...setting interaction dummy variables for the 'all treatment' scenario.
#> Step 2: Generating posterior predictions...
#> ... (reconstructing baseline hazard and getting linear predictors)...
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Step 3: Calculating marginal effects...
#> Done.
#> 
#> --- Calculating specific subgroup effects... ---
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
#> # A tibble: 6 × 4
#>   Subgroup     Median CI_Lower CI_Upper
#>   <chr>         <dbl>    <dbl>    <dbl>
#> 1 Overall       1.02     0.978    1.04 
#> 2 subgroup: S1  1.01     0.976    1.04 
#> 3 subgroup: S2  0.927    0.893    0.944
#> 4 subgroup: S3  1.07     1.03     1.10 
#> 5 region: A     1.02     0.978    1.04 
#> 6 region: B     1.01     0.976    1.03 
#> 
#> $response_type
#> [1] "survival"
#> 
#> $ci_level
#> [1] 0.95
#> 
#> $trt_var
#> [1] "trt"
#> 
#> attr(,"class")
#> [1] "subgroup_summary"
```
