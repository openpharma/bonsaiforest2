
# bonsaiforest2 <a href="https://github.com/openpharma/bonsaiforest2"><img src="man/figures/logo.png" alt="bonsaiforest2 website" align="right" height="138"/></a>

## Overview

The `bonsaiforest2` package is used for Bayesian shrinkage estimation of
subgroup treatment effects in randomized clinical trials. It supports
both one-way models (random effects) and fixed effects modeling
approaches for estimating treatment-by-subgroup interactions, with
built-in support for continuous, binary, time-to-event (Cox), and count
outcomes. The package allows the usage of state-of-the-art shrinkage
priors including Regularized Horseshoe and R2D2, combined with
standardization (G-computation) to provide interpretable marginal
treatment effects. By leveraging `brms` and `Stan`, `bonsaiforest2`
provides a practical tool for obtaining more stable and reliable
subgroup effect estimates in exploratory and confirmatory analyses.

## Installation

You can install the latest stable release of `bonsaiforest2` from
[GitHub](https://github.com/openpharma/bonsaiforest2):

``` r
# install.packages("remotes")
remotes::install_github("openpharma/bonsaiforest2@v0.1.0-beta")
```

Or install the development version from the `main` branch:

``` r
remotes::install_github("openpharma/bonsaiforest2")
```

## Quick Start Example

This example demonstrates Bayesian shrinkage estimation of treatment
effects across subgroups using a **fixed effects model** with a
Regularized Horseshoe prior, applied to the `shrink_data` package
dataset with three subgrouping variables.

``` r
library(bonsaiforest2)

# 1. Load the shrink_data package dataset
data(shrink_data)

# 2. Fit fixed effects model with heterogeneous shrinkage
fit_fixed <- run_brms_analysis(
  data = shrink_data,
  response_type = "continuous",
  response_formula = y ~ trt,
  unshrunk_terms_formula = ~ x1 + x2 + x3,
  shrunk_prognostic_formula = NULL,
  shrunk_predictive_formula = ~ 0 + trt:x1 + trt:x2 + trt:x3,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 5)",
  shrunk_prognostic_prior = NULL,
  shrunk_predictive_prior = "horseshoe(1)",
  chains = 2, iter = 1000, warmup = 500,
  backend = "cmdstanr", refresh = 0
)
```

``` r
# 3. Extract and visualize marginal treatment effects by subgroup
subgroup_effects <- summary_subgroup_effects(fit_fixed)
print(subgroup_effects$estimates)
#> # A tibble: 9 × 4
#>   Subgroup Median  CI_Lower CI_Upper
#>   <chr>     <dbl>     <dbl>    <dbl>
#> 1 x1: a     0.489  0.247       0.724
#> 2 x1: b     0.113 -0.257       0.458
#> 3 x2: a     0.388  0.0938      0.643
#> 4 x2: b     0.411  0.182       0.688
#> 5 x2: c     0.270  0.00723     0.508
#> 6 x3: a     0.343  0.0939      0.593
#> 7 x3: b     0.404  0.165       0.740
#> 8 x3: c     0.381  0.132       0.691
#> 9 x3: d     0.315 -0.000216    0.565
plot(subgroup_effects)
```

<img src="man/figures/README-example-summary-1.png" width="100%" />
