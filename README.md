
# bonsaiforest2 <a href="https://github.com/openpharma/bonsaiforest2"><img src="man/figures/logo.png" alt="bonsaiforest2 website" align="right" height="138"/></a>

## Overview

The `bonsaiforest2` package is used for Bayesian shrinkage estimation of
subgroup treatment effects in randomized clinical trials. It supports
both One-Variable-at-a-Time (OVAT) and Global modeling approaches for
estimating treatment-by-subgroup interactions, with built-in support for
continuous, binary, time-to-event (Cox), and count outcomes. The package
implements state-of-the-art shrinkage priors including Regularized
Horseshoe and R2D2, combined with standardization (G-computation) to
provide interpretable marginal treatment effects. By leveraging `brms`
and `Stan`, `bonsaiforest2` provides a practical tool for obtaining more
stable and reliable subgroup effect estimates in exploratory analyses.

## Installation

**UPDATE TO USUAL INSTALLATION** You can install the development version
of `bonsaiforest2` from its GitLab repository:

``` r
# install.packages("remotes") 
remotes::install_github("openpharma/bonsaiforest2")
```

## Example

This example demonstrates the usage of `bonsaiforest2` for subgroup
treatment effect estimation across multiple overlapping subgroups (Age,
Region, and Biomarker) using a Global Modeling approach with a
Regularized Horseshoe prior.

``` r
library(bonsaiforest2)

# 1. Simulate trial data
set.seed(42)
n <- 200
trial_data <- data.frame(
  outcome  = rnorm(n),
  trt      = factor(sample(c("Control", "Active"), n, replace = TRUE)),
  age_cat  = factor(sample(c("<65", ">=65"), n, replace = TRUE)),
  region   = factor(sample(c("US", "EU", "Asia"), n, replace = TRUE)),
  biomarker = factor(sample(c("Pos", "Neg"), n, replace = TRUE))
)

# 2. Fit a Global Model with default priors
fit <- run_brms_analysis(
  data = trial_data,
  response_type = "continuous",
  response_formula = outcome ~ trt,
  unshrunk_terms_formula = ~ age_cat + region + biomarker, 
  shrunk_predictive_formula = ~ 0 + trt:age_cat + trt:region + trt:biomarker, 
  sigma_ref = sd(trial_data$outcome),  # Use observed sd or protocol value
  chains = 2, iter = 1000, warmup = 500 #
)

# 3. Derive Marginal Treatment Effects
subgroup_effects <- summary_subgroup_effects(
  brms_fit = fit
)
```

``` r
# Print summary
print(subgroup_effects)
#> $estimates
#> # A tibble: 7 × 4
#>   Subgroup         Median CI_Lower CI_Upper
#>   <chr>             <dbl>    <dbl>    <dbl>
#> 1 age_cat: <65    0.128    -0.241     0.516
#> 2 age_cat: >=65  -0.184    -0.566     0.219
#> 3 region: Asia   -0.109    -0.544     0.366
#> 4 region: EU     -0.349    -0.828     0.101
#> 5 region: US      0.387    -0.0988    0.918
#> 6 biomarker: Neg -0.0508   -0.426     0.258
#> 7 biomarker: Pos -0.00256  -0.326     0.346
#> 
#> $response_type
#> [1] "continuous"
#> 
#> $ci_level
#> [1] 0.95
#> 
#> $trt_var
#> [1] "trt"
#> 
#> attr(,"class")
#> [1] "subgroup_summary"

# 4. Visualize Results
plot(subgroup_effects)
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
