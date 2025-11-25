# Quickstart

## 1 Introduction

The `bonsaiforest2` package consists of 3 core functions which are
typically called in sequence:

1.  [`run_brms_analysis()`](https://openpharma.github.io/bonsaiforest2/reference/run_brms_analysis.md) -
    Prepares the model formula and fits the Bayesian model using `brms`.
2.  [`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md) -
    Calculates the marginal overall and subgroup treatment effects.
3.  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) - Creates a
    forest plot from the summary object.

This example makes use of Bayesian modeling, which requires the
installation of the [`brms`](https://paulbuerkner.com/brms/) package and
a working Stan installation (e.g., via
[`cmdstanr`](https://mc-stan.org/cmdstanr/)).

``` r
install.packages("brms")
install.packages("cmdstanr")
cmdstanr::install_cmdstan()
```

## 2 The Data

We will use a simulated example dataset representing a clinical trial
for blood pressure. The relevant endpoint is the change in Systolic
Blood Pressure (SBP) from baseline (`sbp_change`).

We consider an **analysis model** where we want to find the treatment
effect (`trt`) on `sbp_change`. The model will adjust for `baseline_sbp`
as a **prognostic** variable (predictor of the outcome) and explore
`region` and `comorbidity` as **predictive** variables (potential
treatment effect modifiers).

First, let’s load the libraries and create the data.

``` r
# Load the main package
library(bonsaiforest2)

# Load other required packages
library(brms)
library(dplyr)
```

``` r
# Create the example data
set.seed(123)
n_patients <- 200
continuous_data <- data.frame(
  id = 1:n_patients,
  sbp_change = rnorm(n_patients, mean = -5, sd = 10),
  trt = sample(0:1, n_patients, replace = TRUE),
  baseline_sbp = rnorm(n_patients, mean = 140, sd = 15),
  region = factor(sample(c("USA", "EU", "APAC"), n_patients, replace = TRUE)),
  comorbidity = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.4, 0.6)))
)

# `bonsaiforest2` expects the treatment variable to be a factor
continuous_data$trt <- factor(continuous_data$trt, levels = c(0, 1))

print(head(continuous_data))
#>   id sbp_change trt baseline_sbp region comorbidity
#> 1  1 -10.604756   0     129.2714    USA          No
#> 2  2  -7.301775   0     128.7097     EU         Yes
#> 3  3  10.587083   0     125.9219   APAC         Yes
#> 4  4  -4.294916   0     124.2123     EU          No
#> 5  5  -3.707123   0     133.4426   APAC         Yes
#> 6  6  12.150650   0     144.9677    USA          No
```

## 3 `run_brms_analysis()`

The
[`run_brms_analysis()`](https://openpharma.github.io/bonsaiforest2/reference/run_brms_analysis.md)
function is the first step. It builds the model formula and fits the
`brms` model, applying shrinkage priors to exploratory terms.

**Key Arguments:**

- `data`: Your input `data.frame`.
- `response_formula_str`: Defines the outcome and main treatment effect
  (e.g., `"outcome ~ trt"`).
- `response_type`: `"continuous"`, `"binary"`, `"count"`, or
  `"survival"`.
- `shrunk_predictive_formula_str`: Treatment interactions to estimate
  with shrinkage (e.g., `"~ region:trt + sex:trt"`).
- `unshrunk_prognostic_formula_str`: Main effects (not treatment
  interactions) to estimate without strong shrinkage (e.g., `"~ age"`).
- `shrunk_prognostic_formula_str`: Main effects to estimate with
  shrinkage.
- `unshrunk_predictive_formula_str`: Treatment interactions to estimate
  without strong shrinkage (e.g., `"~ prespecified_marker:trt"`).
- `prognostic_effect_priors`, `predictive_effect_priors`: Lists
  specifying priors for shrunk/unshrunk terms.
- `stanvars`: Optional object from
  [`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
  for custom Stan code (e.g., for hierarchical priors).
- `...`: Other arguments passed directly to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  (like `chains`, `iter`, `cores`).

For this example, we will:

- Define the response as `sbp_change ~ trt`.
- Adjust for `baseline_sbp` as an **unshrunk prognostic** variable.
- Explore `region` and `comorbidity` as **shrunk predictive** variables.
  We also include `~ 1` in the shrunk predictive formula to apply
  shrinkage to the overall treatment effect (intercept of the
  interactions).

``` r
# iter and chains are set low for a quick example.
# For a real analysis, use more iterations (e.g., iter = 2000).
continuous_model_fit <- run_brms_analysis(
  data = continuous_data,
  response_formula_str = "sbp_change ~ trt",
  response_type = "continuous",
  unshrunk_prognostic_formula_str = "~ baseline_sbp",
  # Shrink intercept (1), region interaction, and comorbidity interaction
  shrunk_predictive_formula_str = "~ 1 + region:trt + comorbidity:trt",
  chains = 1, iter = 200, warmup = 100, cores = 1,
  refresh = 0, backend = "cmdstanr" # Use cmdstanr if available
)
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: There aren't enough warmup iterations to fit the 
#> Chain 1          three stages of adaptation as currently configured. 
#> Chain 1          Reducing each adaptation stage to 15%/75%/10% of 
#> Chain 1          the given number of warmup iterations: 
#> Chain 1            init_buffer = 15 
#> Chain 1            adapt_window = 75 
#> Chain 1            term_buffer = 10 
#> Chain 1 finished in 2.3 seconds.

print(continuous_model_fit)
#>  Family: gaussian 
#>   Links: mu = identity 
#> Formula: sbp_change ~ unprogeffect + shpredeffect 
#>          unprogeffect ~ baseline_sbp + trt + region + comorbidity
#>          shpredeffect ~ trt_regionAPAC + trt_regionEU + trt_regionUSA + trt_comorbidityNo + trt_comorbidityYes + 0
#>    Data: data (Number of observations: 200) 
#>   Draws: 1 chains, each with iter = 200; warmup = 100; thin = 1;
#>          total post-warmup draws = 100
#> 
#> Regression Coefficients:
#>                                 Estimate Est.Error l-95% CI u-95% CI Rhat
#> unprogeffect_Intercept              1.41      4.83    -7.09    12.29 1.50
#> unprogeffect_baseline_sbp          -0.04      0.04    -0.11     0.03 1.29
#> unprogeffect_trt0                  -1.20      2.55    -6.11     4.20 1.00
#> unprogeffect_regionEU              -1.08      1.99    -5.42     2.79 1.00
#> unprogeffect_regionUSA             -1.74      1.97    -5.50     2.19 1.01
#> unprogeffect_comorbidityYes         1.26      1.79    -1.83     4.64 1.06
#> shpredeffect_trt_regionAPAC        -0.57      1.68    -3.93     3.26 1.02
#> shpredeffect_trt_regionEU           0.38      1.93    -3.03     4.64 1.03
#> shpredeffect_trt_regionUSA          0.55      1.45    -1.73     3.63 0.99
#> shpredeffect_trt_comorbidityNo      0.85      1.74    -1.47     5.28 1.01
#> shpredeffect_trt_comorbidityYes    -0.23      1.71    -3.57     3.07 0.99
#>                                 Bulk_ESS Tail_ESS
#> unprogeffect_Intercept                 2       32
#> unprogeffect_baseline_sbp              4       22
#> unprogeffect_trt0                     48       55
#> unprogeffect_regionEU                 49       71
#> unprogeffect_regionUSA                76       78
#> unprogeffect_comorbidityYes          104       52
#> shpredeffect_trt_regionAPAC           68       52
#> shpredeffect_trt_regionEU             74       60
#> shpredeffect_trt_regionUSA            68       81
#> shpredeffect_trt_comorbidityNo        63       78
#> shpredeffect_trt_comorbidityYes       98       45
#> 
#> Further Distributional Parameters:
#>       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
#> sigma     9.46      0.49     8.56    10.38 1.06      102       59
#> 
#> Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
#> and Tail_ESS are effective sample size measures, and Rhat is the potential
#> scale reduction factor on split chains (at convergence, Rhat = 1).
```

## 4 `summary_subgroup_effects()`

The next step is to use the fitted model to generate interpretable
subgroup effects. This is done with
[`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md).

This function uses **G-computation** to estimate the *marginal*
treatment effect for each subgroup (e.g., “what is the average effect
for all patients in the ‘EU’ region?”), averaging over all other
covariates like `baseline_sbp` and `comorbidity`.

Its main inputs are the `brms_fit` object from the previous step and the
`original_data`.

``` r
continuous_summary <- summary_subgroup_effects(
  brms_fit = continuous_model_fit,
  original_data = continuous_data, # Pass the original data
  trt_var = "trt",
  response_type = "continuous"     # Must match fitting
  # subgroup_vars = "auto" is the default and finds all interactions
)
#> --- Calculating specific subgroup effects... ---
#> Step 1: Creating counterfactual datasets...
#> `subgroup_vars` set to 'auto'. Detecting from model data...
#> ...detected subgroup variable(s): region, comorbidity
#> Step 2: Generating posterior predictions...
#> ... (predicting expected outcomes)...
#> Step 3: Calculating marginal effects...
#> ... processing region
#> ... processing comorbidity
#> Done.

print(continuous_summary)
#> $estimates
#> # A tibble: 5 × 4
#>   Subgroup         Median CI_Lower CI_Upper
#>   <chr>             <dbl>    <dbl>    <dbl>
#> 1 region: APAC       1.16   -2.46      4.76
#> 2 region: EU         1.81   -1.38      5.75
#> 3 region: USA        2.38   -1.05      6.09
#> 4 comorbidity: No    2.00   -0.525     4.52
#> 5 comorbidity: Yes   1.14   -2.39      4.61
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
```

## 5 `plot()`

Finally, the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
function takes the `subgroup_summary` object and creates a forest plot
for easy interpretation.

``` r
# Generate and display the plot
plot(continuous_summary, title = "Continuous: Subgroup Effects on SBP Change")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Quickstart_files/figure-html/ex-plot-1.png)

This plot displays the marginal treatment effect (mean difference in
`sbp_change`) for the overall population and for each subgroup level.
The estimates are “shrunk” towards the overall effect, providing more
stable results than fitting separate models for each subgroup.

## 6 Code

We report below all the code for the main workflow presented in this
vignette.

``` r
library(bonsaiforest2)
library(brms)
library(dplyr)
library(ggplot2)

# --- 2. The Data ---
set.seed(123)
n_patients <- 200
continuous_data <- data.frame(
  id = 1:n_patients,
  sbp_change = rnorm(n_patients, mean = -5, sd = 10),
  trt = sample(0:1, n_patients, replace = TRUE),
  baseline_sbp = rnorm(n_patients, mean = 140, sd = 15),
  region = factor(sample(c("NA", "EU", "APAC"), n_patients, replace = TRUE)),
  comorbidity = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.4, 0.6)))
)
continuous_data$trt <- factor(continuous_data$trt, levels = c(0, 1))


# --- 3. run_brms_analysis() ---
# (Use iter = 2000, chains = 4 for a real analysis)
continuous_model_fit <- run_brms_analysis(
  data = continuous_data,
  response_formula_str = "sbp_change ~ trt",
  response_type = "continuous",
  unshrunk_prognostic_formula_str = "~ baseline_sbp",
  shrunk_predictive_formula_str = "~ 1 + region:trt + comorbidity:trt",
  chains = 1, iter = 200, warmup = 100, cores = 1,
  refresh = 0, backend = "cmdstanr"
)


# --- 4. summary_subgroup_effects() ---
continuous_summary <- summary_subgroup_effects(
  brms_fit = continuous_model_fit,
  original_data = continuous_data,
  trt_var = "trt",
  response_type = "continuous"
)


# --- 5. plot() ---
plot(continuous_summary, title = "Continuous: Subgroup Effects on SBP Change")
```
