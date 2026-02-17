# Quickstart

## 1 Introduction

The `bonsaiforest2` package consists of 3 core functions which are
typically called in sequence:

1.  [`run_brms_analysis()`](https://openpharma.github.io/bonsaiforest2/reference/run_brms_analysis.md) -
    Prepares the model formula and fits the Bayesian model using `brms`.
2.  [`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md) -
    Calculates the marginal subgroup treatment effects.
3.  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) - Creates a
    forest plot from the summary object.

The package enables the implementation of both **global shrinkage
models** ([Wolbers et al. 2025](#ref-wolbers2025using)), that estimate
all prognostic and predictive effects in a single model, and **one-way
shrinkage models** ([Wang et al. 2024](#ref-wang2024bayesian)), that
estimate the predictive effects using a different model for each
subgrouping variable.

This vignette demonstrates how to use the package to fit and compare
these different modeling formulas. You’ll learn how to:

- Fit one-way shrinkage models (one model per subgroup variable)
- Fit a global shrinkage model (all subgroup variables in one model)
- Generate summaries of subgroup treatment effects
- Visualize and compare results from different model specifications

This example makes use of Bayesian modeling, which requires the
installation of the [`brms`](https://paulbuerkner.com/brms/) package and
a working Stan installation (e.g., via
[`cmdstanr`](https://mc-stan.org/cmdstanr/)).

## 2 The Data

We will use a simulated example dataset representing a clinical trial
for blood pressure. The relevant endpoint is the change in Systolic
Blood Pressure (SBP) from baseline (`sbp_change`).

We consider a model where we want to find the treatment effect (`trt`)
on `sbp_change`. The model will adjust for `baseline_sbp` as a
**prognostic** variable (predictor of the outcome) and explore
**multiple subgroup variables** as **predictive** variables (potential
treatment effect modifiers):

- `region`: Geographic region (USA, EU, APAC)
- `comorbidity`: Presence of comorbidities (Yes, No)
- `age_group`: Age category (\< 50, 50-65, \> 65)
- `sex`: Biological sex (M, F)
- `diabetes`: Diabetes status (Yes, No)

First, let’s load the libraries and create the data.

``` r
# Load the main package
library(bonsaiforest2)
library(brms)
```

``` r
# Create the example data with multiple subgroup variables
set.seed(123)
n_patients <- 300

continuous_data <- data.frame(
  id = 1:n_patients,
  sbp_change = rnorm(n_patients, mean = -5, sd = 10),
  trt = sample(0:1, n_patients, replace = TRUE),
  baseline_sbp = rnorm(n_patients, mean = 140, sd = 15),
  region = factor(sample(c("USA", "EU", "APAC"), n_patients, replace = TRUE)),
  comorbidity = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.4, 0.6))),
  age_group = factor(sample(c("<50", "50-65", ">65"), n_patients, replace = TRUE, prob = c(0.3, 0.4, 0.3))),
  sex = factor(sample(c("M", "F"), n_patients, replace = TRUE)),
  diabetes = factor(sample(c("Yes", "No"), n_patients, replace = TRUE, prob = c(0.3, 0.7)))
)

continuous_data$trt <- factor(continuous_data$trt, levels = c(0, 1))

print(head(continuous_data))
#>   id sbp_change trt baseline_sbp region comorbidity age_group sex diabetes
#> 1  1 -10.604756   1     161.4560   APAC          No     50-65   M      Yes
#> 2  2  -7.301775   1     155.6994     EU          No       <50   F       No
#> 3  3  10.587083   1     146.5293     EU          No       <50   F       No
#> 4  4  -4.294916   0     150.7277     EU          No       <50   F       No
#> 5  5  -3.707123   0     153.7576   APAC         Yes     50-65   M       No
#> 6  6  12.150650   0     100.0862     EU          No     50-65   F      Yes
```

## 3 One-way Shrinkage Models

This section demonstrates how to fit separate models, each examining one
subgroup variable at a time. In one-way shrinkage models, we specify
shrunk predictive effects using random effects notation (e.g.,
`~ (1 + trt || subgroup)`), which allows treatment effects to vary by
subgroup with automatic hierarchical regularization.

### 3.1 One-way Model: Region Only

``` r
# Fit model with only region as subgroup variable using one-way shrinkage approach
# Random effects notation (||) estimates varying treatment effects by region
oneway_region <- run_brms_analysis(
  data = continuous_data,
  response_type = "continuous",
  response_formula = sbp_change ~ trt,
  unshrunk_terms_formula = ~ baseline_sbp,
  shrunk_predictive_formula = ~ (1 + trt || region),
  chains = 2, iter = 2000, warmup = 1000, cores = 2,
  refresh = 0, backend = "cmdstanr"
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 2 finished in 50.7 seconds.
#> Chain 1 finished in 53.2 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 52.0 seconds.
#> Total execution time: 53.4 seconds.

summary_oneway_region <- summary_subgroup_effects(brms_fit = oneway_region)
print(summary_oneway_region)
#> $estimates
#> # A tibble: 3 × 4
#>   Subgroup     Median CI_Lower CI_Upper
#>   <chr>         <dbl>    <dbl>    <dbl>
#> 1 region: APAC  1.30    -0.892     4.21
#> 2 region: EU    1.36    -0.888     4.27
#> 3 region: USA   0.944   -1.89      3.70
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

### 3.2 One-way Model: Comorbidity Only

``` r
oneway_comorbidity <- run_brms_analysis(
  data = continuous_data,
  response_type = "continuous",
  response_formula = sbp_change ~ trt,
  unshrunk_terms_formula = ~ baseline_sbp,
  shrunk_predictive_formula = ~ (1 + trt || comorbidity),
  chains = 2, iter = 2000, warmup = 1000, cores = 2,
  refresh = 0, backend = "cmdstanr"
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 finished in 52.9 seconds.
#> Chain 2 finished in 53.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 53.1 seconds.
#> Total execution time: 53.4 seconds.

summary_oneway_comorbidity <- summary_subgroup_effects(brms_fit = oneway_comorbidity)
```

### 3.3 One-way Model: Age Group Only

``` r
oneway_age <- run_brms_analysis(
  data = continuous_data,
  response_type = "continuous",
  response_formula = sbp_change ~ trt,
  unshrunk_terms_formula = ~ baseline_sbp,
  shrunk_predictive_formula = ~ (1 + trt || age_group),
  chains = 2, iter = 2000, warmup = 1000, cores = 2,
  refresh = 0, backend = "cmdstanr"
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 finished in 53.5 seconds.
#> Chain 2 finished in 53.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 53.6 seconds.
#> Total execution time: 53.8 seconds.

summary_oneway_age <- summary_subgroup_effects(brms_fit = oneway_age)
```

### 3.4 One-way Model: Sex Only

``` r
oneway_sex <- run_brms_analysis(
  data = continuous_data,
  response_type = "continuous",
  response_formula = sbp_change ~ trt,
  unshrunk_terms_formula = ~ baseline_sbp,
  shrunk_predictive_formula = ~ (1 + trt || sex),
  chains = 2, iter = 2000, warmup = 1000, cores = 2,
  refresh = 0, backend = "cmdstanr"
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 2 finished in 51.9 seconds.
#> Chain 1 finished in 52.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 52.3 seconds.
#> Total execution time: 52.7 seconds.

summary_oneway_sex <- summary_subgroup_effects(brms_fit = oneway_sex)
```

### 3.5 One-way Model: Diabetes Only

``` r
oneway_diabetes <- run_brms_analysis(
  data = continuous_data,
  response_type = "continuous",
  response_formula = sbp_change ~ trt,
  unshrunk_terms_formula = ~ baseline_sbp,
  shrunk_predictive_formula = ~ (1 + trt || diabetes),
  chains = 2, iter = 2000, warmup = 1000, cores = 2,
  refresh = 0, backend = "cmdstanr"
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 2 finished in 50.9 seconds.
#> Chain 1 finished in 52.4 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 51.6 seconds.
#> Total execution time: 52.5 seconds.

summary_oneway_diabetes <- summary_subgroup_effects(brms_fit = oneway_diabetes)
```

### 3.6 One-way Models: Visualizing All Models

You can combine and visualize results from multiple models using
[`combine_summaries()`](https://openpharma.github.io/bonsaiforest2/reference/combine_summaries.md):

``` r
# Combine all one-way shrinkage models
combined_oneway <- combine_summaries(list(
  "Region" = summary_oneway_region,
  "Comorbidity" = summary_oneway_comorbidity,
  "Age Group" = summary_oneway_age,
  "Sex" = summary_oneway_sex,
  "Diabetes" = summary_oneway_diabetes
))

plot(combined_oneway, title = "One-way Shrinkage Models: All Subgroup Variables")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Quickstart_files/figure-html/compare-all-oneway-1.png)

## 4 Global Shrinkage Model

This section demonstrates how to fit a single model that includes all
subgroup variables simultaneously using the global shrinkage approach.
All treatment-by-subgroup interactions are estimated in one unified
model with strong regularization (Horseshoe prior) applied to the
interaction terms.

### 4.1 Global Shrinkage Model: All Subgroups

``` r
# Fit a single model with ALL subgroup variables simultaneously
# - Unshrunk terms: baseline_sbp with weak priors (will be reference-coded)
# - Shrunk prognostic effects: subgroup main effects with strong regularization (one-hot coded)
# - Shrunk predictive effects: treatment interactions with strong regularization (one-hot coded)
global_shrinkage_model <- run_brms_analysis(
  data = continuous_data,
  response_type = "continuous",
  response_formula = sbp_change ~ trt,
  unshrunk_terms_formula = ~ baseline_sbp,
  shrunk_prognostic_formula = ~ 0 + region + comorbidity + age_group + sex + diabetes,
  shrunk_predictive_formula = ~ 0 + trt:region + trt:comorbidity + trt:age_group + trt:sex + trt:diabetes,
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 0.5)",
  chains = 2, iter = 2000, warmup = 1000, cores = 2,
  refresh = 0, backend = "cmdstanr"
)
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 finished in 5.8 seconds.
#> Chain 2 finished in 6.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 6.1 seconds.
#> Total execution time: 6.5 seconds.
```

### 4.2 Global Shrinkage Model: Summary of Subgroup Effects

Use
[`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md)
to generate marginal treatment effects for each subgroup. The function
automatically extracts all treatment interactions from the fitted model:

``` r
global_summary <- summary_subgroup_effects(brms_fit = global_shrinkage_model)
#> Using trt_var from model attributes: trt
#> Using response_type from model attributes: continuous
#> --- Calculating specific subgroup effects... ---
#> Using data from model attributes
#> Step 1: Identifying subgroups and creating counterfactuals...
#> `subgroup_vars` set to 'auto'. Detecting from model...
#> Model data has 300 rows and 9 columns
#> Column names: id, sbp_change, trt, baseline_sbp, region, comorbidity, age_group, sex, diabetes
#> Treatment variable: 'trt'
#> All coefficient names:
#> unshrunktermeffect_Intercept
#> unshrunktermeffect_baseline_sbp
#> unshrunktermeffect_trt
#> shprogeffect_regionAPAC
#> shprogeffect_regionEU
#> shprogeffect_regionUSA
#> shprogeffect_comorbidityNo
#> shprogeffect_comorbidityYes
#> shprogeffect_age_group<50
#> shprogeffect_age_group>65
#> shprogeffect_age_group50M65
#> shprogeffect_sexF
#> shprogeffect_sexM
#> shprogeffect_diabetesNo
#> shprogeffect_diabetesYes
#> shpredeffect_trt:regionAPAC
#> shpredeffect_trt:regionEU
#> shpredeffect_trt:regionUSA
#> shpredeffect_trt:comorbidityNo
#> shpredeffect_trt:comorbidityYes
#> shpredeffect_trt:age_group<50
#> shpredeffect_trt:age_group>65
#> shpredeffect_trt:age_group50M65
#> shpredeffect_trt:sexF
#> shpredeffect_trt:sexM
#> shpredeffect_trt:diabetesNo
#> shpredeffect_trt:diabetesYes
#> Looking for treatment interactions with pattern: 'trt:'
#> Found 12 treatment interaction coefficients
#> Treatment interaction coefficients found:
#> shpredeffect_trt:regionAPAC
#> shpredeffect_trt:regionEU
#> shpredeffect_trt:regionUSA
#> shpredeffect_trt:comorbidityNo
#> shpredeffect_trt:comorbidityYes
#> shpredeffect_trt:age_group<50
#> shpredeffect_trt:age_group>65
#> shpredeffect_trt:age_group50M65
#> shpredeffect_trt:sexF
#> shpredeffect_trt:sexM
#> shpredeffect_trt:diabetesNo
#> shpredeffect_trt:diabetesYes
#> Detected subgroup variable 'region' from coefficient 'shpredeffect_trt:regionAPAC'
#> Detected subgroup variable 'region' from coefficient 'shpredeffect_trt:regionEU'
#> Detected subgroup variable 'region' from coefficient 'shpredeffect_trt:regionUSA'
#> Detected subgroup variable 'comorbidity' from coefficient 'shpredeffect_trt:comorbidityNo'
#> Detected subgroup variable 'comorbidity' from coefficient 'shpredeffect_trt:comorbidityYes'
#> Detected subgroup variable 'age_group' from coefficient 'shpredeffect_trt:age_group<50'
#> Detected subgroup variable 'age_group' from coefficient 'shpredeffect_trt:age_group>65'
#> Detected subgroup variable 'age_group' from coefficient 'shpredeffect_trt:age_group50M65'
#> Detected subgroup variable 'sex' from coefficient 'shpredeffect_trt:sexF'
#> Detected subgroup variable 'sex' from coefficient 'shpredeffect_trt:sexM'
#> Detected subgroup variable 'diabetes' from coefficient 'shpredeffect_trt:diabetesNo'
#> Detected subgroup variable 'diabetes' from coefficient 'shpredeffect_trt:diabetesYes'
#> Checking for random effects parameters...
#> Retrieved 58 total parameters from model
#> Using regex pattern: '^r_(.+)__[^\[]+\[[^,]+,trt\]'
#> Found 0 matching random effect parameters
#> No random effect parameters matching the pattern were found
#> ...detected subgroup variable(s): region, comorbidity, age_group, sex, diabetes
#> Step 2: Generating posterior predictions...
#> ... detected Fixed Effects (Colon model). Predicting with re_formula = NA.
#> ... (predicting expected outcomes)...
#> Step 3: Calculating marginal effects...
#> ... processing region
#> ... processing comorbidity
#> ... processing age_group
#> ... processing sex
#> ... processing diabetes
#> Done.

# Print the summary of subgroup-specific treatment effects
print(global_summary)
#> $estimates
#> # A tibble: 12 × 4
#>    Subgroup         Median CI_Lower CI_Upper
#>    <chr>             <dbl>    <dbl>    <dbl>
#>  1 region: APAC      1.32    -1.47      4.20
#>  2 region: EU        0.934   -1.63      3.76
#>  3 region: USA       0.503   -2.32      3.03
#>  4 comorbidity: No   0.624   -1.80      3.00
#>  5 comorbidity: Yes  1.38    -1.34      4.18
#>  6 age_group: <50    1.46    -1.40      4.53
#>  7 age_group: >65    0.893   -2.04      3.61
#>  8 age_group: 50-65  0.616   -1.90      3.27
#>  9 sex: F            0.394   -2.19      2.85
#> 10 sex: M            1.38    -1.04      4.02
#> 11 diabetes: No      0.428   -2.02      2.83
#> 12 diabetes: Yes     1.81    -0.921     4.89
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

### 4.3 Global Shrinkage Model: Visualization

Use the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
function to create a forest plot from the summary object:

``` r
plot(global_summary, title = "Global Shrinkage Model: All Subgroup Variables")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Quickstart_files/figure-html/global-plot-1.png)

## 5 Comparing Multiple Models in One Plot

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) function
supports comparing multiple models side-by-side. Pass a named list of
`subgroup_summary` objects to create a comparative forest plot.

### 5.1 Example: Comparing One-way vs Global Shrinkage Models

``` r
# Combine summaries for comparison
combined <- combine_summaries(list(
  "One-way Shrinkage" = combined_oneway,
  "Global Shrinkage" = global_summary
))

# Plot the comparison
plot(combined, title = "Comparing One-way vs Global Shrinkage Models")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Quickstart_files/figure-html/compare-models-1.png)

Wang, Yun, Wenda Tu, William Koh, James Travis, Robert Abugov, Kiya
Hamilton, Mengjie Zheng, Roberto Crackel, Pablo Bonangelino, and Mark
Rothmann. 2024. “Bayesian hierarchical models for subgroup analysis.”
*Pharmaceutical Statistics* 23: 1065–83.

Wolbers, Marcel, Mar Vázquez Rabuñal, Ke Li, Kaspar Rufibach, and Daniel
Sabanés Bové. 2025. “Using shrinkage methods to estimate treatment
effects in overlapping subgroups in randomized clinical trials with a
time-to-event endpoint.” *Statistical Methods in Medical Research*,
1–12.
