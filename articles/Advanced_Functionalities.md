# Advanced Functionalities

This vignette demonstrates advanced functionalities for Bayesian
shrinkage estimation, including specifying an **offset for count data**,
customizing **prior distributions**, and using **stratification** to
handle heterogeneity in nuisance parameters.

## 1 Handling Exposure in Count Outcomes with Offsets

For count data (like disease exacerbations or event rates), the observed
count often depends on the **exposure time** or **follow-up time**.
Using an **offset** variable is the standard statistical method to
properly account for this time variation by modeling the event **rate**
(count per unit time) instead of the raw count.

In `bonsaiforest2`, you can include the offset directly in the
`response_formula` using the standard `brms` syntax:
`y + offset(log_exposure_time) ~ trt`.

### 1.1 Example 1: Count Outcome with Offset (Disease Exacerbations)

This scenario models exacerbation counts using a Negative Binomial
distribution and explicitly accounts for the patient’s exposure time.

*Scenario*: Modeling count outcomes with exposure-time adjustment.
Adjust for `baseline` (unshrunk prognostic effect) and multiple
exploratory `biomarkers` (shrunk prognostic effects). Overall treatment
effect on the rate scale.

``` r
# Data Simulation 
set.seed(789)
library(bonsaiforest2)
set.seed(789)
n <- 200

# Create biomarker data with simple naming
biomarker_data <- as.data.frame(matrix(rnorm(n * 8), ncol = 8))
names(biomarker_data) <- paste0("B", 1:8)

# Generate count data with treatment effect heterogeneity
count_data <- data.frame(
  y = rnbinom(n, size = 1.5, mu = 3),
  trt = rep(0:1, length.out = n),
  baseline = rnorm(n, 10, 2),
  log_exposure_time = log(runif(n, 0.5, 1.5))
)

# Add biomarkers and account for heterogeneous treatment effects
count_data <- cbind(count_data, biomarker_data)
count_data$y <- rnbinom(n, size = 1.5, mu = exp(0.5 + 0.1 * count_data$baseline + 
                        count_data$trt * (0.3 + 0.2 * (count_data$B1 > 0))))

# Create formula string for all biomarkers
shrunk_prog_str <- paste("~ 0 +", paste(names(biomarker_data), collapse = " + "))
```

``` r
# Model Fitting with Offset
count_model_fit <- run_brms_analysis(
  data = count_data,
  # Include offset(log_exposure_time) directly in the response formula
  response_formula = y + offset(log_exposure_time) ~ trt,
  response_type = "count",
  unshrunk_terms_formula = ~ baseline,
  shrunk_prognostic_formula = as.formula(shrunk_prog_str),
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 2.5)",
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 1)",
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Step 1: Preparing formula and data...
#> 
#> Step 2: Fitting the brms model...
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 2 finished in 2.8 seconds.
#> Chain 1 finished in 3.1 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 2.9 seconds.
#> Total execution time: 3.2 seconds.
#> Warning: 3 of 1000 (0.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
#> Loading required namespace: rstan
#> 
#> Analysis complete.
```

## 2 Customizing Prior Distributions

The choice of prior is central to Bayesian shrinkage. `bonsaiforest2`
provides sensible defaults, but it allows for full customization using
dedicated prior arguments: `intercept_prior`, `unshrunk_prior`,
`shrunk_prognostic_prior`, and `shrunk_predictive_prior`.

### 2.1 Prior Specification Mechanics

Priors are specified using separate parameters for each component.

| Prior Component       | Parameter Name            | Default Prior                 | Notes                                                                                                                                                    |
|:----------------------|:--------------------------|:------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Intercept**         | `intercept_prior`         | `NULL` (brms default)         | Optional custom prior                                                                                                                                    |
| **Unshrunk Terms**    | `unshrunk_prior`          | `NULL` (brms default)         | Optional custom prior                                                                                                                                    |
| **Shrunk Prognostic** | `shrunk_prognostic_prior` | `horseshoe(scale_global = 1)` | Strong shrinkage for fixed effects (colon syntax); For random effects (pipe-pipe syntax), automatically uses `normal(0, 1)` on SD scale                  |
| **Shrunk Predictive** | `shrunk_predictive_prior` | `horseshoe(scale_global = 1)` | Strong shrinkage for fixed effects (colon syntax); For random effects (pipe-pipe syntax), can be specified via `set_prior("normal(0, 1)", class = "sd")` |

### 2.2 Practical Examples of Prior Setting

The following examples demonstrate how to customize priors using the new
API. We use a synthetic survival dataset to show different prior
specification strategies, from simple defaults to advanced hierarchical
structures.

#### 2.2.1 Dataset Generation and Model Preparation

This dataset mimics a clinical trial with treatment effects that vary by
subgroups.

``` r
# Load library
library(bonsaiforest2)

# 1. Create Sample Data with treatment heterogeneity
set.seed(123)
n <- 200

sim_data <- data.frame(
  time = round(runif(n, 1, 100)),
  status = sample(0:1, n, replace = TRUE),
  trt = rep(0:1, length.out = n),
  age = rnorm(n, 50, 10),
  X1 = factor(sample(c("A", "B"), n, replace = TRUE)),
  X2 = factor(sample(c("A", "B", "C"), n, replace = TRUE))
)

# Ensure variables are factors
sim_data$trt <- factor(sim_data$trt, levels = c(0, 1))
sim_data$X1 <- as.factor(sim_data$X1)
sim_data$X2 <- as.factor(sim_data$X2)

# 2. Prepare the model formula
prepared_model <- prepare_formula_model(
  data = sim_data,
  response_formula = Surv(time, status) ~ trt,
  shrunk_predictive_formula = ~ 0 + trt:X1,
  unshrunk_terms_formula = ~ age,
  response_type = "survival"
)
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Note: Marginality principle not followed - interaction term 'X1' is used without its main effect. Consider adding 'X1' to prognostic terms for proper model hierarchy.
```

#### 2.2.2 Example 2: Using Default Priors

The simplest approach: use the package defaults that adapt to your
outcome type (horseshoe(1) for shrunk effects, brms defaults for
others).

``` r
fit_ex3 <- fit_brms_model(
  prepared_model = prepared_model,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 2.5)",
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 1)",
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 2 finished in 2.1 seconds.
#> Chain 1 finished in 2.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 2.2 seconds.
#> Total execution time: 2.4 seconds.

# View the priors that were automatically set
cat("\n=== Priors Used ===\n")
#> 
#> === Priors Used ===
print(fit_ex3[["prior"]])
#>                        prior class    coef group resp dpar              nlpar
#>  horseshoe(scale_global = 1)     b                               shpredeffect
#>  horseshoe(scale_global = 1)     b trt:X1A                       shpredeffect
#>  horseshoe(scale_global = 1)     b trt:X1B                       shpredeffect
#>               normal(0, 2.5)     b                         unshrunktermeffect
#>               normal(0, 2.5)     b     age                 unshrunktermeffect
#>               normal(0, 2.5)     b     trt                 unshrunktermeffect
#>                 dirichlet(1) sbhaz                                           
#>  lb ub tag       source
#>                    user
#>            (vectorized)
#>            (vectorized)
#>                    user
#>            (vectorized)
#>            (vectorized)
#>                 default
```

#### 2.2.3 Example 3: Using the R2D2 Shrinkage Prior

The R2D2 prior is useful when you want to control the *global* shrinkage
via the coefficient of determination (\\R^2\\) rather than a scale
parameter. This is often more interpretable for stakeholders.

``` r
# Use a custom R2D2 prior for the shrunk predictive effects (interactions)

fit_ex4 <- fit_brms_model(
  prepared_model = prepared_model,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 2.5)",
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "R2D2(mean_R2 = 0.5, prec_R2 = 1)",
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 2 finished in 1.6 seconds.
#> Chain 1 finished in 2.2 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 1.9 seconds.
#> Total execution time: 2.3 seconds.
#> Warning: 19 of 1000 (2.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.

cat("\n=== Priors Used ===\n")
#> 
#> === Priors Used ===
print(fit_ex4[["prior"]])
#>                             prior class    coef group resp dpar
#>  R2D2(mean_R2 = 0.5, prec_R2 = 1)     b                        
#>  R2D2(mean_R2 = 0.5, prec_R2 = 1)     b trt:X1A                
#>  R2D2(mean_R2 = 0.5, prec_R2 = 1)     b trt:X1B                
#>                    normal(0, 2.5)     b                        
#>                    normal(0, 2.5)     b     age                
#>                    normal(0, 2.5)     b     trt                
#>                      dirichlet(1) sbhaz                        
#>               nlpar lb ub tag       source
#>        shpredeffect                   user
#>        shpredeffect           (vectorized)
#>        shpredeffect           (vectorized)
#>  unshrunktermeffect                   user
#>  unshrunktermeffect           (vectorized)
#>  unshrunktermeffect           (vectorized)
#>                                    default
```

#### 2.2.4 Example 4: Custom Hierarchical Prior (Advanced)

This example demonstrates injecting raw Stan code using `stanvars`. This
is necessary if you want to implement a hierarchical structure that
`brms` does not support natively, such as estimating a shared variance
parameter across coefficients.

``` r
# 1. Define new hyperparameters in Stan
stanvars_full_hierarchical <- brms::stanvar(
  scode = "  real mu_pred;\n  real<lower=0> sigma_pred;\n",
  block = "parameters"
) +
  # Add priors for these parameters
  brms::stanvar(
    scode = "  // Priors on the hierarchical parameters\n  target += normal_lpdf(mu_pred | 0, 4); \n  target += normal_lpdf(sigma_pred | 0, 1) - normal_lccdf(0 | 0, 1); \n",
    block = "model"
  )

# 2. Create prior that references the Stan variables
prior_full_hierarchical <- brms::set_prior("normal(mu_pred, sigma_pred)")

# 3. Pass both to fit_brms_model
fit_ex5 <- fit_brms_model(
  prepared_model = prepared_model,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 2.5)",
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = prior_full_hierarchical,
  stanvars = stanvars_full_hierarchical,
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 finished in 2.9 seconds.
#> Chain 2 finished in 3.0 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 2.9 seconds.
#> Total execution time: 3.1 seconds.
#> Warning: 59 of 1000 (6.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.

# View the used priors
cat("\n=== Priors Used ===\n")
#> 
#> === Priors Used ===
print(fit_ex5[["prior"]])
#>                        prior class    coef group resp dpar              nlpar
#>  normal(mu_pred, sigma_pred)     b                               shpredeffect
#>  normal(mu_pred, sigma_pred)     b trt:X1A                       shpredeffect
#>  normal(mu_pred, sigma_pred)     b trt:X1B                       shpredeffect
#>               normal(0, 2.5)     b                         unshrunktermeffect
#>               normal(0, 2.5)     b     age                 unshrunktermeffect
#>               normal(0, 2.5)     b     trt                 unshrunktermeffect
#>                 dirichlet(1) sbhaz                                           
#>  lb ub tag       source
#>                    user
#>            (vectorized)
#>            (vectorized)
#>                    user
#>            (vectorized)
#>            (vectorized)
#>                 default
print(fit_ex5[["stanvars"]])
#> [[1]]
#> [[1]]$name
#> [1] ""
#> 
#> [[1]]$sdata
#> NULL
#> 
#> [[1]]$scode
#> [1] "  real mu_pred;\n  real<lower=0> sigma_pred;\n"
#> 
#> [[1]]$block
#> [1] "parameters"
#> 
#> [[1]]$position
#> [1] "start"
#> 
#> [[1]]$pll_args
#> character(0)
#> 
#> 
#> [[2]]
#> [[2]]$name
#> [1] ""
#> 
#> [[2]]$sdata
#> NULL
#> 
#> [[2]]$scode
#> [1] "  // Priors on the hierarchical parameters\n  target += normal_lpdf(mu_pred | 0, 4); \n  target += normal_lpdf(sigma_pred | 0, 1) - normal_lccdf(0 | 0, 1); \n"
#> 
#> [[2]]$block
#> [1] "model"
#> 
#> [[2]]$position
#> [1] "start"
#> 
#> [[2]]$pll_args
#> character(0)
#> 
#> 
#> attr(,"class")
#> [1] "stanvars"
```

#### 2.2.5 Example 5: Coefficient-Specific Priors

You can set different priors for specific coefficients by passing a
`brmsprior` object (created with
[`c()`](https://rdrr.io/r/base/c.html)).

*Scenario*: Set a general prior for all unshrunk terms, but use a
tighter prior specifically for treatment-biomarker interactions.

*Key Technique*: Use
[`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md)
to discover the exact coefficient names that will appear in the model.

``` r
# 1. Run prepare_formula_model
prepared_model <- prepare_formula_model(
  data = sim_data,
  response_formula = Surv(time, status) ~ trt,
  shrunk_predictive_formula = ~ 0 + trt:X1,
  unshrunk_terms_formula = ~ age + trt*X2,
  response_type = "survival"
)
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Note: Marginality principle not followed - interaction term 'X1' is used without its main effect. Consider adding 'X1' to prognostic terms for proper model hierarchy.

# 2. Inspect the Results
# A. The generated brms formula object
print(prepared_model$formula)
#> time | cens(1 - status) + bhaz(Boundary.knots = c(0.02, 99.98), knots = c(27, 51, 73), intercept = FALSE) ~ unshrunktermeffect + shpredeffect 
#> unshrunktermeffect ~ 0 + age + trt * X2 + trt
#> shpredeffect ~ 0 + trt:X1

# B. The processed terms
print(prepared_model$stan_variable_names$X_unshrunktermeffect)
#> [1] "age"     "trt"     "X2A"     "X2B"     "X2C"     "trt:X2B" "trt:X2C"

cat("=== Prior Strategy ===\n")
#> === Prior Strategy ===
cat("General unshrunk prior: normal(0, 5)\n")
#> General unshrunk prior: normal(0, 5)
cat("Specific for trt:X2 interactions: normal(0, 1)\n\n")
#> Specific for trt:X2 interactions: normal(0, 1)

# Create combined prior object
# IMPORTANT: Use EXACT coefficient names from the prepared data
# These names appear after contrast encoding in prepare_formula_model
unshrunk_priors_combined <- c(
  brms::set_prior("normal(0, 5)", class = "b"),  # General
  brms::set_prior("normal(0, 1)", class = "b", coef = "trt:X2B"),
  brms::set_prior("normal(0, 1)", class = "b", coef = "trt:X2C")
)

# Fit the model
fit_specific <- fit_brms_model(
  prepared_model = prepared_model,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = unshrunk_priors_combined,  # Pass the combined object
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 1)",
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 finished in 1.5 seconds.
#> Chain 2 finished in 2.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 2.1 seconds.
#> Total execution time: 2.8 seconds.

# View the used priors
cat("\n=== Priors Used ===\n")
#> 
#> === Priors Used ===
print(fit_specific[["prior"]])
#>                        prior class    coef group resp dpar              nlpar
#>  horseshoe(scale_global = 1)     b                               shpredeffect
#>  horseshoe(scale_global = 1)     b trt:X1A                       shpredeffect
#>  horseshoe(scale_global = 1)     b trt:X1B                       shpredeffect
#>                 normal(0, 5)     b                         unshrunktermeffect
#>                 normal(0, 5)     b     age                 unshrunktermeffect
#>                 normal(0, 5)     b     trt                 unshrunktermeffect
#>                 normal(0, 1)     b trt:X2B                 unshrunktermeffect
#>                 normal(0, 1)     b trt:X2C                 unshrunktermeffect
#>                 normal(0, 5)     b     X2A                 unshrunktermeffect
#>                 normal(0, 5)     b     X2B                 unshrunktermeffect
#>                 normal(0, 5)     b     X2C                 unshrunktermeffect
#>                 dirichlet(1) sbhaz                                           
#>  lb ub tag       source
#>                    user
#>            (vectorized)
#>            (vectorized)
#>                    user
#>            (vectorized)
#>            (vectorized)
#>                    user
#>                    user
#>            (vectorized)
#>            (vectorized)
#>            (vectorized)
#>                 default
```

#### 2.2.6 Example 6: Hierarchical Priors with Shared Variance (stanvars)

Create correlated priors by sharing a common variance parameter
estimated from the data.

*Use case*: When you believe treatment effects across subgroups should
be exchangeable, you can pool information by giving them a shared
variance.

``` r
cat("\n=== Hierarchical Prior Structure ===\n")
#> 
#> === Hierarchical Prior Structure ===
cat("Instead of independent priors per coefficient, we pool information through a shared scale:\n")
#> Instead of independent priors per coefficient, we pool information through a shared scale:
cat("  tau ~ half-normal(0, 1)  [shared scale parameter]\n")
#>   tau ~ half-normal(0, 1)  [shared scale parameter]
cat("  beta_i ~ N(0, tau) for each coefficient i\n")
#>   beta_i ~ N(0, tau) for each coefficient i
cat("This creates exchangeability: coefficients are ~similar but with adaptive variation.\n\n")
#> This creates exchangeability: coefficients are ~similar but with adaptive variation.

# Step 1: Declare tau as a parameter to be estimated
tau_parameter <- brms::stanvar(
  scode = "  real<lower=0> biomarker_tau;  // Shared scale parameter\n",
  block = "parameters"
)

# Step 2: Add prior for tau (using normal truncated to positive values with constraint)
tau_prior <- brms::stanvar(
  scode = "  biomarker_tau ~ normal(0, 1);  // Hyperprior on the scale\n",
  block = "model"
)

# Combine stanvars
hierarchical_stanvars <- tau_parameter + tau_prior

# Step 3: Create priors referencing the shared variance parameter
# Note: We identified these coefficient names using prepare_formula_model above
unshrunk_priors_hier <- c(
  brms::set_prior("normal(0, 5)", class = "b"),  # General
  brms::set_prior("normal(0, biomarker_tau)", class = "b", coef = "trt:X2B"),
  brms::set_prior("normal(0, biomarker_tau)", class = "b", coef = "trt:X2C")
)
```

``` r
# Step 4: Fit the hierarchical model
fit_hier <- fit_brms_model(
  prepared_model = prepared_model,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = unshrunk_priors_hier,
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 1)",
  stanvars = hierarchical_stanvars,
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)

# View the used priors
cat("\n=== Priors Used ===\n")
print(fit_hier[["prior"]])
print(fit_hier[["stanvars"]])
```

## 3 Advanced Model Parameterization

This section demonstrates advanced features for controlling how your
model represents subgroups.

### 3.1 Example 7: Custom Contrast Encoding for Shrunk Terms

By default, `bonsaiforest2` uses one-hot encoding (all factor levels, no
reference category) for shrunk interaction terms specified with
`~ 0 + ...` syntax, enabling full hierarchical shrinkage of all levels.
For unshrunk terms, treatment contrasts (k-1 dummy variables with a
reference level) are used by default. However, you can manually set
custom contrasts in your data using
[`contrasts()`](https://rdrr.io/r/stats/contrasts.html), and the library
will preserve your choice throughout the analysis.

*Scenario*: You want to use custom contrast coding for subgroup
variables, such as Helmert contrasts or sum contrasts, instead of the
default treatment contrasts. This is useful for specific hypothesis
structures or when comparing non-orthogonal effects.

*Key Technique*: Set contrasts using
[`contrasts()`](https://rdrr.io/r/stats/contrasts.html) before calling
[`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md).
The library will preserve your choice.

``` r
# Create sample data with multiple subgroups
set.seed(456)
n <- 250
sample_data_contrast <- data.frame(
  id = 1:n,
  trt = rep(0:1, length.out = n),
  y = rnorm(n, 50, 10),
  baseline = rnorm(n, 65, 10),
  X1 = factor(sample(c("A", "B", "C"), n, replace = TRUE))
)

# Set Helmert contrasts for the X1 variable
contrasts(sample_data_contrast$X1) <- contr.helmert(3)
cat("Contrast matrix for X1:\n")
#> Contrast matrix for X1:
print(contrasts(sample_data_contrast$X1))
#>   [,1] [,2]
#> A   -1   -1
#> B    1   -1
#> C    0    2

# Prepare model 
prepared_custom_contrast <- prepare_formula_model(
  data = sample_data_contrast,
  response_formula = y ~ trt,
  unshrunk_terms_formula = ~ baseline,
  shrunk_predictive_formula = ~ 0 + trt:X1, 
  response_type = "continuous"
)
#> Note: Marginality principle not followed - interaction term 'X1' is used without its main effect. Consider adding 'X1' to prognostic terms for proper model hierarchy.

# Observe costum contrast in the data
str(prepared_custom_contrast$data)
#> 'data.frame':    250 obs. of  5 variables:
#>  $ id      : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ trt     : int  0 1 0 1 0 1 0 1 0 1 ...
#>  $ y       : num  36.6 56.2 58 36.1 42.9 ...
#>  $ baseline: num  62.5 65.7 67.5 52.7 64.4 ...
#>  $ X1      : Factor w/ 3 levels "A","B","C": 2 3 2 3 1 2 3 3 2 3 ...
#>   ..- attr(*, "contrasts")= num [1:3, 1:2] -1 1 0 -1 -1 2
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:3] "A" "B" "C"
#>   .. .. ..$ : NULL

# Fit the model
fit_custom_contrast <- fit_brms_model(
  prepared_model = prepared_custom_contrast,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 2.5)",
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 1)",
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 finished in 1.2 seconds.
#> Chain 2 finished in 1.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 1.3 seconds.
#> Total execution time: 1.4 seconds.

estimate_custom_contrast <- summary_subgroup_effects(fit_custom_contrast)
#> --- Calculating specific subgroup effects... ---
#> Step 1: Identifying subgroups and creating counterfactuals...
#> ...detected subgroup variable(s): X1
#> Step 2: Generating posterior predictions...
#> Step 3: Calculating marginal effects...
#> Done.

print(estimate_custom_contrast)
#> $estimates
#> # A tibble: 3 × 4
#>   Subgroup Median CI_Lower CI_Upper
#>   <chr>     <dbl>    <dbl>    <dbl>
#> 1 X1: A      2.15   -0.778     4.74
#> 2 X1: B      2.35   -0.697     5.41
#> 3 X1: C      1.95   -0.948     4.67
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

### 3.2 Example 8: Prediction on New Data

After fitting a model, you can use `brms::predict()` to generate
predictions for new patients.

*Scenario*: We’ve fitted a treatment effect model and want to predict
outcomes for new patients with different characteristics.

``` r
# Use the fitted model from Example 7
# Generate new data for prediction
set.seed(789)
new_patients <- data.frame(
  id = 1:20,
  trt = rep(0:1, length.out = 20),
  baseline = rnorm(20, 65, 10),
  X1 = factor(rep(c("A", "B", "C"), length.out = 20))
)

# IMPORTANT: Set the same contrasts as in Example 7
# The prediction data must use the same encoding as the training data
contrasts(new_patients$X1) <- contr.helmert(3)

# Prepare new data
prepared_custom_contrast <- prepare_formula_model(
  data = new_patients,
  response_formula = y ~ trt,
  unshrunk_terms_formula = ~ baseline,
  shrunk_predictive_formula = ~ 0 + trt:X1, 
  response_type = "continuous"
)
#> Note: Marginality principle not followed - interaction term 'X1' is used without its main effect. Consider adding 'X1' to prognostic terms for proper model hierarchy.
#> Warning: Could not extract Stan variable names. Error: The following variables can neither be found in 'data' nor in 'data2':
#> 'y'

# Generate predictions using brms::predict
# This returns posterior predictive samples
predictions <- predict(fit_custom_contrast, newdata = prepared_custom_contrast$data, summary = TRUE)

head(predictions)
#>      Estimate Est.Error     Q2.5    Q97.5
#> [1,] 50.21508  10.44761 30.06733 70.44563
#> [2,] 42.31451  10.07955 22.91911 62.45141
#> [3,] 48.14653  10.12869 28.95442 68.88111
#> [4,] 51.22246  10.65935 31.23163 71.00125
#> [5,] 47.47709  10.35556 26.90996 67.09338
#> [6,] 48.67175  10.19170 28.09372 68.26898
```

## 4 Stratification for Nuisance Parameters

In many trials, parameters like the observation error variance
(\\\sigma^2\\ for continuous outcomes) or the baseline hazard function
(\\h_0(t)\\ for survival outcomes) are known to vary by site, country,
or other factors. Stratification models this known heterogeneity by
fitting these nuisance parameters separately for each level of a
grouping variable.

Use the `stratification_formula` argument to define the grouping
factor(s).

### 4.1 Example 9: Stratified Continuous Model

*Scenario*: We model outcomes where the observation error **standard
deviation** \\\sigma\\ differs by site. Stratification uses `brms`
distributional formulas to estimate separate residual variance
parameters for each site, allowing for heterogeneity in measurement
noise or outcome variability across sites.

``` r
set.seed(42)
n <- 250
sigma_by_site <- c(A = 6, B = 12, C = 18)
sample_data_strat_cont <- data.frame(
    id = 1:n,
    site = factor(sample(c("A", "B", "C"), n, replace = TRUE))
)
sample_data_strat_cont$trt <- rep(0:1, length.out = n)
sample_data_strat_cont$baseline <- rnorm(n, mean = 100, sd = 8)
sample_data_strat_cont$X1 <- factor(sample(c("A", "B"), n, replace = TRUE))

# Generate outcome with treatment heterogeneity and site-specific noise
noise <- rnorm(n, mean = 0, sd = sigma_by_site[sample_data_strat_cont$site])
sample_data_strat_cont$y <- 50 + 
                            sample_data_strat_cont$trt * (-8 + 2 * (as.numeric(sample_data_strat_cont$X1) - 1)) +
                            0.2 * (sample_data_strat_cont$baseline - 100) +
                            noise
```

``` r
# Model Fitting with Stratified Model

fit_continuous_stratified <- run_brms_analysis(
  data = sample_data_strat_cont,
  response_formula = y ~ trt,
  response_type = "continuous",
  unshrunk_terms_formula = ~ baseline,
  shrunk_predictive_formula = ~ 0 + trt:X1,
  stratification_formula = ~ site,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 2.5)",
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 1)",
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Step 1: Preparing formula and data...
#> Applying stratification: estimating sigma by 'site'.
#> Note: Marginality principle not followed - interaction term 'X1' is used without its main effect. Consider adding 'X1' to prognostic terms for proper model hierarchy.
#> 
#> Step 2: Fitting the brms model...
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 2 finished in 2.6 seconds.
#> Chain 1 finished in 3.0 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 2.8 seconds.
#> Total execution time: 3.1 seconds.
#> 
#> Analysis complete.
```

``` r
strat_continuous_summary <- summary_subgroup_effects(
  brms_fit = fit_continuous_stratified
  # All parameters automatically extracted!
)
#> --- Calculating specific subgroup effects... ---
#> Step 1: Identifying subgroups and creating counterfactuals...
#> ...detected subgroup variable(s): X1
#> Step 2: Generating posterior predictions...
#> Step 3: Calculating marginal effects...
#> Done.

print(strat_continuous_summary$estimates)
#> # A tibble: 2 × 4
#>   Subgroup Median CI_Lower CI_Upper
#>   <chr>     <dbl>    <dbl>    <dbl>
#> 1 X1: A     -6.30    -8.75    -3.83
#> 2 X1: B     -4.87    -7.13    -2.54
plot(strat_continuous_summary, title = "Stratified Continuous: Subgroup Effects")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Advanced_Functionalities_files/figure-html/unnamed-chunk-15-1.png)

### 4.2 Example 10: Stratified Survival Model

*Scenario*: We model a time-to-event outcome with stratified baseline
hazard by country using a piecewise exponential model. Stratification
via `stratification_formula = ~ country` allows the baseline hazard to
differ by country, accommodating regional differences in baseline risk
or standard of care.

``` r
set.seed(123)
n <- 250
lambda_by_country <- c(A = 0.01, B = 0.03)
surv_data_strat <- data.frame(
    id = 1:n,
    country = factor(sample(c("A", "B"), n, replace = TRUE)),
    trt = rep(0:1, length.out = n),
    age = rnorm(n, 65, 10),
    X1 = factor(sample(c("A", "B"), n, replace = TRUE))
)

# Linear predictor with treatment heterogeneity
lp <- (as.numeric(surv_data_strat$trt) - 1) * -0.6 + 
      (surv_data_strat$age - 65) * 0.03 +
      (as.numeric(surv_data_strat$trt) - 1) * (as.numeric(surv_data_strat$X1) - 1) * -0.5

# Generate event times
u <- runif(n)
lambda_vec <- lambda_by_country[surv_data_strat$country]
gamma <- 1.5 # Weibull shape parameter
true_event_time <- (-log(u) / (lambda_vec * exp(lp)))^(1/gamma)
censoring_time <- 60 # Administrative censoring time
surv_data_strat$event_status <- ifelse(true_event_time <= censoring_time, 1, 0)
surv_data_strat$event_time <- pmin(true_event_time, censoring_time)
```

``` r
# Model Fitting with Stratified Baseline Hazard
fit_surv_stratified <- run_brms_analysis(
  data = surv_data_strat,
  response_formula = Surv(event_time, event_status) ~ trt,
  response_type = "survival",
  unshrunk_terms_formula = ~ age,
  shrunk_predictive_formula = ~ 0 + trt:X1,
  stratification_formula = ~ country,
  intercept_prior = "normal(0, 5)",
  unshrunk_prior = "normal(0, 2.5)",
  shrunk_prognostic_prior = "horseshoe(scale_global = 1)",
  shrunk_predictive_prior = "horseshoe(scale_global = 1)",
  chains = 2, iter = 1000, warmup = 500, cores = 2, refresh = 0, backend = "cmdstanr"
)
#> Step 1: Preparing formula and data...
#> Response type is 'survival'. Modeling the baseline hazard explicitly using bhaz().
#> Applying stratification: estimating separate baseline hazards by 'country'.
#> Note: Marginality principle not followed - interaction term 'X1' is used without its main effect. Consider adding 'X1' to prognostic terms for proper model hierarchy.
#> 
#> Step 2: Fitting the brms model...
#> Fitting brms model...
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 finished in 5.8 seconds.
#> Chain 2 finished in 6.6 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 6.2 seconds.
#> Total execution time: 6.7 seconds.
#> Warning: 12 of 1000 (1.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
#> 
#> Analysis complete.
```

``` r
strat_surv_summary <- summary_subgroup_effects(
  brms_fit = fit_surv_stratified
  # All parameters automatically extracted!
)
#> --- Calculating specific subgroup effects... ---
#> Step 1: Identifying subgroups and creating counterfactuals...
#> ...detected subgroup variable(s): X1
#> Step 2: Generating posterior predictions...
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Step 3: Calculating marginal effects...
#> Done.

print(strat_surv_summary$estimates)
#> # A tibble: 2 × 4
#>   Subgroup Median CI_Lower CI_Upper
#>   <chr>     <dbl>    <dbl>    <dbl>
#> 1 X1: A     0.467    0.359    0.586
#> 2 X1: B     0.498    0.383    0.642
plot(strat_surv_summary, title = "Stratified Survival: Subgroup Effects (AHR)")
#> Preparing data for plotting...
#> Generating plot...
#> Done.
```

![](Advanced_Functionalities_files/figure-html/ex6-summary-plot-1.png)

## 5 Summary

This vignette demonstrated advanced functionalities in `bonsaiforest2`:

1.  **Offset variables** for count outcomes with varying exposure times
2.  **Custom prior specification** with flexible prior constraints
3.  **One-hot encoding** for full factor representation in interactions
4.  **Coefficient-specific priors** for fine-grained control
5.  **Hierarchical priors** with shared variance using `stanvars`
6.  **Stratification** for nuisance parameters that vary by groups

These features provide researchers with powerful tools for complex trial
analyses while maintaining the principled Bayesian shrinkage framework.
